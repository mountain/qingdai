from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Optional, Dict

import numpy as np
import os

from .genes import Genes


class PlantState(Enum):
    SEED = auto()
    GROWING = auto()
    MATURE = auto()
    SENESCENT = auto()
    DEAD = auto()


@dataclass
class PlantReport:
    energy_gain: float
    leaf_area: float
    state: PlantState
    transitioned_to: Optional[PlantState] = None
    seed_count: int = 0
    # Optional: daily reflectance bands R_b (shape [NB]) when band inputs are provided
    reflectance_bands: Optional[np.ndarray] = None


@dataclass
class Plant:
    """
    Minimal M4 individual Plant scaffold:
    - Holds genes + simple biomass bookkeeping (root/stem/leaf), height proxy, energy storage
    - Sub-daily: accumulate daily energy buffer（不做形态大跳转）
    - Daily: 状态机 & 形态学更新（能量投资、GDD、寿命/胁迫触发）
    """
    genes: Genes
    state: PlantState = PlantState.SEED
    age_days: int = 0
    # biomass "energy units"（抽象单位，与 genes.leaf_area_per_energy 配套）
    biomass: Dict[str, float] = field(default_factory=lambda: {"root": 0.0, "stem": 0.0, "leaf": 0.0})
    energy_storage: float = 0.0
    # diagnostics / memory
    gdd_accum: float = 0.0
    water_stress_days: float = 0.0
    # instantaneous geometry proxies
    height: float = 0.0
    leaf_area: float = 0.0
    # per-day energy buffer（J-equivalent proxy，按外部 I_b·A_b·Δλ·dt 累积）
    _E_day_buffer: float = 0.0

    # Parameters（可按需求扩展到 env）：
    height_exponent: float = 0.8   # height ∝ stem^γ
    repro_fraction: float = 0.2    # fraction of daily energy to reproduction when MATURE

    def effective_leaf_area(self) -> float:
        return max(0.0, float(self.leaf_area))

    def is_alive(self) -> bool:
        return self.state not in (PlantState.DEAD,)

    def update_substep(self, I_eff_scalar: float, dt_seconds: float, soil_water_index: Optional[float] = None) -> None:
        """
        Sub-daily accumulation of daily energy buffer.
        I_eff_scalar: 已按带积分后的有效光强（W m^-2）或外部提供的 J/s 等价
        dt_seconds: 子步时长
        soil_water_index: 0..1（可选），用于累积小时级水分胁迫（按天归一）
        """
        if not self.is_alive():
            return
        dE = max(0.0, float(I_eff_scalar)) * float(dt_seconds)
        self._E_day_buffer += dE
        # 水分胁迫累计（以天为单位）
        if soil_water_index is not None:
            if float(soil_water_index) < float(self.genes.drought_tolerance):
                self.water_stress_days += float(dt_seconds) / 86400.0

    def update_substep_bands(
        self,
        I_bands: np.ndarray,
        A_b_genotype: np.ndarray,
        dt_seconds: float,
        delta_lambda: Optional[np.ndarray] = None,
        light_availability: float = 1.0,
        soil_water_index: Optional[float] = None,
    ) -> None:
        """
        子步（谱带）能量累积接口（M1/M2）：
          dE = Σ_b [ I_b[b] · A_b_eff[b] · Δλ_b ] · light_availability · dt
        其中 A_b_eff = A_b_genotype · f_LAI(self.leaf_area)；f_LAI = 1−exp(-k_ext·LAI_proxy)
        说明：
          - 本个体未直接维护 LAI，使用 leaf_area 作为局地自遮蔽代理；k_ext 由 env QD_ECO_LAI_K_EXT 控制（默认 0.4）
          - I_bands 单位建议为 W m^-2 nm^-1（或带平均的等效单位），Δλ_b 为带宽（nm）
        """
        if not self.is_alive():
            return
        # Bands and absorbance
        I_b = np.asarray(I_bands, dtype=float).ravel()
        A_b = np.clip(np.asarray(A_b_genotype, dtype=float).ravel(), 0.0, 1.0)
        # Band widths
        if delta_lambda is None:
            dl = np.ones_like(I_b, dtype=float)
        else:
            dl = np.asarray(delta_lambda, dtype=float).ravel()
            if dl.shape[0] != I_b.shape[0]:
                # fallback to ones if mismatched
                dl = np.ones_like(I_b, dtype=float)
        # LAI proxy from leaf area (m^2) using simple Beer-Lambert-like factor
        try:
            k_ext = float(os.getenv("QD_ECO_LAI_K_EXT", "0.4"))
        except Exception:
            k_ext = 0.4
        LAI_proxy = max(0.0, float(self.leaf_area))
        f_LAI = 1.0 - np.exp(-k_ext * LAI_proxy)
        A_eff = np.clip(A_b * f_LAI, 0.0, 1.0)
        # Energy increment
        lv = float(max(0.0, light_availability))
        dt = float(dt_seconds)
        dE = float(np.sum(I_b * A_eff * dl)) * lv * dt
        if dE > 0.0:
            self._E_day_buffer += dE
        # Water stress accumulation
        if soil_water_index is not None:
            if float(soil_water_index) < float(self.genes.drought_tolerance):
                self.water_stress_days += dt / 86400.0

    def _maybe_transition(self, Ts_day: float, day_length_hours: float) -> Optional[PlantState]:
        """
        状态机转换（最小规则）：
        - SEED → GROWING：GDD≥阈值 & 轻微水分条件满足
        - GROWING → MATURE：叶面积或能量/生物量阈值达到（简化为 leaf_area）
        - MATURE → SENESCENT：连续水分胁迫过长或寿命接近上限
        - 任意 → DEAD：超过寿命硬阈值
        """
        transitioned = None
        # 累计 GDD（以地表温度与日长代理）
        # 这里简单：若 Ts_day>0°C，按 (Ts_day-273.15)+ 假设累积；否则 0
        gdd_today = max(0.0, float(Ts_day) - 273.15) * max(0.0, float(day_length_hours)) / 24.0
        self.gdd_accum += gdd_today

        if self.age_days >= int(self.genes.lifespan_days):
            self.state = PlantState.DEAD
            transitioned = PlantState.DEAD
            return transitioned

        if self.state == PlantState.SEED:
            if (self.gdd_accum >= float(self.genes.gdd_germinate)) and (self.water_stress_days < 1.0):
                self.state = PlantState.GROWING
                transitioned = PlantState.GROWING

        elif self.state == PlantState.GROWING:
            # 以叶面积达到某阈值作为成熟条件（简化）
            if self.leaf_area >= 0.2:  # m^2（任意阈值，可调）
                self.state = PlantState.MATURE
                transitioned = PlantState.MATURE

        elif self.state == PlantState.MATURE:
            # 若持续水分胁迫或接近寿命：进入 SENESCENT
            if (self.water_stress_days >=  float(os.getenv("QD_ECO_STRESS_WATER_DAYS", "7"))) or \
               (self.age_days >= int(0.9 * self.genes.lifespan_days)):
                self.state = PlantState.SENESCENT
                transitioned = PlantState.SENESCENT

        elif self.state == PlantState.SENESCENT:
            # 衰老阶段可在强胁迫下死亡（简化规则）
            if self.water_stress_days >= float(os.getenv("QD_ECO_STRESS_WATER_DAYS", "7")) + 5:
                self.state = PlantState.DEAD
                transitioned = PlantState.DEAD

        return transitioned

    def _apply_allocation(self, E_gain_day: float) -> None:
        """
        将“当日净能量”按基因分配到 root/stem/leaf，更新 height 与 leaf_area。
        """
        if E_gain_day <= 0.0 or not self.is_alive():
            return
        g = self.genes
        # reproduction（MATURE）优先分流一部分
        E_repro = 0.0
        if self.state == PlantState.MATURE and self.repro_fraction > 0.0:
            E_repro = self.repro_fraction * E_gain_day
        E_work = max(0.0, E_gain_day - E_repro)
        # 投资比例（已在 Genes.from_env 归一）
        self.biomass["root"] += g.alloc_root * E_work
        self.biomass["stem"] += g.alloc_stem * E_work
        self.biomass["leaf"] += g.alloc_leaf * E_work
        # height 与 leaf_area（简化）
        self.height = max(0.0, (self.biomass["stem"]) ** self.height_exponent)
        self.leaf_area = max(0.0, self.biomass["leaf"] * g.leaf_area_per_energy)
        # reproduction 暂转换为 storage 或 seed_count（在日接口返回）
        self.energy_storage += E_repro

    def update_one_day(
        self,
        Ts_day: float,
        day_length_hours: float,
        soil_water_index: float,
        I_bands_weighted_scalar: Optional[float] = None,
        *,
        I_bands: Optional[np.ndarray] = None,
        A_b_genotype: Optional[np.ndarray] = None,
        delta_lambda: Optional[np.ndarray] = None,
        light_availability: float = 1.0,
    ) -> PlantReport:
        """
        执行“日级慢路径”：
        - 状态机转换（基于 GDD/水分胁迫/寿命）
        - 能量投资与几何属性更新（height/leaf_area）
        - 返回最小日报（含当日能量、leaf_area、seed_count、可选 R_b）

        计算能量优先级：
        1) 若提供 I_bands（W m^-2 nm^-1）、A_b_genotype（[NB]）与 Δλ_b（nm），则按带积分：
              E_banded = Σ_b I_b[b] · A_eff[b] · Δλ_b · light_availability
           其中 A_eff = A_b_genotype · f_LAI(leaf_area)，f_LAI = 1−exp(−k_ext·leaf_area)，k_ext=QD_ECO_LAI_K_EXT(默认0.4)
        2) 叠加子步缓冲 self._E_day_buffer
        3) 若 I_bands_weighted_scalar 提供（外部等效带积分代理），将其与 1)+2) 相加（容错取 0）

        同时在提供 band 输入时输出当日 R_b=1−A_eff（限幅 [0,1]），便于日级诊断。
        """
        if not self.is_alive():
            return PlantReport(energy_gain=0.0, leaf_area=self.effective_leaf_area(), state=self.state)

        transitioned = self._maybe_transition(Ts_day, day_length_hours)

        # 1) 带积分路径（可选）
        E_banded = 0.0
        R_b: Optional[np.ndarray] = None
        if (I_bands is not None) and (A_b_genotype is not None):
            I_b = np.asarray(I_bands, dtype=float).ravel()
            A_b = np.clip(np.asarray(A_b_genotype, dtype=float).ravel(), 0.0, 1.0)
            if delta_lambda is None:
                dl = np.ones_like(I_b, dtype=float)
            else:
                dl = np.asarray(delta_lambda, dtype=float).ravel()
                if dl.shape[0] != I_b.shape[0]:
                    dl = np.ones_like(I_b, dtype=float)
            # f_LAI from leaf_area（LAI proxy）
            try:
                k_ext = float(os.getenv("QD_ECO_LAI_K_EXT", "0.4"))
            except Exception:
                k_ext = 0.4
            try:
                fmin = float(os.getenv("QD_ECO_LAI_FMIN", "0.05"))
            except Exception:
                fmin = 0.05
            LAI_proxy = max(0.0, float(self.leaf_area))
            f_LAI = max(fmin, 1.0 - np.exp(-k_ext * LAI_proxy))
            A_eff = np.clip(A_b * f_LAI, 0.0, 1.0)
            try:
                k_ext = float(os.getenv("QD_ECO_LAI_K_EXT", "0.4"))
            except Exception:
                k_ext = 0.4
            try:
                fmin = float(os.getenv("QD_ECO_LAI_FMIN", "0.05"))
            except Exception:
                fmin = 0.05
            LAI_proxy = max(0.0, float(self.leaf_area))
            f_LAI = max(fmin, 1.0 - np.exp(-k_ext * LAI_proxy))
            A_eff = np.clip(A_b * f_LAI, 0.0, 1.0)
            lv = float(max(0.0, light_availability))
            E_banded = float(np.sum(I_b * A_eff * dl)) * lv
            # 生成 R_b
            R_b = np.clip(1.0 - A_eff, 0.0, 1.0)

        # 2) 子步累积缓冲
        E_buffer = max(0.0, float(self._E_day_buffer))
        # 清空缓冲
        self._E_day_buffer = 0.0

        # 3) 外部等效代理（可选）
        E_proxy = 0.0 if (I_bands_weighted_scalar is None) else max(0.0, float(I_bands_weighted_scalar))

        # 当日总能量
        E_gain_day = E_banded + E_buffer + E_proxy

        # 应用形态投资
        self._apply_allocation(E_gain_day)

        # 水分胁迫按日规则：若当日水分指数良好则缓解
        if soil_water_index >= self.genes.drought_tolerance:
            self.water_stress_days = 0.0

        # 简化繁殖：MATURE 且 E_repro>0（已进 storage），折算为 seed_count（能量/常数）
        seed_energy = 1.0  # 占位常数（未来从 genes.seed_energy 读取）
        seed_count = 0
        if self.state == PlantState.MATURE and self.energy_storage > 0.0:
            seed_count = int(self.energy_storage / seed_energy)
            # 保留残余
            self.energy_storage = self.energy_storage - seed_count * seed_energy

        # 老化
        self.age_days += 1

        return PlantReport(
            energy_gain=E_gain_day,
            leaf_area=self.effective_leaf_area(),
            state=self.state,
            transitioned_to=transitioned,
            seed_count=seed_count,
            reflectance_bands=R_b,
        )
