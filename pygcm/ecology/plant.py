from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Optional, Dict

import numpy as np

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
        I_bands_weighted_scalar: float,
    ) -> PlantReport:
        """
        执行“日级慢路径”：
        - 状态机转换（基于 GDD/水分胁迫/寿命）
        - 能量投资与几何属性更新（height/leaf_area）
        - 返回最小日报（含当日能量、leaf_area、seed_count）
        I_bands_weighted_scalar: 外部带积分（Σ I_b·A_b·Δλ）的等效日能量或其代理（已按日累计）
        """
        if not self.is_alive():
            return PlantReport(energy_gain=0.0, leaf_area=self.effective_leaf_area(), state=self.state)

        transitioned = self._maybe_transition(Ts_day, day_length_hours)

        # 将子步累计的能量与外部提供的带积分代理合并（以外部为主，buffer 为补充）
        E_gain_day = max(0.0, float(I_bands_weighted_scalar)) + max(0.0, float(self._E_day_buffer))
        # 清空缓冲
        self._E_day_buffer = 0.0

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
        )
