# 文档 13（v2）：虚拟生命内核：植物模型（Virtual Life Kernel: Plant Model）

状态（2025‑09‑25）
- [x] v1：个体/FSM/光谱吸收与反射的概念稿
- [x] v2：接口规范化、每日步协议、能量学与形态投资、与 14/15/06/04 的交叉引用
- [ ] M2：参考实现（Genes/Plant 类与最小单元测试）
- [ ] M3：与 PopulationManager 的端到端联调（见 docs/15）
- [ ] M4：参数扫描与默认组固化（与项目 015）

实现注记（2025‑09‑25）
- 代码已提供最小 `Plant`（固定带反射、子步能量累计），当前未在主循环中被调用；小时级反照率回写通过 `EcologyAdapter`（M1）完成，详见 docs/15 与 `scripts/run_simulation.py`。

关联与交叉引用
- 项目 015（v2）：涌现生态学与光谱演化动力学（projects/015-adaptative-spectral-biogeography.md）
- 文档 14：适应性光谱物理学（docs/14-adaptive-spectroscopy.md）
- 文档 15：生态系统框架与 PopulationManager（docs/15-ecology-framework.md）
- 文档 06：能量收支框架（短波/长波与反照率耦合，docs/06-energy-framework.md）
- 文档 04：运行参数目录（本章新增 QD_ECO_* 将在 04 中补列，docs/04-runtime-config.md）

0. 目标与范围
- 定义“虚拟生命内核”的个体级模型：基因（Genes）、植物个体（Plant）、有限状态机（FSM）、内部记忆。
- 指定每日时间步的输入（生态微环境）、能量学与形态投资、输出（反射光谱/带反照率、叶面积/覆盖等）。
- 与光谱引擎（文档 14）对齐：以离散波段 I_b 与吸收 A_b 表示，避免逐纳米积分的高成本。
- 与生态框架（文档 15）对齐：由 PopulationManager 负责竞争与聚合，Plant 专注个体行为。

1. 数据结构与 API

1.1 Genes（物种/基因型定义，建议以 dataclass 表示）
- identity: str ∈ {"grass","tree",...}
- morphology:
  - root_investment_ratio: float ∈ [0,1]
  - root_depth_ratio: float ∈ [0,1]（深根 vs 匍匐，参与水分竞争权重）
  - stem_investment_ratio: float ∈ [0,1]
  - leaf_investment_ratio: float ∈ [0,1]
  - leaf_area_per_energy: float（m²/能量单位）
- spectral_absorption_curve: List[Peak]（见文档 14；多高斯峰参数 center[nm], width[nm], height）
- physiology:
  - drought_tolerance: float（0..1 的土壤湿度阈值或等效指标）
- life cycle:
  - lifespan_in_days: int
  - maturity_threshold: float（生物量或能量阈）
  - seed_energy: float（每颗种子分配的能量包）

1.2 PlantState（FSM）
- SEED, GROWING, MATURE, SENESCENT, DEAD

1.3 Plant（个体）
- genes: Genes
- age_in_days: int
- state: PlantState
- energy_storage: float
- biomass: dict{"root":float,"stem":float,"leaf":float}
- height: float（由 stem 生物量与物种参数换算）
- internal_memory:
  - accumulated_warmth: float（积温）
  - water_stress_days: int（连续水分胁迫日数）
  - last_day_energy_gain: float（诊断）
- methods（关键）：
  - update_one_day(env: EcologicalEnvironment) -> PlantDailyReport
  - update_substep(env_inst: EcologicalEnvironmentInstant, dt_seconds: float) -> SubstepReport
  - get_reflectance_bands(bands: SpectralBands) -> np.ndarray[NB]（R_b,i）
  - effective_leaf_area() -> float
  - is_alive() -> bool

1.4 EcologicalEnvironment（由 PopulationManager 构造）
- weather（当日代表值或平均）：
  - Ts, Ta, wind10, RH/q proxy, precip_daily, soil_water_index ∈ [0,1]
  - spectral_bands: I_b[NB]（地表光谱/短波带强度；见文档 14）
  - day_length_hours, seasonal_phase（可选）
- competition inputs（个体维度）：
  - light_availability ∈ [0,1]（冠层遮挡后到达该个体冠层顶的可用比例）
  - water_share ∈ [0,1]（在本格可用水的份额，依据根系权重）
- params（生态/物候/应激阈值，读取自 QD_ECO_*）

1.5 PlantDailyReport（供 PopulationManager 聚合与诊断）
- energy_gain: float（W·m⁻²·day 或统一能量单位）
- reflectance_bands: R_b[NB]（离散短波带反射率）
- leaf_area: float（m²）
- alive: bool
- transitioned_to: Optional[PlantState]
- reproduction:
  - seed_count: int
  - per_seed_energy: float

2. 每日时间步协议（Pseudo-code）

```
def update_one_day(self, env):
    # 1) 更新内在记忆
    self.internal_memory.accumulated_warmth += f_warmth(env.Ts, env.day_length_hours)
    if env.soil_water_index < self.genes.drought_tolerance:
        self.internal_memory.water_stress_days += 1
    else:
        self.internal_memory.water_stress_days = 0

    # 2) 状态转换
    self._maybe_transition(env)

    # 3) 计算当日能量收入（按光谱带）
    #    effective_light = I_b * light_availability
    I_eff = env.spectral_bands * env.light_availability  # shape [NB]
    A_b = absorbance_from_genes(self.genes.spectral_absorption_curve, env.spectral_bands_info)  # [NB], 见文档14
    E_gain = sum(I_eff * A_b * Δλ_b)  # 文档14定义带宽 Δλ_b 与单位

    # 4) 当日能量分配
    alloc = self._allocation_strategy(env, E_gain)  # root/stem/leaf/repro/storage
    self._apply_allocation(alloc)                    # 更新 biomass/height/energy_storage

    # 5) 叶片反射（带）
    R_b = 1.0 - np.clip(A_b, 0.0, 1.0)

    # 6) 衰老或死亡清理
    if self.age_in_days >= self.genes.lifespan_in_days:
        self.state = PlantState.DEAD

    # 7) 汇总报告
    return PlantDailyReport(
        energy_gain=E_gain,
        reflectance_bands=R_b,
        leaf_area=self.effective_leaf_area(),
        alive=self.is_alive(),
        transitioned_to=self.state,
        reproduction=self._reproduction_today(E_gain, alloc)
    )
```

2.1 状态转换规则（示例）
- SEED → GROWING：accumulated_warmth ≥ QD_ECO_SEED_GERMINATE_GDD；且 soil_water_index ≥ 阈值
- GROWING → MATURE：biomass_total ≥ maturity_threshold
- MATURE → SENESCENT：日长/季相或连续水分胁迫 ≥ QD_ECO_STRESS_WATER_DAYS
- 任意 → DEAD：age ≥ lifespan 或能量/生物量下限触发

2A. 时级子步协议（Sub‑daily）
用于与物理积分步（分钟级/小时级）协同，累积当日能量与胁迫，不进行形态投资与生命周期大跳转（除非设置硬触发）。典型实现：

```
def update_substep(self, env_inst, dt_seconds):
    # env_inst: EcologicalEnvironmentInstant（瞬时天气/光谱带与竞争输入）
    # 1) 即时能量增量（按光谱带）
    I_eff = env_inst.spectral_bands * env_inst.light_availability      # [NB]
    A_b = absorbance_from_genes(self.genes.spectral_absorption_curve, env_inst.spectral_bands_info)  # [NB]
    dE = np.sum(I_eff * A_b * Δλ_b) * (dt_seconds / 86400.0)           # 换算到“当日能量单位”的分数

    # 2) 缓冲累积（当日）
    self.internal_memory.last_day_energy_gain += dE

    # 3) 水分胁迫累计（小时计数）
    if env_inst.soil_water_index < self.genes.drought_tolerance:
        self.internal_memory.water_stress_days += dt_seconds / 86400.0  # 以天为单位累积
    # 注：是否采用“连续小时计数”由调用方/适配层决定

    # 4) 时级通常不做形态投资与状态转换（避免高频波动）
    # 可选硬触发：若极端胁迫/灾害，则提前置入 SENESCENT/DEAD（谨慎使用）

    # 5) 返回子步报告（可用于时级诊断）
    return SubstepReport(
        dE=dE,
        light_availability=float(np.mean(env_inst.light_availability))  # 或保留逐带
    )
```

说明
- A_b（吸收带）与叶面积/反射带 R_b 在小时级通常近似不变；因此时级地表反照率 A_b^surface 的变化主要由 I_b(t) 与云/几何驱动。
- 若需要即时将反照率回耦给短波，可由 PopulationManager 在 step_subdaily 内按缓存策略重算/复用 A_b^surface（见 docs/15 v2.1 §5.2）。

3. 能量学与形态投资

3.1 光谱带积分（与文档 14 一致）
- E_gain ≈ Σ_b I_b · A_b · Δλ_b
  - I_b：地表短波第 b 带辐照度（W·m⁻²·nm⁻¹ 等价），由 forcing/14 提供
  - A_b：本个体吸收率（由基因的吸收峰生成器在该带宽内积分/抽样）
  - Δλ_b：带宽（nm）
- 数值注意：使用带平均代替逐 nm 采样；Δλ_b 与权函数由 14 定义。

3.2 投资策略
- 默认按基因比例将“净盈余能量”投资于 root/stem/leaf（loss/maintenance 简化为常数或忽略项）：
  - invest_root = k_root · E_gain
  - invest_stem = k_stem · E_gain
  - invest_leaf = k_leaf · E_gain
- height = h0 + c_stem · biomass["stem"]^γ（γ∈[0.5,1.0]）
- leaf_area = biomass["leaf"] · leaf_area_per_energy
- reproduction（MATURE 时）：
  - E_repro = α_repro · E_gain（或由 storage 支持）
  - seed_count = floor(E_repro / genes.seed_energy)

3.3 水分胁迫与回收
- 长期水分胁迫（water_stress_days ≥ QD_ECO_STRESS_WATER_DAYS）→ 优先从叶/茎回收能量入 storage，叶面积下降，进入 SENESCENT
- SENESCENT 期停止净增长，能量用于维护/回收

4. 光谱与反照率输出（与文档 14 的接口）
- 个体反射带 R_b = 1 − A_b（限幅 [0,1]）
- PopulationManager 聚合得到地表带反照率：
  - A_b^surface = (Σ_i R_b,i · LAI_i + Albedo_soil_b · Area_soil) / Area_total
  - 注：LAI 或 leaf_area 权重由 QD_ECO_LAI_ALBEDO_WEIGHT 控制（见 15）
- GCM 短波使用带反照率 A_b^surface（文档 06），并可降维为 α_total（若需）

5. 生态微环境与竞争（与文档 15）
- light_availability 由“分层光照衰减模型”（类 Beer-Lambert）计算：按 height 从高至低，逐层衰减
- water_share 由根系生物量与深度权重在“当日可用水”上按比例分配
- Plant 不直接决定二者，均由 PopulationManager 提供

6. 环境变量（将收录至文档 04）
主开关/步长
- QD_ECO_ENABLE（默认 0）
- QD_ECO_DT_DAYS（默认 1.0）

物候与记忆
- QD_ECO_SEED_GERMINATE_GDD（默认 80）发芽积温阈值
- QD_ECO_STRESS_WATER_DAYS（默认 7）连续胁迫进入衰老阈值

投资与增长
- QD_ECO_ALLOC_ROOT/STEM/LEAF（可覆盖基因投资比例；默认由 genes 决定）
- QD_ECO_HEIGHT_EXPONENT（γ，默认 0.8）
- QD_ECO_REPRO_FRACTION（α_repro，默认 0.2）

存活与维护（可选）
- QD_ECO_MAINT_COST（默认 0；如启用则从 E_gain 抵扣）
- QD_ECO_MIN_LEAF_AREA（默认 0）小于则进入 SENESCENT/DEAD

输出与耦合
- QD_ECO_LAI_ALBEDO_WEIGHT（默认 1.0）
- QD_ECO_DIAG（默认 1）

7. 诊断与最小验收
- PlantDailyReport：E_gain、R_b、leaf_area、state 变更、seed_count
- 个体进化行为可解释：干旱时 root 增长占比更有利，叶面积减少
- 端到端（与 15 聚合）：
  - 反照率反馈后，SW_sfc 减少/增加与植被覆盖、云量变化的方向一致
  - Ecology 开启后，能量闭合（文档 06）长期无系统漂移

8. 运行示例（占位）
启用生态（每日步），仅生成个体诊断（不反馈短波）
```bash
export QD_ECO_ENABLE=1
export QD_ECO_ALBEDO_COUPLE=0
export QD_ECO_DT_DAYS=1
python3 -m scripts.run_simulation
```

启用生态并反馈带反照率（与项目 015 一致）
```bash
export QD_ECO_ENABLE=1
export QD_ECO_ALBEDO_COUPLE=1
export QD_ECO_SPECTRAL_BANDS=8
export QD_ECO_SPECTRAL_RANGE_NM=380,780
python3 -m scripts.run_simulation
```

9. 变更记录（Changelog）
- 2025‑09‑25：v2 规范化：数据结构/API、每日步协议、能量学与形态投资、光谱带接口与反照率输出、环境变量建议；与 14/15/06 交叉引用对齐
