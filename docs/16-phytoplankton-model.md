# **文档 16（草案）：虚拟生命内核：浮游生物模型 (Virtual Life Kernel: Phytoplankton Model)**

本文档基于文档 13（植物模型）和文档 14（适应性光谱物理学），定义了一个适用于水生环境中漂浮、自养生命（如蓝藻、硅藻等）的个体级模型。

**0. 目标与范围**

  - 在现有 `Plant` 模型基础上，定义一个 `Phyto` 模型，代表单个或菌落浮游生物。
  - 模型需包含其生命周期、漂浮/沉降特性、以及独特的光谱学特征。
  - 与生态框架（文档 15）对接，由 `PopulationManager` 负责其在水体中的平流（漂流）、竞争与聚合效应。

**1. 数据结构与 API 调整**

核心思想是重用大部分 `Plant` 模型的数据结构，但对特定字段进行重新诠释或替换，以反映浮游生物的特性。

**1.1 `Genes` (物种/基因型定义)**

浮游生物的基因将不再关注根、茎、叶的投资，而是关注细胞结构、浮力、营养盐利用效率和光合色素。

  - `identity: str ∈ {"cyanobacteria", "diatom", ...}`
  - `morphology`:
      - **`cell_wall_investment_ratio: float ∈ [0,1]`**: 用于防御或结构的生物量投资比例。
      - **`buoyancy_factor: float`**: 决定在水中的相对浮力。正值上浮，负值下沉。可用于模拟蓝藻的气囊调节。
      - **`surface_area_per_biomass: float`**: m²/生物量单位。替代 `leaf_area_per_energy`，描述其光合作用和营养吸收的有效表面积。
  - `spectral_absorption_curve: List[Peak]`: **核心特性，保持不变**。蓝藻可以定义在红光和蓝绿光区域有吸收峰（对应叶绿素 a 和藻胆素）。
  - `physiology`:
      - **`nutrient_affinity: dict{"N": float, "P": float, "Fe": float, ...}`**: 对不同营养盐（氮、磷、铁等）的吸收能力或 Michaelis-Menten 参数。
      - **`nutrient_stress_threshold: float`**: 触发休眠或衰老的营养盐浓度阈值。
  - `life_cycle`:
      - `lifespan_in_days: int`
      - `maturity_threshold: float` (生物量阈值)
      - **`division_energy: float`**: 替代 `seed_energy`，指细胞分裂（无性繁殖）所需的能量包。

**1.2 `PhytoState` (FSM)**

有限状态机保持不变：`ACTIVE` (生长与分裂), `DORMANT` (孢子/休眠), `DEAD`。这比陆生植物的 `SEED` -\> `GROWING` -\> `MATURE` 更简单直接。

**1.3 `Phyto` (个体)**

  - `genes: Genes`
  - `age_in_days: int`
  - `state: PhytoState`
  - `energy_storage: float`
  - **`biomass: float`**: 简化为单一浮点数，不再区分 root/stem/leaf。
  - **`depth: float`**: 在水体中的深度（米），替代 `height`。由 `PopulationManager` 根据浮力因子、水流和湍流混合进行更新。
  - `internal_memory`:
      - `accumulated_warmth: float`
      - **`nutrient_stress_days: int`**: 替代 `water_stress_days`。
  - `methods` 接口保持一致，但实现细节不同。

**1.4 `AquaticEnvironment` (由 `PopulationManager` 构造)**

这是对 `EcologicalEnvironment` 的水生版本扩展。

  - `hydrodynamics` (水动力环境):
      - `water_temperature: float`
      - `current_velocity: tuple(u, v, w)`: 水流速度（用于平流漂移）。
      - `turbulence_mixing_coeff: float`: 湍流混合系数，影响垂直位置。
  - `biogeochemistry` (生物地球化学):
      - `nutrient_concentration: dict{"N": float, "P": float, ...}`: 环境营养盐浓度。
      - `spectral_bands_at_depth: I_b[NB]`: **关键输入**。到达该个体深度的光谱带强度，已经过水体和上层生物的衰减。
  - `competition_inputs`:
      - `light_availability` 隐含在 `spectral_bands_at_depth` 中，由 `PopulationManager` 根据所有 `Phyto` 的深度和密度计算光衰减后得出。
      - `nutrient_share` 由 `PopulationManager` 根据个体周围的营养浓度和吸收动力学计算。

**2. 每日时间步协议 (调整后的 `update_one_day`)**

伪代码的核心逻辑不变，但环境变量和内部计算发生变化。

```python
def update_one_day(self, env: AquaticEnvironment):
    # 1) 更新内在记忆 (积温、营养胁迫)
    # ... 逻辑类似，但使用水温和营养浓度
    if env.nutrient_concentration["N"] < self.genes.nutrient_stress_threshold:
        self.internal_memory.nutrient_stress_days += 1
    else:
        self.internal_memory.nutrient_stress_days = 0

    # 2) 状态转换 (ACTIVE <-> DORMANT <-> DEAD)
    # ...

    # 3) 计算当日能量收入（核心光谱学部分不变）
    #    I_eff_at_depth = env.spectral_bands_at_depth
    A_b = absorbance_from_genes(self.genes.spectral_absorption_curve, ...) #
    E_gain = sum(env.spectral_bands_at_depth * A_b * Δλ_b) #

    # 4) 当日能量分配 (简化)
    #    主要用于 biomass 增长、维持和分裂繁殖 (reproduction)
    # ...

    # 5) 生成报告 (同 PlantDailyReport)
    #    报告自身的反射光谱带 R_b = 1.0 - A_b
    #    报告生物量 biomass 和存活状态
    #    报告分裂产生的 "子细胞" 数量和能量
    # ...
```

**3. 光谱学与水体反照率**

这里完全复用**文档 14**的框架，但应用的场景从“地表”变成了“水体表面”。

  - **个体光谱**: `Phyto` 个体根据其基因定义的 `spectral_absorption_curve` 计算出吸收带 `A_b` 和反射带 `R_b`。
  - **水体反照率聚合**: `PopulationManager` 需要聚合所有 `Phyto` 个体的光学特性，并结合水体本身的光学特性，来计算整个水域的表面反照率 `A_b^surface`。
      - `A_b^surface` 主要由水体自身的反射率（与水的纯净度、风速等有关）和水体中浮游生物的后向散射（backscattering）共同决定。
      - 一个简化的聚合公式可以是：
        `A_b^surface = Albedo_water_b + Σ_i (Backscatter_b,i · Biomass_i · concentration_at_surface_i)`
          - `Backscatter_b,i` 是个体 `i` 在波段 `b` 的后向散射率，可以近似为 `R_b,i` 的一个函数。
          - 该反照率将反馈给能量收支模块（文档 06），影响整个行星的能量平衡。

**4. `PopulationManager` 的新职责**

`PopulationManager` 需要进行较大扩展，以管理水生环境。

  - **平流输运 (Advection)**: 根据 `AquaticEnvironment` 提供的水流速度，在每个时间步移动网格内的 `Phyto` 种群。这就是“漂流”的实现。
  - **垂直动态 (Vertical Dynamics)**: 根据每个 `Phyto` 的 `buoyancy_factor` 和环境的湍流混合，更新其 `depth`。
  - **光场模拟 (Light Field Simulation)**: 基于 Beer-Lambert 定律，根据水中所有 `Phyto` 的密度和深度分布，计算光在水下的衰减，为每个个体提供准确的 `spectral_bands_at_depth`。
  - **养分循环 (Nutrient Cycling)**: 管理和更新每个网格的 `nutrient_concentration`，处理消耗和再循环。

**5. 蓝藻基因示例**

```python
# 示例：一种适应浑浊水体、有藻胆素的蓝藻基因
cyanobacteria_genes = Genes(
    identity="cyanobacteria_A",
    morphology={
        "buoyancy_factor": 0.1,  # 略微正浮力
        "surface_area_per_biomass": 150.0
    },
    spectral_absorption_curve=[
        {"center": 440, "width": 40, "height": 0.8}, # 叶绿素 a 蓝光区
        {"center": 620, "width": 30, "height": 0.6}, # 藻蓝蛋白
        {"center": 680, "width": 25, "height": 0.9}  # 叶绿素 a 红光区
    ],
    physiology={
        "nutrient_affinity": {"N": 0.8, "P": 0.6},
        "nutrient_stress_threshold": 0.1
    },
    life_cycle={
        "lifespan_in_days": 10,
        "maturity_threshold": 100.0,
        "division_energy": 20.0
    }
)
```

这个 `Phyto` 模型的设计，充分利用了您现有框架的模块化和通用性，特别是强大的光谱学引擎。通过对数据结构和环境输入的针对性调整，可以构建一个功能丰富且与气候模型紧密耦合的虚拟浮游生态系统。