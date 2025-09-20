# Project 009: 行星水循环闭合（E–P–R–冰/雪/大气水库）

本项目在 P006 能量收支、P007 平板海洋+海冰与 P008 湿度（q）的基础上，建立最小但自洽的**行星水循环闭合**框架：蒸发（E）→ 大气水汽与输送 → 凝结/降水（P）→ 地表汇（海/陆/冰/雪）→ 径流回海（R）→ 相变反馈。目标是实现**质量守恒与能量一致**（潜热）并提供可诊断的全球收支。

## 0. 状态（2025-09-20）
- [x] M1：水库定义与全球诊断（CWV/海冰/陆地贮水/雪）及 E–P 时空积分
- [x] M2：陆地“桶”模型与径流（R）返回海洋；降水相态分配（雨/雪）
- [ ] M3：与 P006/P007/P008 的能量–水量一致性联测（LH/LH_release、海冰相变）
- [ ] M4：长期平均守恒校核与参数标定（R 时标、雪融参数、相态阈值）
- [ ] M5：输出与可视化（E、P、R、CWV、h_ice、W_land、S_snow 的时间序列与地图）

### 实现状态更新（2025-09-20）
- 新增模块：`pygcm/hydrology.py`
  - `HydrologyParams` 与环境变量加载（`QD_RUNOFF_TAU_DAYS`、`QD_WLAND_CAP`、`QD_SNOW_THRESH`、`QD_SNOW_MELT_RATE`、`QD_RHO_W`、`QD_WATER_DIAG`）
  - `partition_precip_phase`（雨/雪相分配，按地表温度阈值）
  - `snow_step`（陆地积雪库与融雪通量）
  - `update_land_bucket`（陆地“桶”模型与线性径流；可选上限溢流）
  - `diagnose_water_closure`（全球面积加权 E–P–R 与水库质量闭合诊断）
- 集成位置：`scripts/run_simulation.py`
  - 在主循环中接入：相分配→陆雪库更新→陆桶与径流→`WaterDiag` 全局诊断输出
  - 使用湿度模块提供的 `E_flux`/`P_cond_flux` 保持与 P006 潜热（`LH`/`LH_release`）能量一致性
  - 径流回海以“记账”闭合，不引入复杂海洋水库
- 快速验证：短程运行可见 `WaterDiag` 打印；风场/湿度自旋后 ⟨E⟩、⟨P⟩、⟨R⟩ 将逐步非零

关联项目（交叉引用）：
- P006 能量收支：潜热通量 LH（地表）与 LH_release（大气）能量一致；TOA/SFC/ATM 诊断统一
- P007 平板海洋+海冰：海冰厚度/掩膜的相变质量与能量一致性；海/陆/冰面蒸发差异
- P008 湿度（q）：E 作为 q 源，超饱和凝结为 P，二者决定大气水库（Column Water Vapor, CWV）变化

## 1. 目标
- 定义并演化最小水库集合：大气水汽（CWV）、海冰（h_ice）、陆地贮水（W_land，可含土壤/地表水总和）、积雪（S_snow，可选）。
- 为陆地降水建立“桶”模型与**径流 R**返海闭合，支持多时标排水。
- 在全局与长期平均下实现**水量守恒**；与 P006 的能量守恒（潜热）一致。
- 输出标准化诊断，支撑参数标定与场景对比。

## 2. 水库与连续方程（最小闭合）

记网格面积权重为 w = cos(lat)。对每个时间步 dt，有：

### 2.1 大气水汽（CWV）
- 近地层单层近似下，CWV ~ ρ_a h_mbl q（或直接以网格平均 q 代表趋势）
- 整体收支（面积权重 w）：
  d⟨CWV⟩/dt ≈ ⟨E⟩ − ⟨P⟩ − ⟨∇·(flux_q)⟩
  - 在全星积分上，平流散度项的面积积分→0（周期边界与无通量极边），故有
  - d⟨CWV⟩/dt ≈ ⟨E⟩ − ⟨P⟩

### 2.2 海冰（h_ice；P007）
- 冻结/融化以能量为主导（见 P007 §2.3），对应质量变化
- 海冰质量变化对水量的贡献：
  d⟨M_ice⟩/dt = ρ_i ⟨dh_ice/dt⟩
  - 冻结（dh_ice>0）：由海洋转入冰库
  - 融化（dh_ice<0）：由冰库回海洋

### 2.3 陆地贮水（W_land，桶模型）
- 定义陆地网格的综合贮水量（mm 水层或 kg/m²），受降水与蒸发与径流控制：
  dW_land/dt = P_land − E_land − R
  - R = W_land / τ_runoff（最简线性排水）或分解快/慢径流
  - 可选：上限容量 W_cap，超过部分按“快流”立刻外排

### 2.4 积雪（S_snow，可选）
- 简化相态分配：当 T_s < T_snow_thresh 时，P→P_snow 累积为 S_snow；T_s≥阈值时按融雪速率 M_snow 融化并进入 W_land 或直接径流 R
- dS_snow/dt = P_snow − M_snow

### 2.5 全局水量守恒（期望）
- 全星积分（面积权重 w）随时间（长期平均）应满足：
  d/dt [⟨CWV⟩ + ⟨M_ice/ρ_w⟩ + ⟨W_land⟩ + ⟨S_snow⟩] ≈ ⟨E⟩ − ⟨P⟩ − ⟨R⟩
- 稳态长期平均：LHS→0，因而 ⟨E⟩ ≈ ⟨P⟩ + ⟨R⟩

## 3. 与 P006/P007/P008 的一致性

- 能量一致性（P006）：
  - 地表：LH = L_v·E 从地表能量中扣除；SFC_net 包含 −LH
  - 大气：LH_release = L_v·P_cond 加热大气；ATM_net 包含 +LH_release
  - 期望长期平均：⟨LH⟩ ≈ ⟨LH_release⟩
- 海冰相变（P007）：
  - 由净能量确定冻结/融化，对应质量通量与反照率切换；冰覆盖处 E 显著减弱
- 湿度 q（P008）：
  - E 进入 q 源项；超饱和→P_cond 为 q 汇；与降水诊断融合（双通道：动力散度与超饱和）

## 4. 设计与接口

### 4.1 新增/扩展模块（建议）
- 新建 `pygcm/hydrology.py`（或集成至 humidity.py）：
  - `update_land_bucket(W_land, P_land, E_land, dt, params) -> (W_next, R)`
  - `partition_precip_phase(P, T_s, T_thresh) -> (P_rain, P_snow)`
  - `snow_step(S_snow, P_snow, melt_rate, dt) -> (S_next, melt_flux)`
  - `diagnose_water_closure(q, E, P, R, h_ice, W_land, S_snow, grid) -> dict`（全局/分区收支）
- 扩展输出记录器：累计 ⟨E⟩、⟨P⟩、⟨R⟩、各水库变化、闭合误差

### 4.2 主循环接入次序（建议）
1) 云/降水诊断（P003/P005；与 P008 的 P_cond 融合）
2) 反照率与辐射分量（P006），动力步
3) 湿度步（P008）：E、LH；P_cond、LH_release；更新 q
4) 地表/大气能量步（P006）：使用 LH/LH_release；海冰相变（P007）
5) 水文步（本项目）：P 的相态分配（雨/雪）；更新 W_land/S_snow，计算 R 并将其加回“海洋水库”（记账）
6) 诊断输出：TOA/SFC/ATM 与 E–P–R–水库变化与闭合误差

## 5. 任务拆解
- T1（M1）诊断与记录：实现 `diagnose_water_closure`，打印/保存 ⟨E⟩、⟨P⟩、⟨R⟩、d⟨CWV⟩/dt、d⟨M_ice⟩/dt、d⟨W_land⟩/dt、d⟨S_snow⟩/dt 与闭合误差
- T2（M2）陆地桶+径流：`update_land_bucket` 与线性/双时标排水；阈值容量可选
- T3（M3）相态与雪：按 T_s 阈值分配 P_rain/P_snow；融雪入 W_land 或直接 R
- T4（M4）一致性联测：与 P006/P007/P008 联动运行，校核长期平均守恒与能量一致
- T5（M5）标定与输出：扫描 τ_runoff、T_thresh、melt_rate 等，固化稳定默认；完善图像与 NetCDF 输出

## 6. 参数与环境变量（建议默认）
- 径流与桶：
  - QD_RUNOFF_TAU_DAYS（默认 10）：线性径流时标（天）
  - QD_WLAND_CAP（默认 None/不启用）：桶容量（mm）
- 相态与雪：
  - QD_SNOW_THRESH（默认 273.15 K）：雨/雪阈值
  - QD_SNOW_MELT_RATE（默认 5 mm/day 等效）：温和融雪率（可温度调制，首版常数）
- 诊断：
  - QD_WATER_DIAG（默认 1）：打印/保存闭合核算
- 与 P006/P007/P008 对齐的参数：L_v、h_mbl/E2Q 因子、α_ice、T_freeze 等使用对方模块参数

## 7. 诊断与验收标准
- 守恒性：
  - [ ] ⟨LH⟩ ≈ ⟨LH_release⟩（能量潜热一致，容差阈值可配）
  - [ ] 长期平均 d/dt 水库总量 ≈ 0，且 ⟨E⟩ ≈ ⟨P⟩ + ⟨R⟩（水量守恒）
- 物理合理性：
  - [ ] 海陆差异：海上 E 大、陆地 E 小、冰上 E ≪；降水与地形/动力有一致结构
  - [ ] 冰/雪：冷区积雪/结冰季节性合理；夏季融化增强 R
- 输出：
  - [ ] 提供 E、P、R、CWV、h_ice、W_land、S_snow 的地图与时间序列
  - [ ] 日志中打印全局闭合与关键参数

## 8. 与其他项目的交叉引用
- P006（能量）—— 使用其 LH/LH_release；水量闭合诊断与 TOA/SFC/ATM 统一输出
- P007（海冰/海洋）—— 冰相变质量与能量一致；海/陆/冰蒸发权重；径流入海闭合
- P008（湿度）—— E/P 决定 q 的源汇；CWV 诊断直接来源于 q

## 9. 运行示例
```bash
# 启用能量/湿度/水文诊断
export QD_ENERGY_W=1
export QD_ENERGY_DIAG=1
export QD_HUMIDITY_DIAG=1
export QD_WATER_DIAG=1

# 外部地形与反照率底图（可选）
export QD_TOPO_NC=$(ls -t data/*.nc | head -n1)
export QD_USE_TOPO_ALBEDO=1

# 水文参数（示例，温和默认）
export QD_RUNOFF_TAU_DAYS=10      # 陆地径流时标（天）
export QD_SNOW_THRESH=273.15      # 雨/雪阈值（K）
export QD_SNOW_MELT_RATE=5        # 融雪速率（mm/day）
# export QD_WLAND_CAP=50           # 陆面“桶”容量（mm，可选，不设则禁用上限）

# 可选地形降水增强
# export QD_OROG=1
# export QD_OROG_K=7e-4

python3 -m scripts.run_simulation
```

---

注：P009 的实现重点在于“**闭合**与**一致**”。在不引入复杂陆面/河网模型的前提下，用最小参数集实现守恒与可解释的行星水循环，并与 P006–P008 共同构成青黛 GCM 的一体化“能量—水量”框架。
