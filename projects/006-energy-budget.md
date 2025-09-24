# Project 006: 物理引擎升级——从牛顿冷却到动态能量收支框架

本文档将把现有的热力学“牛顿冷却”方案升级为“显式的动态能量收支”框架，并给出清晰的接口与落地路径，以便后续可扩展到更复杂过程（海冰、陆面、水循环等）。本版文档遵循项目其他文档的一致结构：状态、目标、设计与接口、任务拆解、约束与一致性、测试与验收、运行示例、后续扩展，并在末尾给出我的修订建议。

## 状态（2025-09-20）
- [x] 里程碑 1：辐射核心就绪（短波/长波模块 + 地表能量收支，边界层通量暂置 0）
- [x] 里程碑 2：引入感热通量（SH）并稳定耦合
- [x] 里程碑 3：闭合潜热环（LH，基于鲍文比/经验代理，与降水潜热释放一致）
- [ ] 里程碑 4：能量闭合诊断（TOA/SFC/大气收支）与长期守恒验证
- [ ] 里程碑 5：参数标定与默认参数组（与 P005 地形/反照率、云方案一致）

### 实现状态更新（2025-09-20）
- 已完成：M1 + M2 + M3
  - M1：辐射核心（短波/长波）+ 地表能量方程（SH/LH=0）接入主循环。
  - M2：引入感热通量（SH），按海陆差异的 Bowen 比方案，稳定耦合到地表能量收支。
  - M3：闭合潜热环（LH 与 LH_release），其中地表潜热 LH = L_v·E 由湿度模块（P008）提供，大气潜热释放 LH_release = L_v·P_cond 注入大气能量路径。
- 主要变更：
  - 新增：`pygcm/energy.py`
    - `shortwave_radiation` / `longwave_radiation`
    - `boundary_layer_fluxes`（SH/LH 的块体公式与 Bowen 比）
    - `compute_atmos_height_tendency` 与 `integrate_atmos_energy_height`（将 SW_atm/LW_atm/SH/LH_release 转为大气位势高度 h 的增温通量）
  - 修改：
    - `pygcm/dynamics.py`：在 `time_step` 中使用上述能量模块；调用 `boundary_layer_fluxes` 计算 SH；结合湿度模块输出的 `LH` 与 `LH_release`，通过 `integrate_atmos_energy_height` 将（SW_atm + LW_atm + SH + LH_release）注入大气能量路径（更新 h）。
    - `scripts/run_simulation.py`：传入 `albedo`，并维持能量/湿度诊断打印。
- 运行开关与参数（补充）：
  - `QD_ENERGY_W`（0..1，默认 0）：混合权重；=1 完全启用能量收支耦合（含大气能量注入）。
  - `QD_ENERGY_DIAG`（默认 1）：每≈200 步打印 TOA/SFC/ATM 收支诊断。
  - `QD_T_FLOOR`（默认 150 K）：夜侧温度下限；`QD_SW_A0/QD_SW_KC`、`QD_LW_EPS0/QD_LW_KC`、`QD_CS` 可调。
  - `QD_CH`、`QD_BOWEN_LAND`、`QD_BOWEN_OCEAN`：边界层/鲍文比参数（M2）。
  - `QD_ATM_H`、`QD_RHO_A`：大气单层有效厚度与密度（用于将能量通量转为 h 的倾向，M3）。
- 快速验证（短时、不出图）：
  ```
  # 单位：行星日
  QD_ENERGY_W=1 QD_ENERGY_DIAG=1 QD_SIM_DAYS=0.05 QD_PLOT_EVERY_DAYS=100 python3 -m scripts.run_simulation
  ```
- 后续里程碑：M4 能量闭合诊断（TOA/SFC/ATM）与长期守恒验证；M5 参数标定与默认参数组。

### 实现状态更新（2025-09-21）
- 新增：温室系数自动微调（Autotune）与持久参数
  - 在 `pygcm/energy.py` 增加 `autotune_greenhouse_params(params, diag)`，基于全局能量诊断 `TOA_net` 对长波系数作温和校正：
    - 若 `TOA_net > 0`（行星增能），略减小 `lw_eps0/lw_kc` 以提升 OLR；
    - 若 `TOA_net < 0`（行星失能），略增大 `lw_eps0/lw_kc` 以降低 OLR。
  - 步长（可配）：`QD_TUNE_RATE_EPS`（默认 5e-5）、`QD_TUNE_RATE_KC`（默认 2e-5）；范围约束：`lw_eps0∈[0.30,0.98]`、`lw_kc∈[0.0,0.80]`。
  - 在 `scripts/run_simulation.py` 中按 `QD_ENERGY_TUNE_EVERY`（默认 50 步）调用，开关 `QD_ENERGY_AUTOTUNE=1` 开启；能量参数 `eparams` 启动时加载并在运行期持续持有与微调。
- 与 P011（海洋）耦合更新（见 P011 文档）：
  - 海–气热交换：`Q_net = SW_sfc − LW_sfc − SH − LH` 默认作为海表垂直通量（`QD_OCEAN_USE_QNET=1`）进入海洋混合层能量；
  - 海–气动量交换：风应力使用相对风 `|V_a − U_o| (V_a − U_o)` 形式，避免对随风同向的强海流过度加速。
- 可视化/诊断：
  - 状态图新增“行星日累计降水”（mm/day）；温度面板（Ts/Ta/SST）统一色表与刻度；风/洋流流线统一速度色标；
  - 维持 `compute_energy_diagnostics` 打印 TOA/SFC/ATM 收支以观察向平衡收敛趋势。

## 1. 背景与目标

- 背景：当前 GCM 使用牛顿冷却（T 向 Teq 弛豫）近似，虽稳定高效，但将辐射-地表-云-大气的能量交换黑箱化，难以扩展与诊断。
- 目标：以显式能量通量为核心，分别为“大气单层（Ta）”与“地表单层（Ts）”建立能量平衡方程，并通过短波（SW）、长波（LW）、边界层通量（SH/LH）联结动力与热力过程，形成可诊断、可扩展、可守恒的热力框架。

## 2. 架构设计与数据流

### 2.1 状态与通量
- 状态变量：
  - 大气温度 Ta（可由厚度场 h 代理或单独诊断）
  - 地表温度 Ts（已有）
  - 云量 C（已有）
- 通量分量：
  - 短波：SW_atm（大气吸收）、SW_sfc（地表吸收）
  - 长波：LW_atm（作用于大气的净长波）、LW_sfc（作用于地表的净长波）
  - 边界层通量：SH（感热）、LH（潜热）
  - 降水潜热释放：LH_release（回馈大气）

### 2.2 主循环数据流（替代牛顿冷却位置）
1) 天文强迫：计算两星入射 I = ISR_A + ISR_B  
2) 诊断：由动力场求降水 P；由 P 与代理项更新云量 C（沿用 P003/P005 流程）  
3) 短波：给定 I 与总反照率 α，分配 SW_atm、SW_sfc  
4) 长波：灰体单层大气近似，计算 LW_atm、LW_sfc（含云/温室增强）  
5) 边界层：计算 SH、LH（块体公式 + 鲍文比/Bowen Ratio 方案）  
6) 地表能量收支：C_s dTs/dt = SW_sfc − LW_sfc − SH − LH  
7) 大气能量收支：c_p ρ H dTa/dt = SW_atm + LW_atm + SH_from_sfc + LH_release  
8) 动力步：对 Ta（或与 h 的联系）进行平流与耦合，进入下一个时间步

## 3. 模块与接口（建议新文件：pygcm/energy.py）

为保证清晰与可测试性，建议将能量过程模块化，提供纯函数接口。

### 3.1 短波辐射 shortwave_radiation(I, albedo, cloud, params)
输入：
- I：层顶入射（W/m²），来自两星几何（已实现）
- albedo：总反照率场（由 base_albedo_map、冰-反照率、云反照率融合）
- cloud：云量 C
- params：短波参数（如大气短波吸收系数 A_sw0、云短波吸收增益 k_sw_cloud）

输出：
- SW_atm, SW_sfc（W/m²），满足 I = SW_atm + 反射 + SW_sfc

示例参数化（简化占位，待标定）：
- 反射：R = I · α_total
- 大气吸收：SW_atm = I · A_sw0 · (1 + k_sw_cloud · C)
- 地表吸收：SW_sfc = I − R − SW_atm

### 3.2 长波辐射 longwave_radiation(Ts, Ta, cloud, params)
输入：
- Ts、Ta：地表/大气温度（K）
- cloud：云量 C
- params：大气发射率 ε_a0、云增强 k_lw_cloud、夜侧温度下限 T_floor（可选）

输出：
- LW_atm、LW_sfc（W/m²）

灰体单层大气近似（示意）：
- ε_eff = clip(ε_a0 + k_lw_cloud · C, 0, 1)
- 上向到太空：OLR = ε_eff σ Ta⁴ + (1 − ε_eff) σ Ts⁴
- 向下长波：DLR = ε_eff σ Ta⁴
- 地表净长波：LW_sfc = DLR − σ Ts⁴
- 大气净长波：LW_atm = （地表向上被大气吸收部分 + 大气向下/向上发射之净效应）  
  注：实现时按能量守恒构造 LW_atm，让（SW/LW/SH/LH）全局守恒。

### 3.3 边界层通量 boundary_layer_fluxes(Ts, Ta, u, v, land_mask, params)
输入：
- Ts、Ta、u/v、land_mask
- params：C_H（感热交换系数）、B（鲍文比）、ρ、c_p、风速权重等

输出：
- SH（W/m²）、LH（W/m²）
- SH = ρ c_p C_H |V| (Ts − Ta)
- LH = SH / B（B>0），对 B→0 做限幅；B 可因地类/干湿度代理/纬度而变

### 3.4 地表能量收支 integrate_surface_energy(Ts, SW_sfc, LW_sfc, SH, LH, dt, params)
- C_s dTs/dt = SW_sfc − LW_sfc − SH − LH
- 对 Ts 进行时间推进。对 stiff 项（LW）可采用半隐式/指数积分以提升稳定性。

### 3.5 大气能量收支 integrate_atmos_energy(Ta, SW_atm, LW_atm, SH_from_sfc, LH_release, dt, params)
- c_p ρ H dTa/dt = SW_atm + LW_atm + SH_from_sfc + LH_release
- 与动力核（h、u、v）的关系：当前原型中 Ta 与 h 存在简化关系；可先保持弱耦合（作为诊断/附加温度），逐步收敛到一致框架（见修订建议）。

## 4. 参数与默认值（建议，可通过环境变量）
- 短波：
  - QD_SW_A0（默认 0.06）大气短波基吸收
  - QD_SW_KC（默认 0.20）云短波吸收增益
- 长波：
  - QD_LW_EPS0（默认 0.70）无云大气发射率
  - QD_LW_KC（默认 0.20）云长波增强
  - QD_T_FLOOR（默认 150 K）夜侧温度下限（防止非物理冷塌）
- 边界层：
  - QD_CH（默认 1.5e-3）感热交换系数
  - QD_BOWEN_LAND（默认 0.7）、QD_BOWEN_OCEAN（默认 0.3）
- 地表热容量：
  - QD_CS（默认 2e7 J m⁻² K⁻¹）
- 诊断开关：
  - QD_ENERGY_DIAG（默认 1）输出 TOA/SFC/大气收支与全球积分

## 5. 任务分解（实施路线）
- T1（M1）辐射核心
  - 实现 shortwave_radiation/longwave_radiation，返回 SW_atm/SW_sfc/LW_atm/LW_sfc
  - Ts 能量方程集成（SH=LH=0）
  - 运行行星日尺度，验证昼夜/季节节律
- T2（M2）感热通量
  - 引入 SH（块体公式），风速用 |V|=sqrt(u²+v²)，C_H 常数或随地形/粗糙度调节
  - 稳定性检查与参数整定
- T3（M3）潜热闭环
  - 引入 LH 与 LH_release（与降水过程能量闭合），Bowen 比方案分海陆
  - 与云/降水的一致性（可由 P、C 调节 ε_eff 或 B）
- T4（M4）能量闭合诊断
  - TOA：入射 − 反射 − OLR；SFC：净 SW/LW − SH − LH
  - 全球空间积分 + 时间平均应近守恒；输出 NetCDF/日志
- T5（M5）参数标定与默认组
  - 参数扫描（CH、BOWEN、EPS0、KC 等），确定默认稳定解
  - 与 P005（地形/反照率/云）联测，输出对比图与诊断报表

## 6. 接口与一致性约束
- 能量守恒：TOA 与 SFC 的全球净通量长期平均应接近 0（容差设阈）
- 物理一致性：
  - 云增多 → 短波表面吸收下降、长波温室增强
  - 陆地相对较大的 B（更干燥）→ 感热占比更高
  - 极区与高山（通过 base_albedo、friction、C_H）应体现差异
- 数值稳定性：
  - 允许对 LW、BL 通量使用半隐式或限幅
  - 对 Ts 引入混合/平滑的弱项（已在动力里采用半拉氏/扩散思想）

## 7. 测试与验收标准
- 功能性：
  - [ ] 各模块独立单元测试（输入/输出维度与合理范围）
  - [ ] 集成后可运行 ≥ 若干行星日，出图与日志齐全
- 守恒性：
  - [ ] TOA/SFC/大气能量收支长期平均接近 0（阈值可配，如 |净| < 2 W/m²）
- 物理合理性：
  - [ ] 云增多导致地表冷却/大气增温的趋势
  - [ ] 海陆差异在 SH/LH 中体现；极区夜侧不出现非物理冷塌
- 输出：
  - [ ] 日志中打印参数与能量诊断
  - [ ] 输出 NetCDF/图像（可与 P003/P005 现有图并列）

## 8. 运行示例（建议占位）
```bash
# 开启能量诊断与夜侧温度下限
export QD_ENERGY_DIAG=1
export QD_ENERGY_W=1
export QD_T_FLOOR=160

# 辐射与边界层参数试验
export QD_SW_A0=0.06
export QD_SW_KC=0.20
export QD_LW_EPS0=0.70
export QD_LW_KC=0.25
export QD_CH=1.5e-3
export QD_BOWEN_LAND=0.7
export QD_BOWEN_OCEAN=0.3

# 与 P005 外部地形/反照率耦合
export QD_TOPO_NC=$(ls -t data/*.nc | head -n1)
export QD_USE_TOPO_ALBEDO=1

python3 -m scripts.run_simulation
```

## 9. 后续扩展
- 可选层次提升：Ta 与 h 的一体化（将 Ta 与厚度/位势联系得更物理一致）
- 云光学厚度诊断：由 C 推断 τ_cloud，进入 SW/LW 的透过/发射率
- 地表类型：引入简单土壤湿度/海冰，动态调节 B、α、ε
- 输出与试验管理：与 P006b（数据管线）协作，完善重启/元数据记录

---

## 10. 我的修订建议（在原方案基础上的改进）

1) 分阶段替换牛顿冷却（可渐进启用）  
- 引入“混合权重” w ∈ [0,1]，以 Teq 式与能量收支式加权，为迁移期提供稳定器：  
  T_new = (1−w)·T_newton + w·T_energy（w 随时间或试验切换）  
- 在回归测试稳定后，将 w→1，完全采用能量收支框架。

2) 夜侧温度下限与隐式项  
- 在 LW 与 SFC 集成中引入 T_floor，防止极端夜侧/高原冷塌（尤其双星几何下夜侧长时间无日照）。  
- 对地表能量方程的 σ Ts⁴ 项采用半隐式更新或指数积分，提升时间步长下的稳定性。

3) 云-辐射一致性与可诊断性  
- 用 ε_eff = ε0 + k_lw_cloud·C 与 A_sw_eff = A0 + k_sw_cloud·C 保持“云增多→短波吸收↑（大气）/地表吸收↓、长波温室↑”的一致趋势；  
- 输出每步全局积分：TOA（I − R − OLR）、SFC（SW_sfc − LW_sfc − SH − LH）、ATM（差值），并累计时间平均，作为守恒标尺。

4) 边界层通量更贴近地类  
- C_H 与 B 建议随地表/风速/海陆/纬度调节：  
  C_H = C_H0 · (1 + k_rough · friction_map_norm)，B = B_land/ocean ± f(|lat|, P, Ts)。  
- 降水多/湿度代理（由 P 与 C）→ B 较小，潜热占比上升；干旱内陆 → B 较大。

5) 与现有动力核的一致性路径  
- 短期：Ta 作为诊断温度跟随能量收支与弱耦合平流；  
- 中期：逐渐将 Ta 与 h 的关系在动力核中显式化（例如利用静力近似联系 geopotential 与温度），减少“两个温度”的割裂。

6) 参数命名与默认值（便于实验管理）  
- 统一以 QD_SW_*、QD_LW_*、QD_BL_*（边界层）、QD_CS、QD_ENERGY_* 命名，默认参数“温和，不爆裂”，并在 README 中给出一组稳定默认。

7) 与 P005 的互操作  
- 使用 base_albedo_map 与云量 C 统一进入 α_total；  
- friction_map 可间接影响 C_H（粗糙度代理），并对风速/地形信号回馈 BL 通量。

以上修订旨在在不牺牲当前可运行性的前提下，逐步将模型从经验弛豫推向物理可诊断、能量可守恒的框架，并留足扩展空间（海冰/土壤湿度/云光学厚度）。

## 11. 与 P007 / P008 / P009 的一体化方案（交叉引用）

- P007（平板海洋/海冰, Slab Ocean + Sea Ice）
  - 提供 C_s_map（§2.1）以替代 P006 的标量 C_s，显式体现海/陆/冰的热惯性差异；
  - 在 P006 的短波/长波路径中，按 surface_type 切换 α_water / α_land / α_ice，触发冰-反照率反馈；
  - 相变（冻结/融化）优先消耗/释放净能量（§2.3），与 P006 地表能量方程一致耦合（SW/LW/SH/LH）。
- P008（大气湿度 q 与 E–P–LH 闭环）
  - 蒸发 E = ρ_a C_E |V| (q_sat(T_s) − q) 作为 P006 地表潜热通量 LH = L_v·E 的来源；
  - 超饱和凝结 P_cond 导致大气潜热释放 LH_release = L_v·P_cond，作为 P006 大气能量方程源项；
  - E 与 P_cond 同时驱动 q 的源/汇，云–辐射可在 M4 以湿度/降水进一步调制（光学厚度/发射率）。
- P009（行星水循环闭合）
  - 定义水库（CWV/海冰/陆地贮水/积雪），实现 E−P−R 全球闭合诊断；
  - 陆地桶模型与径流 R 回海洋，长期平均满足 ⟨E⟩ ≈ ⟨P⟩ + ⟨R⟩；
  - 与 P006 的能量闭合共同输出诊断（TOA/SFC/ATM 与 LH/LH_release 一致性）。

集成主循环（建议顺序）：
1) 云/降水诊断（P003/P005，与 P008 的 P_cond 融合）
2) 反照率与辐射（P006；α 由 P005/P007 与云共同决定），计算 ISR/SW/LW
3) 动力步（风/厚度）
4) 湿度步（P008）：E（LH）与 P_cond（LH_release）、q 更新
5) 地表/大气能量步（P006）：含 LH/LH_release；海冰相变（P007）
6) 水文步（P009）：P 的相态分配、陆地桶/径流更新、闭合诊断

## 12. 路线图与里程碑重排（与 P007–P009 联动）

- M1（完成）：辐射核心 + 地表能量方程，SH/LH=0（本文件 §状态）
- M2：引入感热通量 SH（P006 §3.3），并允许 C_H 随地表/粗糙度（P005/P007）调节
- M3：引入潜热环（P008）
  - 地表：LH = L_v·E（E 见 P008 §2.2）进入 P006 地表能量；
  - 大气：LH_release = L_v·P_cond（P008 §2.3）进入 P006 大气能量；
  - 初步校核 ⟨LH⟩ ≈ ⟨LH_release⟩。
- M4：闭合诊断（能量 + 水量）
  - P006：TOA/SFC/ATM 近守恒（|净| < 阈值）；
  - P009：长期平均 ⟨E⟩ ≈ ⟨P⟩ + ⟨R⟩，水库总量 d/dt ≈ 0。
- M5：参数标定与默认组
  - 与 P007（H_mld、α_ice、T_freeze）、P008（C_E、τ_cond、L_v）、P009（τ_runoff、雪阈值/融化率）联合扫描，固化稳定默认。

## 13. 运行示例（全集成占位）

```bash
# 核心能量+诊断
export QD_ENERGY_W=1
export QD_ENERGY_DIAG=1

# P007：平板海洋 / 海冰（建议温和参数）
export QD_MLD_M=50
export QD_ALPHA_WATER=0.08
export QD_ALPHA_ICE=0.60
export QD_T_FREEZE=271.35
export QD_USE_SEAICE=1

# P008：湿度/蒸发/凝结（示例参数）
export QD_CE=1.3e-3
export QD_LV=2.5e6
export QD_Q_INIT_RH=0.5
export QD_TAU_COND=1800
export QD_HUMIDITY_DIAG=1

# P009：水文闭合（桶/径流/雪）
export QD_RUNOFF_TAU_DAYS=10
export QD_SNOW_THRESH=273.15
export QD_SNOW_MELT_RATE=5
export QD_WATER_DIAG=1

# 外部地形与反照率底图（可选）
export QD_TOPO_NC=$(ls -t data/*.nc | head -n1)
export QD_USE_TOPO_ALBEDO=1

python3 -m scripts.run_simulation
```
