# Project 008: 大气湿度（q）与蒸发—凝结—降水闭环

本项目在 P006 能量收支框架之上，引入单层大气的比湿（specific humidity, q）作为新的预报量，按物理过程构建“蒸发（E）→ 水汽输送/变化 → 凝结/降水（P）→ 潜热释放”的闭环，并与 P007（平板海洋与海冰）和 P009（水循环核算）协同，形成能量—水量自洽的最小系统。

## 0. 状态（2025-09-20）
- [x] M1：新增大气湿度场 q（初始化/平流/扩散），诊断 E 与 q_sat(T)
- [x] M2：蒸发 E 的块体公式接入；潜热通量 LH= L_v·E 进入地表能量方程（P006）
- [x] M3：凝结/降水：q 超饱和 → P_cond；潜热释放 LH_release= L_v·P_cond 加热大气（P006）
      （已完成：实现 energy.integrate_atmos_energy_height，并在 dynamics.SpectralModel 中将 LH_release 连同 SW_atm/LW_atm/SH 注入大气能量积分；由 QD_ENERGY_W 控制耦合强度）
- [x] M4：与云/反照率的一致性（湿度/降水调节云光学厚度，已接入辐射与反照率路径）
- [ ] M5：闭环诊断与守恒校核（与 P009）
      （已加入 ⟨E⟩/⟨P_cond⟩/⟨LH⟩/⟨LH_release⟩ 全局均值打印；完整水量/能量闭合验证待 M3/M4 完成）

### 实现状态更新（2025-09-20）
- 新增模块：
  - pygcm/humidity.py
    - HumidityParams（支持环境变量：QD_CE、QD_RHO_A、QD_MBL_H、QD_LV、QD_P0、QD_OCEAN_EVAP_SCALE、QD_LAND_EVAP_SCALE、QD_ICE_EVAP_SCALE、QD_TAU_COND、QD_HUMIDITY_DIAG）
    - q_sat(T, p)：Tetens 公式饱和比湿
    - q_init(Ts, RH0, p0)：基于相对湿度初始化 q（QD_Q_INIT_RH）
    - surface_evaporation_factor(...)：按海/陆/冰调制蒸发因子
    - evaporation_flux(...)：块体公式蒸发通量 E（kg m^-2 s^-1）
    - condensation(...)：超饱和凝结 P_cond 通量及更新后 q
- 动力核集成：
  - pygcm/dynamics.py（SpectralModel）
    - 初始化 q 与湿度参数；新增诊断：E_flux_last、P_cond_flux_last、LH_last、LH_release_last
    - time_step 中计算 E 与 LH，并将 LH 传入地表能量积分（P006：integrate_surface_energy*）
    - 计算 P_cond 与 LH_release，并通过 energy.integrate_atmos_energy_height 将（SW_atm + LW_atm + SH + LH_release）注入到大气（以 h 表示）的能量路径；支持 QD_ATM_H、QD_ENERGY_W 参数
    - 对 q 采用与 T_s 一致的半拉氏平流和温和扩散
- 运行脚本：
  - scripts/run_simulation.py 新增 HumidityDiag 每≈200步打印：⟨E⟩、⟨P_cond⟩、⟨LH⟩、⟨LH_release⟩（面积权重）
- 快速验证（短程，不出图）：
  ```bash
  QD_ENERGY_W=1 QD_ENERGY_DIAG=1 QD_HUMIDITY_DIAG=1 \
  # 单位：行星日
  QD_SIM_DAYS=0.5 QD_PLOT_EVERY_DAYS=1000 \
  python3 -m scripts.run_simulation
  ```
  注：起步若风速小（|V|≈0）则 ⟨E⟩≈0，随风场自旋上升 E 和 LH 会变为正值。

关联项目：
- P006 动态能量收支：地表 LH 与大气 LH_release 的能量一致性与守恒统计
- P007 平板海洋/海冰：开阔海面主导蒸发；海冰与陆地显著抑制 E
- P009 水循环闭合：E-P 的水量守恒与各水库（大气/海洋/海冰）记账

## 1. 目标
- 将 q 作为 GCM 的新状态变量，随时间演化并被风场平流/扩散。
- 采用块体空气动力学公式计算蒸发 E，驱动地表潜热通量与大气水汽源项。
- 通过超过饱和比湿 q_sat(T_a) 的超饱和量诊断凝结/降水与潜热释放。
- 与 P006 的能量方程耦合，使 LH（地表）与 LH_release（大气）互为能量来源/汇；与 P007 的地表类型/海冰状态调节 E。

## 2. 模型设计与公式

### 2.1 湿度预报方程（单层）
- 变量：q（kg/kg），定义在与 T_s 同样的 2D 网格上，代表近地层平均含湿。
- 预报方程（离散形式）：
  dq/dt = 𝒜(q; u, v) + S_Evap − S_Cond + 𝒟(q)
  - 𝒜：半拉氏平流（与 T_s 同法），数值稳定
  - 𝒟：弱扩散（与现有场一致的温和扩散因子）
  - S_Evap：蒸发源（见 2.2）
  - S_Cond：凝结汇（见 2.3）

### 2.2 蒸发 E 与潜热通量 LH
- 块体公式（近地层一层近似）：
  E = ρ_a · C_E · |V| · (q_sat(T_s) − q_a) · S_type
  - |V| = √(u² + v²) （可用近地风）
  - q_a ≡ q（本单层近似）
  - q_sat(T) 可用 Tetens/Clausius-Clapeyron 近似（温和、连续）
  - S_type：地表类型与海冰调制因子（见 P007），海洋=1，海冰 ≪1，陆地中间或小值
- 地表潜热通量：
  LH = L_v · E
  - L_v：汽化潜热，默认 2.5×10⁶ J kg⁻¹（可随温度微调，首版常数）
  - LH 以“自地表向上”为正，进入 P006 的地表能量方程（从地表能量中扣除）

- 作为大气湿度源项：
  S_Evap = E / h_mbl
  - h_mbl：混合边界层有效高度（参数，默认 500–1000 m）；单层近似下用于换算质量通量→比湿倾向
  - 简化：可用单位换算常数 k_E 代替 h_mbl，将 E 直接转为 dq/dt（见参数）

### 2.3 凝结/降水 P 与潜热释放
- 饱和比湿 q_sat(T_a)：
  - 以大气温度 T_a（可用诊断 Ta 或 Ts 的平滑近似）计算
- 超饱和凝结（时间步 dt）：
  P_cond = max(0, q − q_sat(T_a)) / τ_cond
  - τ_cond：凝结时间尺度（默认接近一个时间步或数个时间步，避免耗尽/数值不稳定）
  - q 更新：q ← q − P_cond · dt
- 潜热释放（加热大气）：
  LH_release = L_v · P_cond
  - 作为对大气能量方程的源（P006：integrate_atmos_energy 的输入项）

- 与降水诊断/云量的关系（兼容现有 P003/P005）：
  - 可将 P_cond 融入现有降水诊断作为“物理上限”或与散度法加权融合；
  - 云量可继续由“降水驱动 + 代理源”方式计算，或引入湿度对云的调制（M4）。

## 3. 接口与数据结构（建议）

### 3.1 新增模块（建议）
- 新建 `pygcm/humidity.py`：
  - `q_init(grid, RH0=0.5, T_a0=280K) -> q0`（根据相对湿度与温度初始化）
  - `q_advect_diffuse(q, u, v, dt, kappa_q) -> q_next`
  - `q_sat(T) -> qsat`（连续、数值稳定）
  - `evaporation_flux(Ts, q, |V|, surface_factor, params) -> E`
  - `condensation(q, T_a, dt, params) -> (P_cond, q_next)`
  - 统一参数结构（C_E, L_v, h_mbl/k_E, τ_cond, kappa_q 等）

### 3.2 与现有模块的耦合点
- `pygcm/energy.py`：
  - 地表能量方程 `integrate_surface_energy*` 增加 LH（由本模块 E 提供）
  - 大气能量收支 `integrate_atmos_energy` 增加 LH_release（由本模块 P_cond 提供）
- `pygcm/dynamics.py` 主循环次序（建议）：
  1) 先更新云/降水诊断（保持现方案或与湿度法融合）
  2) 计算反照率与 ISR、辐射分量
  3) 计算风场与动力步
  4) 湿度步：E（LH）与 P_cond（LH_release）、q 的更新（可与云/降水合并步）
  5) 地表/大气能量积分（包含 LH 与 LH_release）
  6) 可选：云的再次融合/平滑

## 4. 任务拆解
- T1（M1）q 场与运移：实现 q 的初始化、半拉氏平流与温和扩散；输出 q 基本统计
- T2（M2）蒸发与 LH：实现 E（块体公式）与 LH= L_v·E，接入 P006 地表能量
- T3（M3）凝结/降水与 LH_release：实现 P_cond 与 LH_release，接入 P006 大气能量
- T4（M4）云与湿度一致性（可选）：湿度/降水调制云光学厚度/发射率参数（P006）
- T5（M5）闭环诊断：与 P009 联测，E-P 与水量守恒、能量闭合验证

## 5. 参数与环境变量（建议默认）
- 湿度/相变：
  - QD_LV（默认 2.5e6 J kg⁻¹）
  - QD_TETENS_A/B（或采用内置公式常数），qd：q_sat 参数
  - QD_Q_INIT_RH（默认 0.5），QD_TA_INIT（默认 280 K）
  - QD_TAU_COND（默认 单步~数步）
- 蒸发：
  - QD_CE（默认 1.3e-3），QD_MBL_H（默认 800 m）或 QD_E2Q_FACTOR（E→dq/dt 换算）
  - QD_OCEAN_EVAP_SCALE（默认 1.0），QD_LAND_EVAP_SCALE（默认 0.2），QD_ICE_EVAP_SCALE（默认 0.05）
- 云-湿度一致性（M4）：
  - QD_CLOUD_COUPLE（默认 1）：开启湿度/降水耦合云光学厚度（cloud_eff）用于辐射与反照率
  - QD_RH0（默认 0.6）：相对湿度阈值，RH > RH0 时云量按增益系数增强
  - QD_K_Q（默认 0.3）：RH 超阈增益系数（作用于 cloud_eff）
  - QD_K_P（默认 0.4）：P_cond 相对强度增益系数（tanh 归一后作用于 cloud_eff）
  - QD_PCOND_REF（默认 P_cond 正值的中位数）：P_cond 归一化参考尺度
- 数值：
  - QD_Q_DIFF（默认 1e-6–1e-5 等量纲因子），温和扩散强度
- 诊断：
  - QD_HUMIDITY_DIAG（默认 1）：打印/输出 q、E、P_cond、LH、LH_release 的全局均值与守恒核算

## 6. 诊断与输出
- 空间—时间统计：〈q〉、〈E〉、〈P_cond〉、〈E−P〉，Column Water Vapor（CWV）近似
- 能量一致性：〈LH〉（SFC）与〈LH_release〉（ATM）长期均衡（与 P006 的 TOA/SFC/ATM 一同打印）
- 海—陆—冰对比：E 与 q 在不同下垫面上的差异（与 P007 surface_type_map 联动）

## 7. 测试与验收
- [ ] E 随 |V| 与 (q_sat(Ts)−q) 增大而增大；海>陆>冰
- [ ] q 超饱和时出现凝结/降水；LH_release 正加热大气
- [ ] 能量一致性：全球时间平均 LH ≈ LH_release（差异在容差内）
- [ ] 与 P009 的水量守恒：∫(E−P) 与水库变化匹配在阈值内

## 8. 与其他项目的交叉引用
- P006（能量收支）：LH 与 LH_release 的双向能量耦合；诊断统一
- P007（平板海洋/海冰）：surface_type_map 调节 E；海冰抑制蒸发
- P009（水循环）：E 与 P_cond 进入水量守恒核算与库容变化

## 9. 运行示例（占位）
```bash
export QD_ENERGY_W=1
export QD_ENERGY_DIAG=1
# 湿度模块参数（示例）
export QD_CE=1.3e-3
export QD_LV=2.5e6
export QD_Q_INIT_RH=0.5
export QD_TAU_COND=1800
export QD_HUMIDITY_DIAG=1

python3 -m scripts.run_simulation
```

---

注：P008 初期采取“单层湿度 + 单层能量”的最小一致性方案，优先保证能量/水量闭合与数值稳定，再逐步扩展云光学厚度的诊断与湿度—云的强一致性。
