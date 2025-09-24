# Project 011: 动态洋流（最小海洋环流）——风驱动浅水海洋与热输送

本项目在 P007 平板海洋（混合层）基础上，引入“会流动的海洋”。目标是在不引入复杂三维密度分层与深海动力的前提下，用单层风驱动浅水海洋模型捕捉行星尺度的水平热量输送，并与现有 P006（能量）、P008（湿度）、P009（水循环）保持一致耦合。

## 0. 状态（2025-09-21）
- [x] 文档与方案定稿（本文件）
- [x] M1：海洋浅水动力核心（uo, vo, η）与数值步进
- [x] M2：风应力强迫与底摩擦、边界条件（海陆掩膜）
- [x] M3：海表温度（SST）平流与垂直净热通量耦合（与 P006/P007 对齐）
- [ ] M4：联调与验证（风带驱动副热带/副极地大环流、能量/水量诊断）
- [ ] M5：参数扫描与默认组固化（C_D、r_bot、K_h、σ₄_ocean 等）

### 实现状态更新（2025-09-21）
- 新增模块 `pygcm/ocean.py`，实现类 `WindDrivenSlabOcean`（M1–M3）：
  - M1（动力核心）：barotropic 浅水动力，预报 `uo, vo, η`，球面有限差分，连续方程 `∂η/∂t = -H ∇·v`。
  - M2（强迫与边界）：风应力 `τ=ρ_a C_D |V_a| V_a`、底摩擦 `r_bot`、海陆边界（陆格速度钳制 0），∇⁴ 超扩散 + 可选 Shapiro 滤波。
  - M3（SST 耦合）：SST 半拉氏平流 + 水平扩散 `K_h ∇²Ts`，可选垂直热通量 `Q_net/(ρ_w c_p,w H)`（默认关闭，避免与 P006 双计）。
- 主循环接入（`scripts/run_simulation.py`）：
  - 新增开关 `QD_USE_OCEAN=1` 启用动态海洋；初始化 `WindDrivenSlabOcean(grid, land_mask, H)`。
  - 每步计算 `Q_net = SW_sfc − LW_sfc − SH − LH`（复用 P006 能量模块），调用 `ocean.step(dt, u_atm, v_atm, Q_net, ice_mask)`。
  - 开阔海（无冰）用 `ocean.Ts` 回写 `gcm.T_s`，与能量/湿度路径一致。
- 数值稳定与保护：
  - 海洋 CFL 子步：`QD_OCEAN_CFL` 控制子步目标，显式子步 `n_sub` 自适应。
  - 反噪：`QD_SIGMA4_OCEAN`、`QD_OCEAN_DIFF_EVERY`、`QD_OCEAN_K4_NSUB`；可选 `QD_OCEAN_SHAPIRO_*`。
  - SST 物理钳制：`QD_TS_MIN`、`QD_TS_MAX`（默认 150–340 K；建议开发期 260–305 K）。
- 烟测运行（安全默认）：
  ```
  QD_USE_OCEAN=1 QD_ENERGY_W=0 QD_USE_SEAICE=0 \
  # 单位：行星日
  QD_SIM_DAYS=0.02 QD_DT_SECONDS=120 QD_OCEAN_CFL=0.4 \
  QD_DIFF_ENABLE=1 QD_FILTER_TYPE=combo QD_SIGMA4=0.02 QD_SHAPIRO_EVERY=6 QD_SHAPIRO_N=2 \
  QD_TS_MIN=260 QD_TS_MAX=305 python3 -m scripts.run_simulation
  ```
  输出包含 `[OceanDiag]`（KE、Umax、η统计、CFL）与 TrueColor 图；模型稳定完成短程自检。

### 增补（2025-09-21 晚）：海–气耦合、可视化与稳定性更新

- 海–气热交换（默认开启）：在 `ocean.step` 中启用 `Q_net/(ρ_w c_p,w H)` 垂直热通量（`QD_OCEAN_USE_QNET=1`），其中  
  `Q_net = SW_sfc − LW_sfc − SH − LH` 来自 P006 能量模块。该通量直接更新海表温度（SST）。  
- 海–气动量交换（相对风）：风应力改为使用相对风  
  `τ = ρ_a C_D |V_a − U_o| (V_a − U_o)`，避免当海流接近大气风速时仍持续注入不合理动量。  
  相关参数：`QD_CD`、`QD_TAU_SCALE`（默认 0.2）、`QD_WIND_STRESS_VCAP`（默认 15 m/s）。  
- 可视化改进：  
  - 状态图使用流线显示大气与洋流，并统一速度色标；三个温度面板（Ts, Ta, SST）统一色表与刻度。  
  - 降水改为“行星日累计”（mm/day），便于与风场/地形对比。  
- 数值稳定与物理保护：  
  - 极区 sponge 阻尼（默认 |φ|≥70° 渐加强化）；极区度量保护（cosφ 下限）。  
  - 海流异常点采用四邻域均值回填（`QD_OCEAN_OUTLIER=mean4`，默认），并保留轻微矢量限幅（`QD_OCEAN_MAX_U` 默认 3 m/s）。  
  - 继续使用 ∇⁴ 选择性耗散与可选 Shapiro 滤波。  
- 默认行为变更：动态海洋默认开启（`QD_USE_OCEAN=1`；显式设 0 可关闭）。  
- 环境变量小结（新增/重要）：  
  - `QD_OCEAN_USE_QNET`（1：启用海–气热交换；默认 1）  
  - `QD_TAU_SCALE`（相对风应力效率，默认 0.2）、`QD_WIND_STRESS_VCAP`（应力风速上限，默认 15）  
  - `QD_OCEAN_OUTLIER=mean4|clamp`（默认 mean4）  
  - `QD_OCEAN_MAX_U`（海流矢量限幅，默认 3 m/s）

## 1. 目标
- 在球面网格上实现单层风驱动浅水海洋（barotropic/shallow-water）：
  - 预报量：海洋速度 uo、vo 与海表高度 η
  - 强迫：近地面大气风产生的风应力 τ_wind，底摩擦 r_bot，尺度选择性耗散
- 将海表温度 Ts 从“仅热力积分”升级为“平流+热力”的示踪物：
  - 水平平流由洋流 (uo, vo) 驱动
  - 垂直净热通量 Q_net（短波/长波/感热/潜热）来自 P006/P007
- 与现有模块一致耦合与闭合诊断（能量与水量）

## 2. 背景与动机（最小原则）
- P007 提供混合层热容量 C_s，使地表具备热惯性，但海洋静止，无法进行水平热再分配。
- 行星尺度气候受海洋输送显著影响：热带–高纬热量输送、东西岸温差（西边界强化）、海冰边界反馈等。
- 本项目仅引入“单层、风驱动、浅水+SST 平流”的最小海洋环流，避免三维密度分层与复杂湍混过程，满足“最小但自洽”的设计哲学。

## 3. 模型设计与方程

### 3.1 状态变量
- uo(lat, lon, t)、vo(lat, lon, t)：海洋水平速度（m s⁻¹）
- η(lat, lon, t)：海表高度异常（m）
- Ts(lat, lon, t)：海表温度（K，示踪物；陆地/海冰区按地表类型处理）

### 3.2 控制方程（球面近似）

- 动量（向量形式）：
  d(vo)/dt = −(vo·∇)vo − f k×vo − g∇η + τ_wind/(ρ_w H) − r_bot·vo − 𝔻₄(vo)

- 连续性：
  ∂η/∂t = −H ∇·vo

- SST 输运（示踪物）：
  ∂Ts/∂t = −(vo·∇)Ts + Q_net/(ρ_w c_p,w H) + K_h ∇²Ts − 𝔻₄(Ts)

说明：
- f：科里奥利参数；g：重力加速度；H：混合层有效厚度（优先复用 P007 的 QD_MLD_M）
- τ_wind：大气风应力；ρ_w、c_p,w：海水密度/比热
- r_bot：底摩擦（线性阻尼）
- 𝔻₄(·)：四阶超扩散（可选，用于抑制网格噪声，参照 P010 思路）

### 3.3 风应力与耦合
- τ_wind = ρ_a C_D |va| va
  - va：近地面大气风（可用 GCM 的 u, v 近地层代表）
  - C_D：拖曳系数（~1.0–2.0×10⁻³，可调）
- Q_net = SW_sfc − LW_sfc − SH − LH（来自 P006，受云/反照率/海冰影响）
- Ts 在海冰覆盖处按 P007 海冰方案处理（抑制蒸发，反照率切换，可能近冻结点钳制）

### 3.4 数值与稳定性（最小实现）
- 网格：复用 SphericalGrid（与大气同网格，避免插值复杂度）
- 时间推进：显式步进（与现有框架一致），满足 CFL：
  - dt_ocean ≤ CFL · min(Δx/|vo|max, √(Δx/gH))，建议 CFL≈0.5
- 耗散：尺度选择性耗散（∇⁴，σ₄_ocean 自适应）+ 低频 Shapiro（可选），参考 P010
- 边界条件：
  - 海/陆边界：陆地像元强钳制 vo=0（或高摩擦区），η 不更新（或强阻尼）
  - 经度周期、极区数值稳健同大气

## 4. 耦合接口与数据流（与 P006/P007/P008/P009 一致）

主循环建议顺序（与现框架一致的小改）：
1) 轨道几何/辐射：计算 I、albedo、SW/LW（P006）
2) 动力学（大气）：u, v, h 等更新（含反噪，P010）
3) 湿度/云/降水：E、P_cond、LH/LH_release（P008），云–辐射一致性
4) 海洋步（本项目）：
   - 计算风应力 τ_wind(u_atm, v_atm)
   - 推进 uo, vo, η（风应力/科氏/压力梯度/底摩擦/超扩散）
   - 对 Ts 做一次“平流子步”：Ts ← Advect(Ts; uo, vo, dt)
5) 地表/大气能量（P006）：
   - 用“已平流后的 Ts”计算表面通量，进行 Q_net 能量积分（海/陆/冰各自参数；海洋使用 C_s_ocean = ρ_w c_p,w H）
   - 大气能量注入 SW_atm/LW_atm/SH/LH_release
6) 水文（P009）：相态分配、陆地桶与径流 R、闭合诊断
7) 诊断与输出：TOA/SFC/ATM 能量闭合、E–P–R 与水库变化、海洋 KE/η/SST 统计

注：SST 的“平流→热力”的算子分裂顺序有助于维持能量路径的直觉一致性（先搬运水，再收支）。

## 5. 代码结构与 API（建议）

新增模块：`pygcm/ocean.py`
- `class WindDrivenSlabOcean`:
  - `__init__(grid, land_mask, H_m, params)`
  - `step(dt, u_atm, v_atm, Q_net, ice_mask=None)`：
    - 计算 τ_wind(u_atm, v_atm)
    - 推进 uo, vo, η（含 r_bot 与 ∇⁴/滤波）
    - 平流 Ts（半拉氏或稳定的双线性流随插值），返回 Ts_adv
  - `diagnostics()`：返回 KE、η_stats、|vo|max、CFL 指标等

与现有脚本对接（`scripts/run_simulation.py`）：
- 读取开关 `QD_USE_OCEAN=1` 则初始化 `WindDrivenSlabOcean`
- 运行环内：获取 GCM 近地风、Q_net、ice_mask，调用 `ocean.step`
- 用 ocean.Ts_adv 进入 P006 地表能量积分（海、陆、冰分类型）

## 6. 参数与环境变量（建议默认）
- 开关与维度：
  - QD_USE_OCEAN（默认 0=关闭；1=开启）
  - QD_OCEAN_DIAG（默认 1，打印海洋诊断）
- 混合层厚度与热容量：
  - QD_MLD_M（默认 50 m，复用 P007；可被 QD_OCEAN_H_M 覆盖）
  - QD_OCEAN_H_M（可选，优先于 QD_MLD_M）
- 风应力与摩擦：
  - QD_CD（默认 1.5e-3）风应力拖曳系数
  - QD_R_BOT（默认 1.0e-6 s⁻¹）底摩擦
- 水平混合与反噪：
  - QD_KH_OCEAN（默认 5.0e3 m² s⁻¹）SST 水平扩散
  - QD_SIGMA4_OCEAN（默认 0.02）∇⁴ 无量纲强度（自适应到 K₄）
  - QD_OCEAN_DIFF_EVERY（默认 1）∇⁴ 施加频率
  - QD_OCEAN_K4_NSUB（默认 1）∇⁴ 子步
  - 备选滤波：QD_OCEAN_SHAPIRO_N（默认 0 关闭）、QD_OCEAN_SHAPIRO_EVERY
- 数值保护：
  - QD_OCEAN_CFL（默认 0.5）建议 CFL
  - QD_OCEAN_MAX_U（默认 3e2 m s⁻¹ 上限钳制，极值保护）

## 7. 任务拆解与里程碑

- M1 动力核心
  - [ ] 实现球面 ∇、∇·、∇² 与 ∇⁴（重用/适配 P010 方法）
  - [ ] uo、vo、η 显式推进与极区稳定化；海/陆边界处理
- M2 强迫与边界
  - [ ] τ_wind 块体公式（近地风→应力）
  - [ ] 底摩擦 r_bot 与海岸阻尼；诊断 CFL 与最大流速
- M3 SST 平流 + 垂直热通量
  - [ ] Ts 平流（半拉氏或稳定双线性）与 K_h 扩散
  - [ ] 与 P006 地表能量积分串联（先平流，后能量）
- M4 联调与验证
  - [ ] 副热带/副极地大环流是否出现（以理想化纬向风带驱动）
  - [ ] 与 P006/P007/P008/P009 闭合诊断一致稳定
- M5 参数标定
  - [ ] 扫描 C_D、r_bot、K_h、σ₄_ocean，固化默认参数组
  - [ ] 输出 KE/η/SST 的统计与图像

## 8. 约束与一致性
- 物理方向性：
  - 风应力驱动“西边界强化”的大环流格局（西岸暖流/东岸寒流雏形）
  - Ts 平流导致热带热量北（南）运，高纬 SST 回升，海冰边界更合理
- 能量一致性（与 P006）：
  - 海洋仅在“地表能量方程”通过 Q_net 与大气热交换
  - Ts 平流不直接改变系统总能量（仅搬运），数值上需避免过度扩散
- 水量一致性（与 P009）：
  - P−E 的净淡水通量目前仅影响能量（通过 LH），对密度不显式耦合
  - 后续如需盐度/密度效应，再扩展（保持最小原则）

## 9. 测试与验收标准
- 功能性：
  - [ ] QD_USE_OCEAN=1 时模型稳定运行若干行星日，输出诊断与图像
  - [ ] 关闭时回退至 P007 原有行为
- 物理合理性：
  - [ ] 在理想化风带强迫下，出现副热带/副极地大环流（η 场闭合环流、KE 合理）
  - [ ] SST 出现热带→中高纬热量输送迹象；东西岸温差方向正确
- 守恒与诊断：
  - [ ] TOA/SFC/ATM 能量诊断长期平均近守恒（阈值如 |净| < 2 W m⁻²）
  - [ ] E–P–R 与水库变化闭合误差在阈值内（与 P009 一致）
- 数值稳定性：
  - [ ] CFL 充足，未见网格噪音或数值爆裂；∇⁴/滤波参数温和可控

## 10. 运行示例
- 开启动态海洋（使用 data 下最新地形/海陆掩膜；温和参数，能量诊断开启）：
```bash
export QD_TOPO_NC=$(ls -t data/*.nc | head -n1)
export QD_USE_TOPO_ALBEDO=1

export QD_USE_OCEAN=1
export QD_OCEAN_DIAG=1
export QD_MLD_M=50
export QD_CD=1.5e-3
export QD_R_BOT=1.0e-6
export QD_KH_OCEAN=5.0e3
export QD_SIGMA4_OCEAN=0.02

export QD_ENERGY_W=1
export QD_ENERGY_DIAG=1

python3 -m scripts.run_simulation
```

- 理想化风带试验（验证副热带/副极地大环流；可在脚本中提供选项或简化风场）：
```bash
export QD_USE_OCEAN=1
export QD_OCEAN_DIAG=1
# 可在 run_simulation 中切换理想化风（占位示例）
# export QD_USE_IDEAL_WIND=1
python3 -m scripts.run_simulation
```

## 11. 与其他项目的交叉引用
- P006（能量收支）：Q_net 提供海–气热交换；SST 影响 SH/LH/辐射；能量诊断统一
- P007（平板海洋/海冰）：复用混合层厚度 H 与海冰抑制蒸发/反照率切换；SST 与海冰边界相互作用
- P008（湿度 q）：SST 改变蒸发 E，影响 q 与降水；LH/LH_release 能量一致
- P009（水循环）：E−P 闭合与水库诊断不变；动态海洋不引入额外水库
- P010（反噪）：海洋同样需要 ∇⁴/滤波以抑制高频噪声（提供独立 ocean_* 参数）

## 12. 后续扩展（不纳入本里程碑）
- 分层/密度（baroclinic，多层浅水或简化密度层）
- 盐度/淡水通量对密度的影响（E−P、融化/冻结）
- 海冰动力（风应力作用于海冰、冰–海相对运动）
- 等密度混合/等熵混合的参数化
- 河流入海与沿岸流（与 P009 更细的径流路径）

---

附注：本项目在不引入三维复杂性的前提下，优先恢复“风驱动水平热输送”这一关键机制。实现后，预计可显著改善热带过暖/高纬偏冷与东西岸差异不足的问题，并为进一步的海冰–海洋–大气反馈研究奠定基础。
