# 7. 平板海洋与海冰、风驱动浅水海洋（Slab Ocean, Sea Ice & Wind-Driven Ocean）

本章整合项目 P007（平板海洋与海冰）、P011（动态洋流：风驱动浅水海洋）与 P012（极点处理）的设计与实现，形成面向 Qingdai GCM 的海–气耦合与数值稳定的一体化说明。

状态（2025-09-21）
- P007：已实现混合层等效热容量地图（C_s_map）与最小海冰热力学（冻结/融化、反照率切换、低热容），与能量/湿度路径一致；诊断完善中
- P011：已实现风驱动浅水海洋（uo, vo, η）与 SST 平流，海–气热/动量耦合，CFL 子步、自适应抑噪；参数标定与回归中
- P012：已在海洋模块接入极点修正（标量/矢量极点一致化），默认启用

关联文档
- docs/06-energy-framework.md：能量收支与 Q_net、SH/LH/LH_release
- docs/08-humidity-and-clouds.md：蒸发 E 与潜热一致性
- docs/04-runtime-config.md：运行时环境变量总表（含海洋/海冰条目）
- projects/007-slab-ocean.md、projects/011-ocean-model.md、projects/012-polar-treatment.md


## 7.1 目标与范围

- 为地表提供真实的热惯性（混合层海洋）与冰–反照率反馈（最小海冰）
- 引入风驱动单层浅水海洋（barotropic）以实现行星尺度的水平热输送（SST 平流）
- 与能量（P006）、湿度（P008）、水文（P009）保持能量—水量一致性与守恒诊断
- 在极点保证数值与物理一致（标量/矢量极点平均回填）


## 7.2 平板海洋（Slab / Mixed-Layer Ocean, MLO）

等效热容量（每单位面积）
- 海洋：C_s,ocean = ρ_w · c_p,w · H_mld  
  ρ_w≈1000 kg m⁻³，c_p,w≈4200 J kg⁻¹ K⁻¹，H_mld（默认 50 m）
- 陆地：C_s,land（默认 3e6 J m⁻² K⁻¹）
- 海冰：C_s,ice（默认 5e6 J m⁻² K⁻¹，薄层近地表）

C_s_map 构建
- 根据 land_mask 生成海/陆热容量地图；若存在海冰掩膜则对冰上采用 C_s,ice
- 地表能量方程见 docs/06：C_s dT_s/dt = SW_sfc − LW_sfc − SH − LH


## 7.3 最小海冰热力学与反照率切换

冻结/融化优先（相变能量）
- 当开阔海接近冻结点且 Q_net < 0：以相变生成冰厚 Δh_i = −Q_net·dt/(ρ_i L_f)
- 有冰区域且 Q_net > 0：优先融化冰，Δh_i 同上取负；T_s 不越过冻结点
- 参数：T_freeze（默认 271.35 K）、ρ_i（917 kg m⁻³）、L_f（3.34e5 J kg⁻¹）

光学冰覆盖
- ice_frac = 1 − exp(−h_ice/H_ref)，H_ref 默认 0.5 m
- 反照率合成：α_base 使用 α_ice 在冰覆盖处替代水面 α_water，再与云量混合（见 docs/06）

蒸发抑制
- 海冰区域蒸发权重显著降低（QD_ICE_EVAP_SCALE），见 docs/08


## 7.4 风驱动浅水海洋（Wind-Driven Slab Ocean）

状态变量
- uo、vo：海洋水平速度（m s⁻¹）
- η：海表高度异常（m）
- Ts_ocean：SST（作为示踪物，海/陆/冰类型下分别处理）

控制方程（球面近似）
- 动量：d(vo)/dt = −(vo·∇)vo − f k×vo − g∇η + τ/(ρ_w H) − r_bot·vo − 𝔻₄(vo)
- 连续：∂η/∂t = −H ∇·vo
- SST：∂Ts/∂t = −(vo·∇)Ts + Q_vert/(ρ_w c_p,w H) + K_h ∇²Ts − 𝔻₄(Ts)

海–气耦合
- 风应力（相对风）：τ = ρ_a C_D |V_a − U_o| (V_a − U_o)，限定风速上限（VCAP）
- 垂直热通量：Q_vert = Q_net = SW_sfc − LW_sfc − SH − LH（默认启用注入海洋）
- 海冰掩膜：冰下热通量按比例减弱（QD_OCEAN_ICE_QFAC）

数值与稳定
- 显式步进 + 自适应子步，满足 CFL（QD_OCEAN_CFL）；风应力效率（QD_TAU_SCALE）
- 选择性抑噪：∇⁴ 超扩散（QD_SIGMA4_OCEAN、QD_OCEAN_K4_NSUB、QD_OCEAN_DIFF_EVERY）+ 可选 Shapiro（QD_OCEAN_SHAPIRO_*）
- 极值保护：速度限幅（QD_OCEAN_MAX_U）、异常像元回填（QD_OCEAN_OUTLIER=mean4|clamp）


## 7.5 极点处理（Polar Corrections，P012）

动机：规则经纬度网格在极点同一物理点对应多个经度样本；标量/矢量需要一致化。

- 标量（如 SST）：对极圈海洋格点取均值，并回填该圈海洋格点
- 矢量（uo, vo）：在统一极点切平面中将各经度的（east, north）分量映射到共同基底后做向量平均，再映射回各自经度的局地基底后回填
- 开关：QD_OCEAN_POLAR_FIX=1（默认启用）；整圈全陆地时跳过


### 7.5.1 极点硬约束（能量路径，P006）

- 动机：在规则经纬度网格上，极点行（物理上单点）在数值上表现为一整圈经度样本；再叠加极点平均/一致化等处理，易将局部的“暖极”数值伪影放大为“南极无冰、周围有冰”的不物理状态。
- 实现：在 energy 模块的 integrate_surface_energy_with_seaice 中加入南极行（lat = −90°）的硬约束：
  - 条件：海洋格点且净热通量为冷却 Q_net < 0 且更新后 Ts_next > t_freeze
  - 动作：将该格点的 Ts_next 钳制为 t_freeze，从而抵消数值性加热伪影，允许随后的相变公式在下一个时间步自然增长海冰。
- 开关：QD_POLAR_FREEZE_FIX（默认 1=启用）；仅作用于南极行（grid.lat 从 −90→+90 递增时的第 0 行），后续可按需扩展到北极。
- 相容性：不直接修改 h_ice；海冰相变、反照率与蒸发路径保持一致，能量与水量闭合诊断不受破坏。

## 7.6 接口与主循环集成（pygcm/ocean.py、scripts/run_simulation.py）

- 初始化：WindDrivenSlabOcean(grid, land_mask, H_m)；与 GCM 共用网格以免插值
- 每步：
  1) 计算风应力 τ(u_atm, v_atm, U_o) 与 Q_net（来自能量模块）
  2) 推进 uo, vo, η（含底摩擦、抑噪、CFL 子步、极区 sponge/极点一致化）
  3) 平流 SST 并应用垂直热通量更新（开阔海）；海冰/陆地按各自规则处理
  4) 将开阔海 SST 回写地表 Ts（供能量/湿度路径使用）
- 诊断：OceanDiag 打印 KE、|U|max、η 统计、CFL、⟨SST⟩ 等


## 7.7 环境变量（节选；详见 docs/04-runtime-config.md）

混合层与热容
- QD_MLD_M（默认 50 m）、QD_CS_LAND（3e6）、QD_CS_ICE（5e6）
- QD_RHO_W（1000）、QD_CP_W（4200）

海/冰反照率与冻结
- QD_ALPHA_WATER（0.08）、QD_ALPHA_ICE（0.60）
- QD_T_FREEZE（271.35 K）、QD_RHO_ICE（917）、QD_LF（3.34e5）、QD_HICE_REF（0.5 m）
- QD_USE_SEAICE（1）

海–气动量/热交换与稳定
- QD_USE_OCEAN（1）、QD_OCEAN_USE_QNET（1）
- QD_CD（1.5e-3）、QD_R_BOT（2e-5 s⁻¹）
- QD_TAU_SCALE（0.2）、QD_WIND_STRESS_VCAP（15 m s⁻¹）
- QD_KH_OCEAN（5.0e3 m² s⁻¹）、QD_SIGMA4_OCEAN（0.02）、QD_OCEAN_DIFF_EVERY（1）、QD_OCEAN_K4_NSUB（1）
- QD_OCEAN_SHAPIRO_N（0 关闭）、QD_OCEAN_SHAPIRO_EVERY
- QD_OCEAN_CFL（0.5）、QD_OCEAN_MAX_U（3 m s⁻¹）、QD_OCEAN_OUTLIER（mean4|clamp）
- QD_OCEAN_POLAR_SPONGE_LAT（70°）、QD_OCEAN_POLAR_SPONGE_GAIN（5e-5 s⁻¹）
- QD_OCEAN_POLAR_FIX（1）


## 7.8 推荐默认与运行示例

启用海洋与海冰、能量诊断
```bash
export QD_TOPO_NC=$(ls -t data/*.nc | head -n1)
export QD_USE_TOPO_ALBEDO=1

export QD_USE_OCEAN=1
export QD_OCEAN_USE_QNET=1
export QD_USE_SEAICE=1
export QD_MLD_M=50

export QD_ENERGY_W=1
export QD_ENERGY_DIAG=1

# 抑噪与极点
export QD_FILTER_TYPE=combo
export QD_SIGMA4_OCEAN=0.02
export QD_OCEAN_DIFF_EVERY=1
export QD_OCEAN_POLAR_FIX=1

python3 -m scripts.run_simulation
```

快速烟测（安全步长/限幅）
```bash
# 单位：行星日
export QD_SIM_DAYS=0.02
export QD_DT_SECONDS=120
export QD_OCEAN_CFL=0.4
export QD_TS_MIN=260
export QD_TS_MAX=305
python3 -m scripts.run_simulation
```


## 7.9 诊断与验收标准（建议）

- 功能性：启用海洋后稳定运行若干“青黛日/年”，输出 OceanDiag 与状态图
- 物理合理性：风带驱动下出现副热带/副极地大环流迹象；SST 出现热带→中高纬热输送与东西岸差异
- 守恒与一致性：TOA/SFC/ATM 能量多年平均近守恒；⟨LH⟩≈⟨LH_release⟩；E–P–R 水量闭合误差在阈值内
- 数值稳定性：CFL 充足，无高频条纹或爆裂；极点一致化后极圈场无经向条纹


## 7.10 变更记录（Changelog）

- 2025‑09‑20：完成海洋热容地图与最小海冰（P007 M1–M3）；SST 平流与 Q_net 注入联通
- 2025‑09‑21：风应力改为相对风；极区 sponge 与极点一致化；OceanDiag 与 TrueColor/状态图改进
- 2025‑09‑21：文档迁移与整合至 docs/07‑ocean‑and‑sea‑ice.md；与 04/06/08 章节交叉引用对齐
