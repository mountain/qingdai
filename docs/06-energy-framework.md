# 6. 动态能量收支框架（Energy Budget Framework）

本章系统化迁移项目 P006（物理引擎升级：从牛顿冷却到动态能量收支）的设计与实现，给出清晰的物理结构、模块接口、诊断方法与推荐参数。完成本章后，可将能量路径作为一等公民接入 GCM 主循环，并与湿度（水汽潜热）、海洋（SST 与 Q_net）、海冰（相变）保持一致耦合。

状态（2025-09-21）
- 已完成（M1–M3）：短波/长波辐射、地表能量方程、感热通量 SH、潜热通量 LH 与潜热释放 LH_release、与动态海洋 Q_net 耦合
- 进行中（M4–M5）：能量闭合诊断稳态验收、参数标定与默认组固化
- 新增：温室参数微调（Autotune）、统一诊断打印，TrueColor/状态图面板统一色表

关联文档
- docs/03-climate-model.md：整体气候框架与强迫
- docs/04-runtime-config.md：运行时环境变量总表（含能量模块条目）
- docs/07-ocean-and-sea-ice.md：混合层海洋/海冰与 Q_net
- docs/08-humidity-and-clouds.md：蒸发 E、凝结 P_cond 与潜热一致性
- projects/006-energy-budget.md：原始项目文档（细节与路线图）


## 6.1 目标与原则

- 用“显式能量通量”替代或补充牛顿冷却，使辐射（SW/LW）、边界层（SH/LH）、潜热释放（LH_release）在地表与大气两条能量方程中自洽出现。
- 可诊断、可守恒、可扩展：统一的 TOA/SFC/ATM 收支打印，长期平均接近守恒；便于逐步加入海冰相变、海洋热输送、云光学厚度等过程。
- 平滑迁移：提供混合权重（QD_ENERGY_W∈[0,1]），支持从 Teq 式逐步切换到能量收支式。


## 6.2 变量、通量与方程

状态变量（核心）
- Ts：地表温度（海/陆/冰类型区分）
- Ta：大气温度（单层诊断/与 h 的代理关系）；实现中由 energy_height 路径间接反映
- C：云量（来自降水驱动与代理源的融合）
- q：比湿（docs/08-humidity-and-clouds.md）

通量分量（W m⁻²）
- 短波：SW_atm（大气吸收）、SW_sfc（地表吸收）
- 长波：LW_atm（作用于大气的净长波）、LW_sfc（地表净长波，=DLR−σTs⁴）
- 边界层：SH（感热）、LH（潜热，=L_v·E）
- 潜热释放：LH_release（=L_v·P_cond，进入大气能量）
- 表面净热：Q_net = SW_sfc − LW_sfc − SH − LH（用于海洋垂直通量）

地表能量方程（每格点）
- C_s dTs/dt = SW_sfc − LW_sfc − SH − LH
  - C_s：等效热容量地图（海/陆/冰不同；见 docs/07-ocean-and-sea-ice.md）
  - 海冰：相变优先（冻结/融化），详见 7 章

大气能量方程（单层近似）
- c_p ρ_a H dTa/dt = SW_atm + LW_atm + SH + LH_release
  - H：单层等效厚度（QD_ATM_H）；实现中还可间接通过位势高度 h 的倾向体现


## 6.3 模块与接口（pygcm/energy.py）

- shortwave_radiation(I, albedo, cloud, params) → (SW_atm, SW_sfc)
  - I：两星几何入射（见 docs/02-orbital-dynamics.md / forcing）
  - 反射 R = I·α_total；大气吸收 SW_atm = I·A_sw0·(1 + k_sw_cloud·C)
- longwave_radiation(Ts, Ta, cloud, params) → (LW_atm, LW_sfc)
  - ε_eff = clip(ε0 + k_lw_cloud·C, 0, 1)，DLR = ε_eff σ Ta⁴，LW_sfc = DLR − σ Ts⁴
- boundary_layer_fluxes(Ts, Ta, u, v, land_mask, params) → (SH, LH_partial?)
  - 计算 SH，LH 由湿度模块蒸发 E 提供（LH = L_v·E）
- integrate_surface_energy(Ts, SW_sfc, LW_sfc, SH, LH, dt, params | C_s_map)
  - Map 版本支持逐格 C_s；海冰相变见 7 章
- integrate_atmos_energy_height(..., SW_atm, LW_atm, SH, LH_release, dt, params)
  - 将通量转为大气能量高度/温度路径的倾向

实现要点
- 夜侧/高原保护：QD_T_FLOOR（默认 150 K）可限制极端冷塌
- 隐式/半隐式更新：对 σTs⁴ 等刚性项采用半隐式处理提升稳定性


## 6.4 反照率合成与云/海冰耦合

总反照率 α_total（示意）
- α_total = α_base(type, Ts)·(1 − C) + α_cloud·C
  - α_base：来自 base_albedo_map（地形/高山/海陆差异；见 docs/05）与海冰 α_ice
  - C：云量（见 docs/08），同时调制 SW/LW 的光学厚度/发射率
- 开关：QD_USE_TOPO_ALBEDO（默认 1）；α_water、α_ice、α_cloud 可在 04 章查表


## 6.5 与湿度/潜热的一致性（见 docs/08）

- 地表潜热 LH = L_v·E（E：块体公式蒸发），从地表能量中扣除
- 大气潜热释放 LH_release = L_v·P_cond（超饱和凝结），进入大气能量
- 期望长期平均：⟨LH⟩ ≈ ⟨LH_release⟩（能量一致）
- 蒸发/凝结与云量：湿度与降水调制云光学厚度（C_eff）→ 反照率/长波路径一致


## 6.6 与海洋/海冰的一致性（见 docs/07）

- 海洋垂直热通量：Q_net/(ρ_w c_p,w H) 注入 SST（QD_OCEAN_USE_QNET=1，默认开）
- 海冰：相变能量优先（冻结/融化），反照率 α_ice 切换，蒸发显著降低
- 海/陆/冰的 C_s、α_base、蒸发权重 S_type 分别处理


## 6.7 诊断与能量闭合

统一诊断打印（建议每 ~200 步）
- [EnergyDiag] TOA_net、SFC_net、ATM_net、⟨Ts⟩、⟨Ta⟩
  - TOA_net = 入射 − 反射 − OLR
  - SFC_net = SW_sfc − LW_sfc − SH − LH
  - ATM_net ≈ TOA_net − SFC_net（或 SW_atm+LW_atm+SH+LH_release 的全球积分）
- 期望长期平均（稳态）：|⟨TOA_net⟩|、|⟨SFC_net⟩|、|⟨ATM_net⟩| → 0（阈值 ~2 W m⁻²）
- 可视化：状态图中统一 Ts/Ta/SST 色表，降水图为“行星日累计（mm/day）”


## 6.8 环境变量（节选，详见 docs/04-runtime-config.md）

主控
- QD_ENERGY_W（0..1，默认 0）：能量路径权重（=1 时完全使用能量收支）
- QD_ENERGY_DIAG（默认 1）：能量诊断打印
- QD_T_FLOOR（默认 150 K）：夜侧温度下限（保护）

短波（SW）
- QD_SW_A0（默认 0.06）大气短波基吸收
- QD_SW_KC（默认 0.20）云短波吸收增益

长波（LW）
- QD_LW_EPS0（默认 0.70）无云大气发射率
- QD_LW_KC（默认 0.20）云长波增强

边界层与鲍文比
- QD_CH（默认 1.5e-3）感热交换系数
- QD_BOWEN_LAND（默认 0.7）、QD_BOWEN_OCEAN（默认 0.3）

大气层等效厚度
- QD_ATM_H（默认 800 m）：将能量源项转为大气温度/高度倾向

自调参（Autotune）
- QD_ENERGY_AUTOTUNE（默认 0 关）：启用后按 TOA_net 微调 ε_lw 与云增益
- QD_ENERGY_TUNE_EVERY（默认 50 步）：调参频率
- 步长：QD_TUNE_RATE_EPS（默认 5e-5）、QD_TUNE_RATE_KC（默认 2e-5）
- 范围约束：ε_lw∈[0.30,0.98]，k_lw_cloud∈[0.0,0.80]


## 6.9 推荐默认与运行示例

温和稳定（能量路径 + 诊断）
```bash
export QD_ENERGY_W=1
export QD_ENERGY_DIAG=1
export QD_T_FLOOR=160

# 短波/长波/边界层（示例）
export QD_SW_A0=0.06
export QD_SW_KC=0.20
export QD_LW_EPS0=0.70
export QD_LW_KC=0.25
export QD_CH=1.5e-3
export QD_BOWEN_LAND=0.7
export QD_BOWEN_OCEAN=0.3

# 外部地形与反照率底图（可选）
export QD_TOPO_NC=$(ls -t data/*.nc | head -n1)
export QD_USE_TOPO_ALBEDO=1

python3 -m scripts.run_simulation
```

自调参与海气耦合（示例）
```bash
export QD_ENERGY_W=1
export QD_ENERGY_AUTOTUNE=1
export QD_ENERGY_TUNE_EVERY=50
export QD_OCEAN_USE_QNET=1
python3 -m scripts.run_simulation
```


## 6.10 数值稳定与调参建议

- 先以混合权重平滑迁移：逐步将 QD_ENERGY_W 从 0 提至 1
- 夜侧/极端冷塌：适度提高 QD_T_FLOOR（例如 155–170 K）防止发散
- SH 与 LH：Bowen 比与 C_H 影响大，可从默认温和值开始；海陆差异必需
- 与反噪（docs/10-numerics-and-stability.md）配合：保持 ∇⁴ 超扩散与 Shapiro 的温和默认，减少条纹伪迹对降水/云/能量的污染
- 与湿度闭合一致性：检查 ⟨LH⟩ ≈ ⟨LH_release⟩ 的长期平均差值（应很小）


## 6.11 验收标准（建议）

- 运行稳定 ≥ 若干“青黛年”，状态图与诊断正常
- 能量闭合：TOA/SFC/ATM 的多年平均净通量 |净| < 2 W m⁻²（阈值可按分辨率/步长调整）
- 潜热一致性：⟨LH⟩ ≈ ⟨LH_release⟩（长期平均）
- 物理合理性：云增多 → 地表短波吸收下降/大气长波增暖；海陆差异体现在 SH/LH 分配
- 输出：日志打印关键参数与诊断；可选保存 NetCDF 汇总


## 6.12 变更记录（Changelog）

- 2025‑09‑20：完成辐射核心与地表能量；引入 SH 与潜热闭环；统一能量诊断  
- 2025‑09‑21：新增温室参数 Autotune；与动态海洋 Q_net 耦合改为默认；可视化/诊断改进  
- 2025‑09‑21：文档迁移与整合至 docs/06‑energy‑framework.md；与 04/07/08 章节交叉引用对齐
