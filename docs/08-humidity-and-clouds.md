# 8. 湿度 q、蒸发–凝结–降水闭环与云–辐射耦合（Humidity, E–P–LH & Clouds）

本章整合项目 P003（云–降水–反照率反馈）与 P008（大气湿度 q 的引入与 E–P–LH 闭环）的设计与实现，确立“湿度—降水—潜热—云—辐射”的一体化路径，并与能量（P006）、海洋/海冰（P007/011）、水循环（P009）保持守恒一致。

状态（2025-09-21）
- P008：已实现单层比湿 q 的初始化、平流与温和扩散；蒸发 E（块体公式）→ LH；超饱和凝结 P_cond → LH_release；与能量路径、海洋/海冰/水循环联动；诊断打印开启
- P003：已将“动力-降水诊断 + 云量融合 + 动态反照率”方案并入主循环，与湿度模块一致性增强（RH、P_cond 对云光学厚增益）

关联文档
- docs/06-energy-framework.md：能量收支、SH/LH、LH_release、Q_net
- docs/07-ocean-and-sea-ice.md：海面蒸发权重（开阔海/海冰/陆地）与 SST 的影响
- docs/09-hydrology-closure.md：水量闭合与 E–P–R
- projects/003-cloud-precipitation-albedo.md、projects/008-humidity.md


## 8.1 目标与原则

- 在单层大气中引入比湿 q 作为新状态量，使蒸发 E、凝结/降水 P_cond 与潜热通量 LH/LH_release 自洽出现
- 建立云–辐射的一致性：云量随降水/湿度增强，削弱短波地表吸收、增强长波温室
- 与能量/海洋/水文的闭环守恒一致：⟨LH⟩≈⟨LH_release⟩，长期 ⟨E⟩≈⟨P⟩+⟨R⟩


## 8.2 状态与方程（单层近地层平均）

湿度预报（离散思路）
- dq/dt = 𝒜(q; u, v) + S_Evap − S_Cond + 𝒟(q)  
  𝒜：半拉氏平流；𝒟：温和扩散  
  S_Evap：蒸发源；S_Cond：凝结汇

蒸发 E（kg m⁻² s⁻¹）与潜热 LH（W m⁻²）
- 块体公式：E = ρ_a C_E |V| (q_sat(T_s) − q_a) · S_type  
  |V| 近地风；q_a≈q（单层）  
  S_type：地表类型调制（海=1、陆<1、冰≪1；见 7 章）  
- LH = L_v · E（进入地表能量负号项；docs/06）

凝结/降水 P_cond（kg m⁻² s⁻¹）与潜热释放 LH_release
- 当 q > q_sat(T_a)：P_cond = (q − q_sat)/τ_cond，q←q − P_cond·dt  
- LH_release = L_v · P_cond（加热大气；docs/06）

饱和比湿 q_sat(T)
- 采用 Tetens/Clausius–Clapeyron 近似；需连续、数值稳定


## 8.3 云量参数化与动态反照率（P003 融合）

动力–降水诊断（与湿度法融合）
- 基线：由风场散度/辐合诊断降水倾向（平滑阈值/连续函数，见 docs/03）  
- 一致性增强：将 P_cond 视为“物理上限/增益”，与散度法融合提升稳健性

云量融合与光学厚
- C_base = C_max · tanh(P/P_ref)（降水驱动），对云场做轻度高斯平滑  
- 背景源：由温度、涡度、锋生代理构成的 S_cloud  
- 记忆–降水–背景融合：Cⁿ⁺¹ = w_mem·Cⁿ + w_P·C_base + w_src·clip(Cⁿ + S_cloud·Δt/6h, 0,1)  
- 湿度–云一致性（增强项，可控）：  
  - RH 增益：RH>RH0 时按 K_Q 增强  
  - 凝结增益：按 P_cond（归一化到 PCOND_REF）以 K_P 增强  
- 反照率合成：α_total = α_base(type, Ts)·(1−C) + α_cloud·C（docs/06）  
- 长波发射：ε_eff = ε0 + k_lw_cloud·C（docs/06）


## 8.4 模块与接口（pygcm/humidity.py 及主循环）

核心函数（示意）
- q_init(Ts, RH0, p0) → q0  
- q_sat(T) → qsat  
- evaporation_flux(Ts, q, |V|, surface_factor, params) → E  
- condensation(q, T_a, dt, params) → (P_cond, q_next)

主循环次序（建议）
1) 云/降水诊断（动力 + 湿度一致性增强）→ C  
2) 反照率与辐射（docs/06）→ SW/LW；动力步  
3) 湿度步：E（→LH）、P_cond（→LH_release）、更新 q  
4) 地表/大气能量步：扣除 LH、加入 LH_release；海冰相变/海洋 Q_net（docs/07）  
5) 水文步：E–P–R 闭合（docs/09）


## 8.5 环境变量（节选；详见 docs/04-runtime-config.md）

湿度与蒸发
- QD_CE（1.3e-3）、QD_LV（2.5e6 J/kg）、QD_Q_INIT_RH（0.5）、QD_MBL_H（800 m 或 E→dq 因子）
- 地表类型蒸发缩放：QD_OCEAN_EVAP_SCALE（1.0）、QD_LAND_EVAP_SCALE（0.2）、QD_ICE_EVAP_SCALE（0.05）
- QD_TAU_COND（~1800 s）、QD_HUMIDITY_DIAG（1 打印诊断）、QD_Q_DIFF（1e-6..1e-5）

云–湿度一致性（P003/P008 融合）
- QD_CLOUD_COUPLE（1）、QD_RH0（0.6）、QD_K_Q（0.3）、QD_K_P（0.4）、QD_PCOND_REF（中位数）
- 云上限/尺度：QD_CMAX（0.95）、QD_PREF（正降水中位数）

反照率（与 docs/06 同步）
- QD_USE_TOPO_ALBEDO、QD_ALPHA_WATER、QD_ALPHA_ICE、QD_TRUECOLOR_CLOUD_*


## 8.6 推荐默认与运行示例

最小一致闭环（能量 + 湿度）
```bash
export QD_ENERGY_W=1
export QD_ENERGY_DIAG=1
export QD_HUMIDITY_DIAG=1

export QD_CE=1.3e-3
export QD_LV=2.5e6
export QD_Q_INIT_RH=0.5
export QD_TAU_COND=1800

# 海/陆/冰蒸发缩放
export QD_OCEAN_EVAP_SCALE=1.0
export QD_LAND_EVAP_SCALE=0.2
export QD_ICE_EVAP_SCALE=0.05

python3 -m scripts.run_simulation
```

云–湿度耦合增强（示例）
```bash
export QD_CLOUD_COUPLE=1
export QD_RH0=0.6
export QD_K_Q=0.3
export QD_K_P=0.4
# PCOND_REF 默认使用正 P_cond 的中位数
python3 -m scripts.run_simulation
```


## 8.7 诊断与验收标准（建议）

- 能量一致性：长期平均 ⟨LH⟩ ≈ ⟨LH_release⟩（差值很小）  
- 水量闭合：配合 docs/09，长期 ⟨E⟩ ≈ ⟨P⟩ + ⟨R⟩  
- 云–辐射：云量上升 → 地表短波吸收下降、长波温室增强；冷区/高地云更常见  
- 数值稳健：开启 P010“combo”反噪后，降水/云条纹显著减弱；q 场平滑与结构合理  
- 输出：HumidityDiag 打印 ⟨E⟩、⟨P_cond⟩、⟨LH⟩、⟨LH_release⟩；TrueColor/状态图可辨别云与降水结构


## 8.8 与其它模块的关系

- docs/06：LH 与 LH_release 的能量一致；SW/LW 的云光学厚/发射率增益  
- docs/07：SST 影响 q_sat(T_s) 与 E；海冰显著抑制蒸发  
- docs/09：E 与 P_cond 进入水量收支与水库变化；R 影响陆面与海洋闭环  
- docs/10：反噪抑制网格尺度噪音，改善降水与云的图像与诊断稳定性


## 8.9 变更记录（Changelog）

- 2025‑09‑20：完成 q 场 + E/LH + P_cond/LH_release；诊断开启  
- 2025‑09‑21：云–湿度耦合增强（RH、P_cond 增益）；与能量/海洋/水文一致性联测  
- 2025‑09‑21：文档迁移与整合至 docs/08‑humidity‑and‑clouds.md；与 04/06/07/09/10 交叉引用对齐
