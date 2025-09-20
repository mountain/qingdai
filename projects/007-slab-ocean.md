# Project 007: 平板海洋（Slab Ocean）与海冰（Sea Ice）混合层

本项目将为青黛 GCM 引入“平板海洋 / 混合层海洋（Slab / Mixed-Layer Ocean, MLO）”与最小海冰热力学，使地表具备真实的热惯性与冰-反照率反馈，并与 P006（能量收支）、P008（湿度/蒸发/凝结）、P009（水循环闭合）形成一致的物理框架。

## 0. 状态（2025-09-20）
- [x] M1：混合层海洋的等效热容量（C_s）地图与地表类型（海/陆）分离的能量积分
      - 已实现：在 energy.py 新增 integrate_surface_energy_map 支持逐网格 C_s；在 dynamics.SpectralModel 接入 C_s_map；
        在 scripts/run_simulation.py 基于 land_mask 构建 C_s_map 并传入（海洋 C_s = ρ_w·c_p,w·H_mld，陆地 C_s = QD_CS_LAND）。
        环境变量：QD_MLD_M、QD_CS_LAND、QD_RHO_W（默认1000）、QD_CP_W（默认4200）。
- [x] M2：最小海冰热力学（冻结/融化能量收支，α_ice 切换，低热容量）
      - 已实现：在 energy.py 新增 integrate_surface_energy_with_seaice（相变优先：先融后冻，能量闭合到 h_ice）；在 dynamics.SpectralModel 引入 h_ice 预报并调用该积分器；反照率使用基于厚度的 ice_frac = 1 − exp(−h_ice/H_ref)（scripts/run_simulation.py）。
      - 参数与默认：QD_USE_SEAICE=1、QD_T_FREEZE=271.35 K、QD_RHO_ICE=917、QD_LF=3.34e5 J/kg、QD_CS_ICE=5e6 J/m²/K、QD_HICE_REF=0.5 m。
- [x] M3：与辐射（反照率/长波）、边界层通量（SH/LH）、蒸发/凝结的联动
- [ ] M4：诊断与验证（延迟效应、季节性相位、能量平衡核算）
- [ ] M5：参数标定（H_mld、α_water/α_ice、冻结点、冰导热/热容等）

### 实现状态更新（2025-09-20）
- 已接通反照率—海冰—云：scripts/run_simulation.py 使用 h_ice → ice_frac 融合到 calculate_dynamic_albedo；cloud_eff_last 驱动辐射/反照率。
- 已接通边界层通量：pygcm/dynamics.SpectralModel 在能量路径调用 boundary_layer_fluxes 计算 SH，并在海冰/海洋/陆地的地表能量积分中使用；LH 由湿度模块蒸发得到并进入地表能量。
- 已接通湿度闭环：pygcm/humidity.py 提供 E 与 P_cond，dynamics.SpectralModel 将 LH = L_v·E 扣除地表能量、LH_release = L_v·P_cond 注入大气能量。设置 QD_ENERGY_W > 0 时，SW_atm/LW_atm/SH/LH_release 经 integrate_atmos_energy_height 加入大气能量路径。
- 运行开关：建议设置 QD_ENERGY_W=1、QD_USE_SEAICE=1、QD_ENERGY_DIAG=1 以启用完整 M3 联动与诊断。

关联项目：
- P006 动态能量收支框架（已完成 M1）
- P008 大气湿度 q 的引入（蒸发/凝结、潜热一致）
- P009 水循环闭合（蒸发→凝结/降水→海冰/融化）

## 1. 目标

- 引入混合层海洋：以固定深度 H_mld 的无流海洋层表征海洋热容量，形成显著的昼夜与季节滞后。
- 实现最小海冰方案：当海面温度 T_s 接近/低于冻结点时，将净能量用于相变而非继续降温；冰面具有高反照率与低热容量。
- 与 P006 的能量收支一致：短波/长波/SH/LH 的地表能量方程在“海/冰/陆”下分别使用相应参数。
- 为 P008 的蒸发源项提供地表控制：开阔海面强蒸发，结冰区域显著减弱（或近零）。
- 在 P009 水量收支中记账：E-P 的海—气交换闭合，冻结与融化的相变能量与质量一致。

## 2. 模型设计与公式

### 2.1 混合层海洋的等效热容量
- 海洋混合层热容量（每单位面积）：
  [C_s]_ocean = ρ_w · c_p,w · H_mld
  - ρ_w ≈ 1000 kg m⁻³；c_p,w ≈ 4200 J kg⁻¹ K⁻¹；H_mld 典型 30–70 m（默认 50 m）
  - 对于陆地，[C_s]_land 取较小值（默认 2e6–5e6 J m⁻² K⁻¹），体现快速响应
- 地表热容量地图：
  C_s_map = where(land_mask==0, [C_s]_ocean, [C_s]_land)

### 2.2 地表能量方程与类型分离
- 继承 P006 表达式：
  C_s dT_s/dt = SW_sfc − LW_sfc − SH − LH
- 对海/陆/冰分别采用不同参数：
  - 反照率：α_water（~0.06–0.10）、α_land（来自 base_albedo_map）、α_ice（~0.5–0.7）
  - 发射率/温室：通过云与大气参数（P006），不在此重复定义
  - C_s：海/陆/冰不同（冰层取低热容）
- 实现上优先采用“每格点 C_s_map”推进 T_s（见 §4 接口建议）。

### 2.3 最小海冰热力学
- 冻结点 T_freeze ≈ 271.35 K（-1.8℃，盐水近似，可配置）
- 融化潜热 L_f ≈ 3.34×10⁵ J kg⁻¹（淡水近似；盐水可适度调整）
- 冰密度 ρ_i ≈ 917 kg m⁻³
- 能量分配策略（单层近似）：
  1) 当地表为开阔海（无冰）且 T_s → T_freeze 且净通量 Q_net < 0：
     - 优先将能量亏损用于“生成冰厚”：Δh_i = −(Q_net · dt) / (ρ_i · L_f)
     - T_s 不低于 T_freeze（或设置窄带缓冲）
  2) 当存在海冰且净通量 Q_net > 0 且 T_s,ice → 273.15 K：
     - 优先用于融化冰：Δh_i = −(Q_net · dt) / (ρ_i · L_f)，不为负
  3) 冰面热容量：采用小 C_s,ice（雪/冰层浅表热容），使冰面对辐射快速响应
  4) 反照率切换：有冰区域使用 α_ice，显著提高反照率，触发冰-反照率反馈
- 细化（后续 M5）：导热与冰下海水耦合、积雪与多层结构暂不引入

## 3. 接口与数据结构（建议）

### 3.1 energy.py 扩展（与 P006 对齐）
- 新增（建议）：
  - integrate_surface_energy_map(Ts, SW_sfc, LW_sfc, SH, LH, dt, C_s_map, params)
  - surface_type_map: {0:ocean_open, 1:land, 2:sea_ice}（或分离 boolean masks）
  - sea_ice_thickness（可选诊断/预报数组），在能量积分前后更新
- 兼容保留：现有 integrate_surface_energy(Ts,...) 支持标量 C_s；map 版本优先

### 3.2 反照率与辐射
- P006 中 calculate_dynamic_albedo 已接受 base_albedo_map；对海/冰：
  - 海：使用 α_water 基底；冰：使用 α_ice 基底
  - 与云量 α_cloud 的混合保持一致
- 与 P005（外部地表）协作：当存在 base_albedo_map 时，海/陆/冰的融合优先基于实际类型

### 3.3 湿度与蒸发（与 P008）
- 开阔海面蒸发率大，E = ρ_a C_E |V| (q_sat(T_s) − q_a)
- 海冰/陆地处 E 显著减弱（可近似 0–微量），通过 surface_type_map 控制
- 产出的 E 同时作为地表潜热（LH）与大气水汽源项（见 P008）

## 4. 任务拆解

- T1（M1）C_s 地图
  - 从 land_mask 与参数 H_mld 构建 C_s_map，接入地表能量推进
  - energy.integrate_surface_energy_map（或在 dynamics 中就地实现 map 版推进）
- T2（M2）海冰最小方案
  - 引入 sea_ice_mask 与 h_ice（诊断/预报）
  - 在能量积分中分配相变能量并更新 h_ice、surface_type、α
- T3（M3）与辐射/BL/蒸发联动
  - 动态 α_water/α_ice 混合进入短波
  - SH/LH：海/陆/冰差异化；E 控制 LH 与 q 源项（对接 P008）
- T4（M4）诊断与验证
  - 延迟效应：T_s 对 ISR 的相位滞后、季节性幅度
  - 能量收支：TOA/SFC/ATM 全局均衡的长期平均
- T5（M5）参数标定
  - 参数扫描：H_mld、α_ice、T_freeze、C_E 等；选择稳定默认组

## 5. 环境变量（建议默认）
- QD_MLD_M（默认 50）：海洋混合层深度（m）
- QD_CS_LAND（默认 3e6）：陆地等效热容量（J m⁻² K⁻¹）
- QD_CS_ICE（默认 5e6）：海冰/积雪等效热容量（J m⁻² K⁻¹，薄冰面层）
- QD_RHO_W（默认 1000）：海水密度（kg m⁻³）用于 C_s_ocean
- QD_CP_W（默认 4200）：海水比热（J kg⁻¹ K⁻¹）用于 C_s_ocean
- QD_ALPHA_WATER（默认 0.08）、QD_ALPHA_ICE（默认 0.60）
- QD_T_FREEZE（默认 271.35 K）、QD_RHO_ICE（默认 917）、QD_LF（默认 3.34e5）
- QD_HICE_REF（默认 0.5 m）：将冰厚转换为光学冰覆盖的 e 折算厚度
- QD_USE_SEAICE（默认 1）：启用海冰方案
- QD_OCEAN_EVAP_SCALE（默认 1.0）、QD_ICE_EVAP_SCALE（默认 0.05）：蒸发缩放（与 P008 配合）

## 6. 诊断与输出
- 海洋/陆地/海冰面积分随时间演变
- h_ice 全局均值/分布、α_total 统计
- 能量收支（与 P006 使用相同的 compute_energy_diagnostics）
- 蒸发 E 与 E−P 场（与 P008/P009 整合）

## 7. 测试与验收
- [ ] 昼夜/季节滞后相位合理（与 H_mld 调整相关）
- [ ] 冰-反照率反馈可触发但不数值爆裂（参数温和）
- [ ] 长期平均 TOA/SFC/ATM 接近守恒（与 P006-M4 共同验证）
- [ ] 海冰生成/融化对 α_total 与 E 的影响方向正确

## 8. 与其他项目的交叉引用
- P006（能量收支）：本项目提供 C_s_map、surface_type 与 α 切换；能量积分中处理相变能。
- P008（湿度 q）：开阔海面的 E 为 q 源项；凝结释放的潜热回到大气；海冰/陆地显著抑制 E。
- P009（水循环）：E−P 的闭合核算；降水回到海洋/陆地；冻结/融化在能量/质量上自洽。

## 9. 运行示例（占位）
```bash
# 启用能量收支 + slab ocean + 海冰最小方案（参数温和）
export QD_ENERGY_W=1
export QD_ENERGY_DIAG=1
export QD_MLD_M=50
export QD_ALPHA_WATER=0.08
export QD_ALPHA_ICE=0.60
export QD_T_FREEZE=271.35
export QD_USE_SEAICE=1

python3 -m scripts.run_simulation
```

---

附注：本项目侧重“热容量 + 相变 + 反照率反馈”的一阶效应，不引入洋流与海冰动力。随着 P006/P008/P009 的推进，海气耦合将具备完整的能量与水量一致性。
