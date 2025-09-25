# 12. 开发者指南：代码架构与模块 API（Developer Guide: Code Architecture & APIs）

本章迁移与扩展项目 P002（核心物理引擎实现）的“软件包结构与模块 API”方案，并结合当前仓库代码现状，形成对 `PyGCM for Qingdai` 的开发者导向文档。内容涵盖包结构、关键模块职责与接口、脚本入口、环境变量与扩展点，便于后续开发与维护。

关联文档
- docs/01–03：天文设定、轨道与气候框架
- docs/04：运行配置与环境变量目录
- docs/05–11：地形/反照率、能量、海洋/海冰、湿度/云、水文闭合、数值稳定、Spin-up；P014：地表路由与湖泊


## 12.1 包结构与职责

仓库关键路径（简要）
```
pygcm/
  ├── __init__.py
  ├── constants.py      # 物理/天文常数与星体参数
  ├── grid.py           # 球面网格与几何/度量
  ├── orbital.py        # 双星-行星系统几何与入射
  ├── forcing.py        # 天文入射 → 局地辐射平衡强迫（接口）
  ├── dynamics.py       # 大气单层浅水/原始方程步进、反噪、能量路径对接
  ├── energy.py         # 短/长波、边界层通量、地表/大气能量积分、Autotune
  ├── physics.py        # 反照率/云-辐射参数化拼装等（物理组合器）
  ├── topography.py     # 地形 I/O（NetCDF）、插值与 base_albedo/friction 生成
  ├── humidity.py       # q 场、蒸发 E、饱和与凝结、潜热释放
  ├── hydrology.py      # 水库/相态/陆地桶/径流/闭合诊断
  ├── routing.py        # 地表水文路由（P014）：河网/湖泊、异步汇流
  └── ocean.py          # 风驱动浅水海洋与 SST 平流、极点处理/稳定化
scripts/
  ├── run_simulation.py     # 主仿真入口（读取 env、初始化模块、主循环）
  ├── generate_topography.py # P004 地形生成
  ├── plot_topography.py     # 地形与属性可视化
  ├── diag_isr.py            # 双星短波分量诊断
  ├── generate_orbit_plots.py|test_orbital_module.py|verify_calculation.py
  └── generate_hydrology_maps.py  # P014 离线预计算（Pit Filling/流向/湖泊识别）
```

设计总览
- “核心环路”按 docs/06/07/08/09/10 的顺序协调：轨道几何/辐射 → 动力步 → 湿度步 → 能量步（地表/大气）→ 海洋步（若先平流再能量亦可）→ 水文闭合（P009）→（到达水文步长时）地表路由/湖泊（P014）→ 诊断/出图
- “数值稳定”在 dynamics/ocean 内通过 ∇⁴/滤波实现，参数集中由环境变量控制（详见 docs/10）
- 跨模块公共度量（如行星半径、自转、天文常量）集中在 constants.py


## 12.2 关键模块 API 摘要

constants.py
- Planet/system constants（示例）
  - R_PLANET, OMEGA, G, SIGMA_SB, M_SUN, L_SUN, AU, …
  - 双星参数（M_A/M_B、L_A/L_B、a_bin）、行星轨道（a_p）

grid.py
- class SphericalGrid(n_lat: int, n_lon: int)
  - lat, lon: 1D 坐标（deg）
  - lat_mesh, lon_mesh: 2D 网格
  - dlat_rad, dlon_rad: 弧度格距
  - metric/cache：cosφ、面积权重 w=cosφ 等
  - coriolis_param(Ω): f=2Ω sinφ（或作为属性预存）

orbital.py
- class OrbitalSystem(params)
  - positions(t) → (x_A,y_A),(x_B,y_B),(x_p,y_p)
  - distances(t) → d_A(t), d_B(t)
  - stellar_flux_total(t) → S_total(t)=L_A/(4πd_A²)+L_B/(4πd_B²)
- 备注：圆轨道/共面假设；用于 forcing 计算层顶/次太阳点几何（docs/02）

forcing.py
- class ThermalForcing(grid, orbital_system, planet_params)
  - calculate_insolation(t) → I(lat,lon,t)
  - calculate_equilibrium_temp(t, α) → T_eq(lat,lon,t)
  - 注：在能量框架启用时，T_eq 可仅作为诊断/对照；主能量路径见 energy.py

physics.py（物理组合器）
- calculate_dynamic_albedo(cloud, T_s, alpha_base, alpha_ice, alpha_cloud) → α_total
  - alpha_base: 标量或 2D 数组（base_albedo_map），优先使用地图
- 其它参数化拼装函数（如云-辐射增益）作为组合层，维持模块间一致性

energy.py
- shortwave_radiation(I, α_total, C, params) → (SW_atm, SW_sfc)
- longwave_radiation(Ts, Ta, C, params) → (LW_atm, LW_sfc)
- boundary_layer_fluxes(Ts, Ta, u, v, land_mask, params) → SH
  - LH 来自 humidity.evaporation（E）
- integrate_surface_energy(Ts, SW_sfc, LW_sfc, SH, LH, dt, *, C_s_map, params) → Ts_next
- integrate_surface_energy_with_seaice(..., h_ice, …) → (Ts_next, h_ice_next, ice_diag)
- integrate_atmos_energy_height(..., SW_atm, LW_atm, SH, LH_release, dt, params) → state_next
- autotune_greenhouse_params(params, diag) → params’（参见 docs/06）

topography.py
- load_topography_from_netcdf(path, grid, regrid="auto") → (elevation, land_mask, base_albedo, friction)
  - 经度周期双线性插值；land_mask 最近邻；范围/NaN 校验
- generate_base_properties(mask, elevation) → (base_albedo, friction)

humidity.py
- q_init(grid, RH0, T_a0|Ts) → q0
- q_sat(T) → qsat
- evaporation_flux(Ts, q, |V|, surface_factor, params) → E
- condensation(q, T_a, dt, params) → (P_cond, q_next)

hydrology.py
- partition_precip_phase(P, T_s, T_thresh) → (P_rain, P_snow)
- snow_step(S_snow, P_snow, melt_rate, dt) → (S_next, M_snow)
- update_land_bucket(W_land, P_land, E_land, dt, τ_runoff, W_cap=None) → (W_next, R)
- diagnose_water_closure(q, E, P, R, h_ice, W_land, S_snow, grid) → dict

ocean.py
- class WindDrivenSlabOcean(grid, land_mask, H_m, params)
  - step(dt, u_atm, v_atm, Q_net, ice_mask=None) → updates uo, vo, η, Ts
  - diagnostics() → {KE, umax, η_stats, CFL, …}
  - 内部：相对风应力、底摩擦、∇⁴/滤波、极点一致化、极区 sponge、SST 平流与热通量注入

routing.py
- class RiverRouting(grid, network_nc_path, dt_hydro_hours=6.0, …)
  - step(R_land_flux, dt_seconds, precip_flux=None, evap_flux=None) → 异步累计并在水文步长时执行汇流/湖泊更新
  - diagnostics() → {flow_accum_kgps, lake_volume_kg|map, ocean_inflow_kgps, mass_closure_error_kg}
  - 数据：flow_to_index/flow_order/lake_mask/lake_id 等由 data/hydrology_network.nc 提供（见 P014）

dynamics.py
- class SpectralModel/Model(grid, friction_map, params, …)
  - time_step(Teq|energy_forcing, dt) → 更新 u, v, h, Ts, C, q 等
  - 位置：对接 energy/humidity/hydrology/ocean 的调用；施加 ∇⁴/滤波；地转钳制/赤道正则化等
  - 具体类名/实现以仓库当前代码为准（能量路径启用时依赖 energy.integrate_*）

scripts/run_simulation.py（主入口）
- 读取环境变量（docs/04），初始化网格/模块/数据源；主循环内依序调用模块
- 打印诊断（EnergyDiag/HumidityDiag/WaterDiag/OceanDiag），按 QD_PLOT_EVERY_DAYS 出图


## 12.3 运行脚本与环境变量（与 docs/04 对齐）

- 运行控制：QD_DT_SECONDS（秒）、QD_TOTAL_YEARS（行星年）|QD_SIM_DAYS（行星日）、QD_PLOT_EVERY_DAYS（行星日）
- 地形与地表：QD_TOPO_NC、QD_USE_TOPO_ALBEDO、QD_OROG、QD_OROG_K
- 能量：QD_ENERGY_W、QD_ENERGY_DIAG、QD_T_FLOOR、QD_SW_*、QD_LW_*、QD_CH、QD_BOWEN_*
- 海洋/海冰：QD_USE_OCEAN、QD_OCEAN_USE_QNET、QD_MLD_M|QD_OCEAN_H_M、QD_CD、QD_R_BOT、QD_KH_OCEAN、QD_SIGMA4_OCEAN、QD_OCEAN_*（见 docs/07）
- 湿度：QD_CE、QD_LV、QD_Q_INIT_RH、QD_TAU_COND、QD_OCEAN/ICE/LAND_EVAP_SCALE、QD_HUMIDITY_DIAG
- 水文：QD_RUNOFF_TAU_DAYS、QD_WLAND_CAP、QD_SNOW_THRESH、QD_SNOW_MELT_RATE、QD_WATER_DIAG
- 反噪：QD_FILTER_TYPE、QD_SIGMA4、QD_K4_*、QD_DIFF_EVERY、QD_SHAPIRO_*、QD_SPEC_*、QD_DIFF_FACTOR


## 12.4 主循环建议顺序（一致性与稳定）

1) 轨道/入射几何与反照率（I、α_total；见 docs/05/06）
2) 动力步（u,v,h 等、地表摩擦、∇⁴/滤波；见 docs/10）
3) 湿度步：E（→LH）、P_cond（→LH_release）、q 的更新（见 docs/08）
4) 地表/大气能量：SW/LW/SH/LH 与海冰相变；Q_net 注入海洋（见 docs/06/07）
5) 海洋步：SST 平流与风应力驱动、极点一致化（见 docs/07）
6) 水文闭合：相态/雪/陆面桶/径流，WaterDiag（见 docs/09）
7) 可视化与日志：统一色表与指标；按需保存 NetCDF/图像


## 12.5 开发与扩展建议

- 参数管理
  - 统一从环境变量加载（docs/04），在 scripts/run_simulation.py 集中打印参数快照与拓扑统计
- 可测试性
  - 各模块提供纯函数/轻类接口，便于单元测试（维度/范围/守恒）
  - 诊断统一：TOA/SFC/ATM、⟨LH⟩ vs ⟨LH_release⟩、E–P–R 与水库变化
- 数值稳定
  - 优先使用“combo”（∇⁴+Shapiro），必要时加谱阻尼；极区注意度量保护与 sponge
- 性能与可视化
  - 限制出图频率；采用面积权重进行全局均值/积分；TrueColor 统一云/冰渲染参数
- 未来扩展
  - 多层大气/海洋、盐度/密度、海冰动力、云光学厚度诊断更物理化、单调性平流方案、重启/元数据标准化导出等


## 12.6 快速检查清单（Developers Checklist）

- [ ] 能量诊断（TOA/SFC/ATM）长期平均 |净| < 阈值
- [ ] ⟨LH⟩ ≈ ⟨LH_release⟩；⟨E⟩ ≈ ⟨P⟩ + ⟨R⟩
- [ ] 无“羽状条纹”；高波数方差比显著下降；大尺度动能保持
- [ ] 极圈一致化与 sponge 启用；SST/流速无经向条带
- [ ] 日志包含 topo 来源、海陆比、反照率/摩擦统计、参数快照

## 12.7 生态模块（P015 M1）：子包与接口

包结构补充
```
pygcm/
  └── ecology/
      ├── __init__.py          # 统一导出 SpectralBands / EcologyAdapter 等
      ├── spectral.py          # 光谱带定义与 TOA→surface 带强度 I_b 计算
      ├── types.py             # WeatherInstant/WeatherDaily 等数据结构
      ├── plant.py             # M1 极简 Plant（固定吸收模板，可后续替换为 Genes/FSM）
      ├── population.py        # M1 PopulationManager（聚合带反照率、缓存/子采样）
      └── adapter.py           # EcologyAdapter：连接主循环与生态子步/日步
```

核心 API（M1）
- SpectralBands / make_bands(nbands, λ0, λ1)：构造带边界与权重
- toa_to_surface_bands(I_total, cloud, bands, mode) → I_b[NB]：按 simple/rayleigh 模式将 gcm.isr 分配到各带
- class PopulationManager.step_subdaily(weather_inst, dt) → Optional[A_b_surface[NB]]：小时级聚合带反照率，按缓存策略低频重算
- class PopulationManager.step_daily(weather_day) → EcologyDailyReport：日级“慢路径”（M1 可简化）
- class EcologyAdapter：封装子步/日步的调度、缓存与回写

主循环接入（与 12.4 的顺序一致性）
1) 降水/云
2) 计算 insA/insB 与 gcm.isr（Forcing）
3) Ecology 子步（若开启且到采样频率）：
   - I_b ← toa_to_surface_bands(gcm.isr, cloud_eff, bands, mode)
   - A_b^surface ← PopulationManager.step_subdaily(...)
   - α_surface_ecology = Σ_b A_b^surface · (I_b / Σ I_b)
   - base_albedo_eff = 将 α_surface_ecology 融合到 base_albedo_map（陆面权重由 QD_ECO_LAI_ALBEDO_WEIGHT 控制）
4) albedo = calculate_dynamic_albedo(cloud, T_s, base_albedo_eff, α_ice, α_cloud, land_mask, ice_frac)
5) Teq 与 Dynamics/Ocean/Hydrology（与既有流程一致）
6) 达到日界时 EcologyAdapter.step_daily(weather_day)

运行控制（与 docs/04 第 11 节一致）
- QD_ECO_ENABLE=1、QD_ECO_SUBDAILY_ENABLE=1、QD_ECO_SUBSTEP_EVERY_NPHYS=1
- QD_ECO_FEEDBACK_MODE=instant、QD_ECO_ALBEDO_COUPLE=1
- QD_ECO_SPECTRAL_BANDS（M1 建议 16）、QD_ECO_TOA_TO_SURF_MODE=rayleigh（可选）
