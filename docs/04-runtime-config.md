# 4. 运行配置与环境变量目录（Runtime Config & Env Catalogue）

本文梳理当前实现所涉及的运行时环境变量（env），覆盖气候模型的主要子系统，并给出推荐的使用场景与调参指导。所有变量均可在运行前通过环境变量导出，或在脚本中设置（例如 `spin-up.sh` 已提供常用组合）。

- 适用版本：本仓库当前实现（含 P004–P014 以及最近动力学/可视化修订）
- 单位：若未特别说明，SI 制；角度单位见各条目

目录：
1) 全局与运行控制  
2) 地形与地表参数（P004/P005）  
3) 能量收支与辐射/边界层（P006）  
4) 平板海洋/海冰与动态海洋（P007/P011）  
5) 湿度 q 与蒸发—凝结—潜热（P008）  
6) 水循环闭合（P009）  
7) 动力学反噪与滤波（P010）  
8) 可视化与诊断（绘图/日志）  
9) 自旋与重启（P013）
10) 路由与湖泊（P014）

---

## 1) 全局与运行控制

- QD_DT_SECONDS（默认 300）积分步长，秒  
- QD_SIM_DAYS（可选）整段运行时长，单位：行星日（不设置则走 QD_TOTAL_YEARS 或默认 5 年）  
- QD_TOTAL_YEARS（推荐）整段运行时长，单位：行星年；优先于 QD_SIM_DAYS  
- QD_PLOT_EVERY_DAYS（默认 10）出图间隔（行星日），减少 I/O 可适当调大  
- QD_TS_MIN / QD_TS_MAX（默认 150/340 K）表面温度钳制（可视化/数值保护）

动力学（浅水动量）选择：
- QD_MOM_SCHEME（默认 geos）  
  - geos：地转松弛方案，数值稳健  
  - primitive：显式原始动量方程 du/dt/dv/dt（含压强梯度、科氏力、线性摩擦）

---

## 2) 地形与地表参数（P004/P005）

- QD_TOPO_NC：外部 NetCDF 地形路径（含 elevation / land_mask / base_albedo / friction），若未设置则使用程序生成的简化地形  
- QD_USE_TOPO_ALBEDO（默认 1）：使用外部 `base_albedo` 融合云/冰的动态反照率  
- QD_OROG（默认 0）：启用地形降水增强（迎风抬升）  
- QD_OROG_K（默认 7e-4）：地形抬升强度系数（单位量纲与风/坡度近似匹配，建议微调）

---

## 3) 能量收支与辐射/边界层（P006）

主开关：
- QD_ENERGY_W（0..1，默认 0）能量收支权重；=1 时完全使用显式能量路径  
- QD_ENERGY_DIAG（默认 1）定期打印 TOA/SFC/ATM 收支诊断  
- QD_T_FLOOR（默认 150 K）夜侧/高原温度下限（地表能量积分保护）
- QD_GH_LOCK（默认 1）固定温室效应因子 g，按 g = 1 − OLR/(σ Ts^4) 强制长波出流与下行长波（见 energy.py）  
- QD_GH_FACTOR（默认 0.40）固定的温室因子 g（地球量级 ≈ 0.4）；当 QD_GH_LOCK=1 时自动禁用 QD_ENERGY_AUTOTUNE，且动力学旧路径的 greenhouse_factor 亦从该值读取

短波与长波（简化参数化，实际实现见 `pygcm/energy.py`）：
- QD_SW_A0（默认 0.06）大气短波基吸收  
- QD_SW_KC（默认 0.20）云短波吸收增强  
- QD_LW_EPS0（默认 0.70）无云大气发射率  
- QD_LW_KC（默认 0.20）云长波增强  
- QD_LW_V2（默认 1）启用 v2 长波（含地表发射率图）

边界层与鲍文比：
- QD_CH（默认 1.5e-3）感热交换系数  
- QD_CP_A（默认 1004.0）大气比热  
- QD_BOWEN_LAND（默认 0.7）、QD_BOWEN_OCEAN（默认 0.3）海陆差异鲍文比

大气能量耦合与自调参：
- QD_ATM_H（默认 800 m）单层大气等效厚度，用于将能量源项转为高度倾向  
- QD_ENERGY_AUTOTUNE（默认 0）温室参数自动微调  
- QD_ENERGY_TUNE_EVERY（默认 50 步）自调参频率

---

## 4) 平板海洋/海冰与动态海洋（P007/P011）

平板海洋/海冰（地表热容量地图、相变）：
- QD_MLD_M（默认 50 m）混合层深度（决定海洋 C_s）  
- QD_CS_LAND（默认 3e6 J/m^2/K）、QD_CS_ICE（默认 5e6）陆/冰等效热容量  
- QD_ALPHA_WATER（默认 0.08）、QD_ALPHA_ICE（默认 0.60）海/冰反照率基底  
- QD_T_FREEZE（默认 271.35 K）冻结点  
- QD_RHO_ICE（默认 917 kg/m^3）、QD_LF（默认 3.34e5 J/kg）冰的密度/潜热  
- QD_HICE_REF（默认 0.5 m）将冰厚转换为光学冰覆盖的 e 折算厚度  
- QD_USE_SEAICE（默认 1）启用海冰最小热力学

动态海洋（WindDrivenSlabOcean）：
- QD_USE_OCEAN（默认 1）启用动态海洋  
- QD_OCEAN_H_M（默认等于 QD_MLD_M）海洋动力学使用的深度  
- QD_CD（默认 1.5e-3）风应力拖曳系数  
- QD_R_BOT（默认 2.0e-5 s^-1）海洋底摩擦  
- QD_RHO_A（默认 1.2）风应力中的空气密度  
- QD_WIND_STRESS_VCAP（默认 15 m/s）风应力速度上限  
- QD_TAU_SCALE（默认 0.2）风应力效率（向浅层动量传递的比例）  
- QD_POLAR_SPONGE_LAT（默认 70°）、QD_POLAR_SPONGE_GAIN（默认 5e-5 s^-1）极区阻尼

海洋混合/反噪/数值保护：
- QD_KH_OCEAN（默认 5.0e3 m^2/s）SST 水平扩散  
- QD_SIGMA4_OCEAN（默认 0.02）∇⁴ 强度（自适应到 K4）  
- QD_OCEAN_K4_NSUB（默认 1）∇⁴ 子步  
- QD_OCEAN_DIFF_EVERY（默认 1）∇⁴ 执行频率  
- QD_OCEAN_SHAPIRO_N（默认 0 关闭）、QD_OCEAN_SHAPIRO_EVERY（默认 8）Shapiro 滤波  
- QD_OCEAN_CFL（默认 0.5）CFL 目标值，用于自动子步数  
- QD_OCEAN_MAX_U（默认 3 m/s）海流速度上限  
- QD_OCEAN_OUTLIER（默认 mean4，可选 clamp）异常海流处理  
- QD_ETA_CAP（默认 5 m）海表高度异常上限

海—气热通量耦合与诊断：
- QD_OCEAN_USE_QNET（默认 1）将 Q_net/(ρ c_p H) 注入 SST  
- QD_OCEAN_ICE_QFAC（默认 0.2）海冰下的垂直热通量比例（相对开阔海）  
- QD_OCEAN_ADV_ALPHA（默认 0.7）SST 半拉氏平流混合权  
- QD_OCEAN_ENERGY_DIAG（默认 1）、QD_OCEAN_DIAG_EVERY（默认 200 步）海洋能量诊断  
- QD_OCEAN_POLAR_LAT（默认 60°）极区带诊断范围

---

## 5) 湿度 q 与蒸发—凝结—潜热（P008）

- QD_CE（默认 1.3e-3）蒸发块体公式系数  
- QD_LV（默认 2.5e6 J/kg）汽化潜热  
- QD_Q_INIT_RH（默认 0.5）初始相对湿度  
- QD_TAU_COND（典型 1800 s）凝结时间尺度（若实现中暴露）  
- QD_MBL_H（默认 800 m）混合边界层厚度（质量换算）  
- QD_OCEAN_EVAP_SCALE（默认 1.0）开阔海蒸发缩放  
- QD_LAND_EVAP_SCALE（默认 0.2）陆地蒸发缩放  
- QD_ICE_EVAP_SCALE（默认 0.05）海冰上蒸发缩放  
- QD_HUMIDITY_DIAG（默认 1）湿度/潜热诊断  
- QD_CLOUD_COUPLE（默认 1）湿度/降水对云光学厚影响  
- QD_RH0（默认 0.6）、QD_K_Q（默认 0.3）相对湿度对云量增益  
- QD_K_P（默认 0.4）、QD_PCOND_REF（默认中位数）凝结对云量增益  
- QD_Q_DIFF（默认 1e-6–1e-5）q 的温和扩散强度（如实现中暴露）

---

## 6) 水循环闭合（P009）

- QD_RUNOFF_TAU_DAYS（默认 10 天）陆面“桶”径流时标  
- QD_WLAND_CAP（可选）桶容量（mm）  
- QD_SNOW_THRESH（默认 273.15 K）雨/雪阈值  
- QD_SNOW_MELT_RATE（默认 5 mm/day）融雪速率  
- QD_WATER_DIAG（默认 1）水量闭合诊断打印

---

## 7) 动力学反噪与滤波（P010）

主控：
- QD_DIFF_ENABLE（默认 1）启用数值抑噪  
- QD_FILTER_TYPE（默认 combo，可选 hyper4 | shapiro | spectral | combo）  
- QD_DIFF_EVERY（默认 1）施加频率  
- QD_DIFF_FACTOR（默认 0.998）温和全局扩散（乘法因子）

超扩散（∇⁴）：
- QD_SIGMA4（默认 0.02）以无量纲 σ₄ 计算 K₄=σ₄·Δx_min⁴/dt  
- QD_K4_U/V/H/Q/CLOUD（若直接给定系数以覆盖自适应）  
- QD_K4_NSUB（默认 1）子步  
- QD_DYN_DIAG（默认 0）打印反噪诊断

Shapiro 与谱带阻：
- QD_SHAPIRO_N（默认 2）阶数  
- QD_SHAPIRO_EVERY（默认 6）频率（步）  
- QD_SPEC_EVERY（默认 0 关闭）、QD_SPEC_CUTOFF（默认 0.75）、QD_SPEC_DAMP（默认 0.5）

---

## 8) 可视化与诊断（绘图/日志）

- QD_PLOT_EVERY_DAYS（默认 10）出图间隔（行星日）  
- QD_PLOT_ISR（默认 0）额外输出双星短波分量图  
- QD_PLOT_PS_MODE（默认 anom）表面气压绘图模式：  
  - anom：压强距平（hPa）= ρ g h / 100  
  - abs：绝对气压（hPa）= (p0 + ρ g h) / 100
- QD_TRUECOLOR_ICE_FRAC（默认 0.15）TrueColor 冰渲染阈值  
- QD_TRUECOLOR_CLOUD_ALPHA（默认 0.60）云不透明度  
- QD_TRUECOLOR_CLOUD_WHITE（默认 0.95）云白度  
- QD_TRUECOLOR_SNOW_BY_TS（默认 0）按温度渲染陆地积雪

能量/湿度/水文/海洋诊断：
- QD_ENERGY_DIAG、QD_HUMIDITY_DIAG、QD_WATER_DIAG、QD_OCEAN_ENERGY_DIAG 见各模块

---

## 9) 自旋与重启（P013）

- QD_RESTART_IN：重启输入（NetCDF，由 `save_restart` 写出字段）  
- QD_RESTART_OUT：重启输出路径（NetCDF）  
- QD_INIT_BANDED（默认 1 in spin-up.sh；代码默认 0）：分带初始 Ts  
- QD_INIT_T_EQ（默认 295 K）、QD_INIT_T_POLE（默认 265 K）分带初始温度端值  
- QD_TOTAL_YEARS / QD_SIM_DAYS：运行时长（建议使用年）  

脚本（两阶段 SOP）：  
- spin-up.sh  
  - PHASE1_YEARS / PHASE2_YEARS（默认 150 / 120）  
  - MLD_P1 / MLD_P2（默认 5 m / 50 m）  
  - RESTART_P1 / RESTART_EQ（默认 `restart_phase1.nc` / `restart_equilibrium.nc`）

---

## 10) 路由与湖泊（P014）

- QD_HYDRO_NETCDF（默认 data/hydrology_network.nc）：离线路由网络（flow_to_index/flow_order/lake_*）NetCDF  
- QD_HYDRO_DT_HOURS（默认 6）：水文步长（小时），达到该累计时长时执行一次全图路由/湖泊水量更新  
- QD_TREAT_LAKE_AS_WATER（默认 1）：湖面在能量/湿度路径上按水体（海洋）处理  
- QD_ALPHA_LAKE（可选）：覆盖湖面基础反照率（不设则与海洋相同）  
- QD_HYDRO_DIAG（默认 1）：打印路由诊断（入海通量、最大流量、质量闭合误差等）  

# 推荐配置与使用指引

场景 A：稳定默认（快速体验/出图）
- export QD_TOTAL_YEARS=0.05（或 QD_SIM_DAYS=50）  
- export QD_USE_OCEAN=1 QD_ENERGY_W=1 QD_ENERGY_DIAG=1  
- export QD_FILTER_TYPE=combo QD_SIGMA4=0.02 QD_SHAPIRO_EVERY=6 QD_SHAPIRO_N=2  
- python3 -m scripts.run_simulation

场景 B：两阶段 Spin-up（推荐用于科学试验）
- chmod +x spin-up.sh  
- ./spin-up.sh  
- 阶段一（MLD=5 m，150 年）→ 生成 restart_phase1.nc  
- 阶段二（MLD=50 m，120 年）→ 生成 restart_equilibrium.nc

场景 C：动力学检验（显式科氏偏转）
- export QD_MOM_SCHEME=primitive  
- export QD_TOTAL_YEARS=0.02  
- python3 -m scripts.run_simulation  
说明：primitive 方案显式包含 PGF 与科氏项，更便于检查风向相对等压线的偏转；数值更敏感，需保持 P010 的反噪默认。

场景 D：高抑噪（减弱条纹/高频伪迹）
- export QD_FILTER_TYPE=combo QD_SIGMA4=0.03 QD_K4_NSUB=2  
- 可叠加：export QD_SPEC_EVERY=6 QD_SPEC_CUTOFF=0.70 QD_SPEC_DAMP=0.5

场景 E：能量与水量闭合诊断
- export QD_ENERGY_W=1 QD_ENERGY_DIAG=1  
- export QD_HUMIDITY_DIAG=1 QD_WATER_DIAG=1  
- 观察 TOA/SFC/ATM 近守恒与 ⟨LH⟩≈⟨LH_release⟩、⟨E⟩≈⟨P⟩+⟨R⟩

---

# 变量相互作用与常见问题

- 表面气压负值：绘图默认显示“距平”（anom），若切换成绝对值需设 QD_PLOT_PS_MODE=abs；确保 h 的定义与公式匹配。  
- 海冰与反照率：`h_ice → ice_frac → α_total/LW emissivity` 路径已贯通；海冰下 `QD_OCEAN_ICE_QFAC` 允许弱耦合，避免极区热力孤立。  
- 动力学稳定性：primitive 方案对步长/反噪更敏感，建议保持 `QD_SIGMA4≈0.02–0.04`、`QD_SHAPIRO_EVERY≈6`。  
- 运行时长优先级：QD_TOTAL_YEARS > QD_SIM_DAYS > 默认（5 年）。Spin-up 建议用年。  
- 外部地形插值：`pygcm/topography.load_topography_from_netcdf` 会按经度周期插值、`land_mask` 用最近邻，注意分辨率差异。

---

# 参考

- P004/005/006/007/008/009/010/011/012/013 项目文档  
- 脚本：`scripts/run_simulation.py`、`spin-up.sh`  
- 主要模块：`pygcm/dynamics.py`、`pygcm/energy.py`、`pygcm/humidity.py`、`pygcm/hydrology.py`、`pygcm/ocean.py`、`pygcm/topography.py`、`pygcm/routing.py`

如需把本目录集成到 README 的“运行 GCM”小节，可在后续提交中将该文档链接加入目录列表。
