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

- QD_HYDRO_ENABLE（默认 1）：开启/关闭在线径流路由模块（RiverRouting）  
- QD_HYDRO_NETCDF（默认 data/hydrology_network.nc）：离线路由网络（flow_to_index/flow_order/lake_*）NetCDF  
- QD_HYDRO_DT_HOURS（默认 6）：水文步长（小时），达到该累计时长时执行一次全图路由/湖泊水量更新  
- QD_TREAT_LAKE_AS_WATER（默认 1）：湖面在能量/湿度路径上按水体（海洋）处理  
- QD_ALPHA_LAKE（可选）：覆盖湖面基础反照率（不设则与海洋相同）  
- QD_HYDRO_DIAG（默认 1）：打印路由诊断（入海通量、最大流量、质量闭合误差等）  

可视化（河网/湖泊叠加，M4）：
- QD_PLOT_RIVERS（默认 1）：在状态图与 TrueColor 上叠加河网与湖泊图层  
- QD_RIVER_MIN_KGPS（默认 1e6）：河网显示阈值（kg/s），仅显示主干河流（可按分辨率/需要调整）  
- QD_RIVER_ALPHA（默认 0.35；TrueColor 中 0.45）：河网叠加透明度  
- QD_LAKE_ALPHA（默认 0.40）：湖泊叠加透明度  

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

## 11) 生态模块（P015：Emergent Ecology & Spectral Dynamics）

实现注记（2025‑09‑25）：已接入 M1 小时级适配器（陆地标量 α 回写）；PopulationManager 与日级慢路径待 M2/M3。
双时序（时级/日级）生态接口，用于植被—辐射即时耦合与日级慢过程（形态投资、生命周期、繁殖/传播）。

主控与步长
- QD_ECO_ENABLE（默认 0）：开启/关闭生态模块（Plant/PopulationManager）。  
- QD_ECO_DT_DAYS（默认 1.0）：生态日级更新步长（天），用于慢路径。  
- QD_ECO_SUBDAILY_ENABLE（默认 1）：启用与物理步对齐的时级子步接口。  
- QD_ECO_SUBSTEP_EVERY_NPHYS（默认 1）：每 N 个物理步调用一次时级子步。  
- QD_ECO_DT_HOURS（默认 自动由 dt_seconds 换算）：用于时级诊断打印。  
- QD_ECO_FEEDBACK_MODE（instant|daily，默认 instant）：反照率写回短波的策略；即时/仅日末。  
- QD_ECO_ALBEDO_COUPLE_FREQ（subdaily|daily，默认 subdaily）：反照率写回频率。

光谱与带宽（与 docs/14 对齐）
- QD_ECO_SPECTRAL_BANDS（默认 8）：短波离散波段数（带分辨）。  
- QD_ECO_SPECTRAL_RANGE_NM（默认 380,780）：可见光范围（nm，下限,上限）。  
- QD_ECO_TOA_TO_SURF_MODE（simple|rayleigh|custom，默认 simple）：层顶→地表光谱调制模式。  
- QD_ECO_RAYLEIGH_T0（默认 0.9）：Rayleigh 模式基础透过率。  
- QD_ECO_RAYLEIGH_LREF_NM（默认 550）：Rayleigh 参照波长（nm）。  
- QD_ECO_RAYLEIGH_ETA（默认 4.0）：Rayleigh 衰减指数。  
- QD_ECO_SOIL_SPECTRUM（路径，可选）：土壤背景光谱文件。  
- QD_ECO_TRUECOLOR_ENABLE（默认 1）：启用基于 CIE 的 TrueColor 渲染（可视化）。  
- QD_ECO_TRUECOLOR_GAMMA（默认 2.2）：TrueColor gamma 校正。

光竞争（冠层）
- QD_ECO_LIGHT_K（默认 0.5）：冠层光衰减系数（Beer-Lambert 型）。  
- QD_ECO_LIGHT_UPDATE_EVERY_HOURS（默认 6）：时级重算冠层/反照率的最小间隔（小时）。  
- QD_ECO_LIGHT_RECOMPUTE_LAI_DELTA（默认 0.05）：LAI 相对变化阈值，超过强制重算。

水竞争（根系）
- QD_ECO_WATER_PRIORITY（root_mass|depth_weighted，默认 root_mass）：根系权重策略。  
- QD_ECO_SOIL_WATER_CAP（可选）：单日或小时可用水上限比例（归一化），超出截断。

形态/投资/物候（个体级，docs/13）
- QD_ECO_SEED_GERMINATE_GDD（默认 80）：发芽积温阈值。  
- QD_ECO_STRESS_WATER_DAYS（默认 7）：连续水分胁迫进入衰老阈值（天）。  
- QD_ECO_ALLOC_ROOT（可选）：覆盖基因的根投资比例（0..1）。  
- QD_ECO_ALLOC_STEM（可选）：覆盖基因的茎投资比例（0..1）。  
- QD_ECO_ALLOC_LEAF（可选）：覆盖基因的叶投资比例（0..1）。  
- QD_ECO_HEIGHT_EXPONENT（默认 0.8）：height ∝ stem_mass^γ 的 γ。  
- QD_ECO_REPRO_FRACTION（默认 0.2）：成熟期用于繁殖的能量比例。  
- QD_ECO_MAINT_COST（默认 0）：维护能量成本（若启用，从日能量中扣除）。  
- QD_ECO_MIN_LEAF_AREA（默认 0）：小于该叶面积阈值进入 SENESCENT/DEAD。

演化/传播（种群级，docs/15）
- QD_ECO_SEED_BANK_MAX（默认 1000）：本地种子库容量上限（溢出丢弃或劣化）。  
- QD_ECO_LONGDIST_FRAC（默认 0.05）：远距离传播比例（余下本地播种）。  
- QD_ECO_MUT_RATE（默认 1e-3）：突变概率（每颗种子）。  
- QD_ECO_INIT_SPECIES（默认 grass,tree）：初始化基因型集合，逗号分隔。

聚合与短波耦合
- QD_ECO_LAI_ALBEDO_WEIGHT（默认 1.0）：LAI/叶面积在反照率聚合中的权重。  
- QD_ECO_ALBEDO_COUPLE（默认 1）：开启生态反照率回写短波带 α。  

诊断与可视化
- QD_ECO_DIAG（默认 1）：生态诊断打印（EcologyDiag）与可选图层。  
- QD_ECO_OPEN（默认 1）：macOS 下首次绘图自动打开生态面板（ecology_day_*.png）。  
- QD_ECO_BANDS_COUPLE（默认 0）：启用带化反照率与短波耦合（日级降维混入）。  
- QD_ECO_TRUECOLOR_VEG（默认 1）：TrueColor 植被颜色叠加（按带反射与日/夜/I_b 调制）。  
- QD_ECO_VEG_MIX_EXP（默认 1.3）：TrueColor 叠加时对冠层因子 f(LAI) 的指数强化，提升高 LAI 时的“绿度主导感”。  

说明
- 时级接口：PopulationManager.step_subdaily 与 Plant.update_substep 仅累积“当日能量/小时胁迫”，不进行形态大跳转；反照率按缓存策略低频重算（配置见上）。  
- 日级接口：PopulationManager.step_daily 与 Plant.update_one_day 执行“慢路径”（形态投资、生命周期更替、繁殖/传播/突变）并生成日末诊断。  
- 反照率耦合：FEEDBACK_MODE=instant 时，时级返回的 A_b^surface 将立即写回短波，下一物理步生效；daily 时仅在日末写回。

# 参考

- P004/005/006/007/008/009/010/011/012/013 项目文档  
- 脚本：`scripts/run_simulation.py`、`spin-up.sh`  
- 主要模块：`pygcm/dynamics.py`、`pygcm/energy.py`、`pygcm/humidity.py`、`pygcm/hydrology.py`、`pygcm/ocean.py`、`pygcm/topography.py`、`pygcm/routing.py`

如需把本目录集成到 README 的“运行 GCM”小节，可在后续提交中将该文档链接加入目录列表。

### P015 M1 快速配置（小时级生态回耦）

当前代码实现了生态模块的小时级子步接口（Sub‑daily），可在每个物理步（dt 秒）按光谱带 I_b(t) 计算植被聚合反照率并即时回写到短波。推荐的最小配置如下（NB=16、每物理步子步、即时回耦）：

```bash
export QD_ECO_ENABLE=1
export QD_ECO_SUBDAILY_ENABLE=1
export QD_ECO_SUBSTEP_EVERY_NPHYS=1     # 每 N 个物理步调用 1 次子步，这里为每步
export QD_ECO_FEEDBACK_MODE=instant     # 子步立即回写用于下一物理步
export QD_ECO_ALBEDO_COUPLE=1           # 开启生态反照率回写
export QD_ECO_SPECTRAL_BANDS=16         # 光谱带数（建议 16）
# 可选：TOA→Surface 光谱调制
export QD_ECO_TOA_TO_SURF_MODE=rayleigh
# 性能与缓存（可选）
export QD_ECO_LIGHT_UPDATE_EVERY_HOURS=6
```

说明
- I_b(t) 由双星入射与瑞利/云调制降维得到；生态返回带反照率 A_b^surface 后按带能量权重降维为标量 α_surface，再与云/海冰在辐射路径中合成最终 albedo。
- 若仅做日级诊断，设 `QD_ECO_SUBDAILY_ENABLE=0` 且 `QD_ECO_FEEDBACK_MODE=daily`，仅在日末写回。
