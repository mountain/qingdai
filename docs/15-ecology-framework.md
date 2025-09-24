# 文档 15（v2.1）：生态系统框架（Ecology Framework）

状态（2025‑09‑25）
- [x] v1：层次架构（个体—种群—全球）、光/水竞争、繁殖与传播、聚合反照率概念
- [x] v2：PopulationManager API、每日步序列、光/水竞争算法与带反照率聚合、环境变量与诊断、与 13/14/06/015 的交叉引用
- [x] v2.1：新增“时级接口（Sub‑daily/Hourly）”与双时序耦合方案（即时反照率回耦 + 日级慢路径）
- [ ] M3：参考实现与单元/端到端测试
- [ ] M4：参数扫描与默认组固化

交叉引用
- 项目 015：生态演化与气候耦合总览（projects/015-adaptative-spectral-biogeography.md）
- 文档 13：虚拟生命内核（Genes/Plant/FSM；新增 update_substep，docs/13-plant-model.md）
- 文档 14：适应性光谱（I_b(t) 即时带强度；日累计器，docs/14-adaptive-spectroscopy.md）
- 文档 06：能量收支（短波每物理步使用地表反照率；建议即时耦合，docs/06-energy-framework.md）
- 文档 04：运行参数目录（新增 QD_ECO_* 将在 04 汇总，docs/04-runtime-config.md）

0. 目标与范围
- 定义生态“中层控制器”PopulationManager，将“个体（Plant）”与“GCM 世界”解耦：
  - 从 GCM/Adapter 接收天气与光谱带 I_b
  - 在格点内执行光/水竞争与个体更新
  - 聚合个体反射带与土壤背景为地表带反照率 A_b^surface，反馈至短波
  - 管理繁殖/传播/变异与种子库，维持演化与扩张
- 新增：双时序接口
  - 时级（Sub‑daily/Hourly）：与物理积分步同步或按子采样频率，进行“快速能量/反照率耦合”（即时向短波反馈）
  - 日级（Daily）：进行“慢路径生物过程”（形态投资、状态转换、繁殖/传播、种子库管理等）

1. 数据结构与 API

1.1 PopulationManager（每格点一个实例；也可全局集中管理）
- 初始化
  - PopulationManager(grid_cell_area, spectral_bands, soil_albedo_bands, params)
  - 内部状态：
    - plants: List[Plant]
    - seed_bank
    - diag/rolling stats
    - caches（见 1.4；用于加速时级调用）
- 时级步（新）
  - step_subdaily(weather_inst: WeatherInstant, dt_seconds: float) -> Optional[np.ndarray[NB]]
    - 输入瞬时天气（见 1.3）；根据性能策略选择是否重算光/水竞争
    - 调用 Plant.update_substep（见文档 13）积累“当日能量/胁迫缓冲”
    - 计算/更新 A_b^surface（可采用缓存/低频重算；见 5.2）
    - 返回 A_b^surface（NB 带）用于“即时短波耦合”，或 None（若本步不重算）
- 日级步（原）
  - step_daily(weather_day: WeatherDaily) -> EcologyDailyReport
    - 按 2 的流程执行慢路径：光/水竞争（可重用缓存）、调用 Plant.update_one_day、清理/繁殖/传播/突变、聚合 A_b^surface、输出诊断
- 聚合与反馈
  - aggregate_surface_albedo(reports or state) -> A_b^surface[NB]
  - summary_metrics() -> dict（LAI、总生物量、年龄结构、主导物种等）

1.2 WeatherDaily（由 Adapter 聚合每日代表值；与 13/14 对齐）
- Ts, Ta, wind10, RH proxy, precip_daily
- soil_water_index ∈ [0,1]
- spectral_bands: I_b_day[NB]（日均或代表值；见文档 14）
- day_length_hours, seasonal_phase（可选）

1.3 WeatherInstant（时级输入；与 14/06 对齐）
- Ts, Ta, wind10, RH proxy, precip_rate
- soil_moisture_proxy 或 soil_water_index_inst（可选：由 P009 桶 + 降水速率近似）
- spectral_bands: I_b(t)[NB]（即时带强度；见文档 14）
- solar_zenith, cloud_proxy（可选，用于快速近似光场）

1.4 缓存与加速（建议）
- canopy_cache：最近一次计算的“分层光衰减权重”与排序，附带时间戳
- albedo_cache：最近一次 A_b^surface 及时间戳
- recompute policy：按 QD_ECO_LIGHT_UPDATE_EVERY_HOURS 或“LAI/结构阈值变化”触发重算

1.5 EcologyDailyReport（回主循环/可视化）
- A_b_surface: np.ndarray[NB]（日末聚合）
- LAI, biomass_total, n_plants, dominant_identity
- energy_gain_mean/std、seed_production、mortality
- 可选：Reflected RGB（TrueColor）

2. 双时序时间步顺序

2.1 物理积分循环（每步 dt_seconds）
1) 轨道/辐射/动力/湿度/能量步（docs/06/07/08/10）
2) 若 QD_ECO_SUBDAILY_ENABLE=1 且到达子采样频率：
   - Adapter 生成 WeatherInstant（I_b(t)、Ts、RH、soil_proxy…）
   - PopulationManager.step_subdaily(weather_inst, dt_seconds)
   - 若返回 A_b^surface：写回短波带反照率用于“下一物理步”
3) 海洋/水文/路由等

2.2 日界（或满 N 小时）触发慢路径
1) Adapter 聚合当日（或周期）WeatherDaily
2) PopulationManager.step_daily(weather_day)
3) 输出 EcologyDailyReport、更新日末 A_b^surface（用于诊断与可视化）
4) 清空当日缓冲（个体累计能量/胁迫），进入新日

3. 资源竞争算法（与 v2 相同，增加“时级策略”）

3.1 光竞争：分层光照衰减（类 Beer‑Lambert）
- 与 v2 相同的层级衰减公式（按 height 排序；I_top,b ← I_top,b · (1 − min(0.95, k_canopy · L_k / A_cell))）
- 时级策略（性能）：
  - 默认：仅在“canopy_cache 过期（超过 QD_ECO_LIGHT_UPDATE_EVERY_HOURS）或 LAI/结构差异超阈值”时重算
  - 其它物理步共用最近一次衰减权重；light_availability 标量近似可用于时级调用

3.2 水竞争：根系权重比例分配
- 与 v2 相同（按 root_mass × f_depth 权重），soil_water_index_inst 替换日均值
- 可选“可用水上限”V_avail_inst（由降水速率 + 桶量代理确定）→ 对 shares 归一后乘上限

3.3 物候与应激桥接（由 Plant 处理）
- 时级：update_substep 累积“当日能量缓冲 + 胁迫小时数”；不做生命周期跳转（除非“硬条件”）
- 日级：update_one_day 应用“投资/形态增长/状态转换/繁殖/传播”

4. 地表带反照率聚合与耦合（与 v2 相同，增加“即时回耦”）
- A_b^surface 公式不变：
  A_b^surface = [ Σ_i (R_b,i · W_i) + Albedo_soil_b · Area_soil ] / Area_total
- 时级（即时回耦）
  - 若 QD_ECO_FEEDBACK_MODE=instant：在 step_subdaily 返回 A_b^surface 后，Adapter 立即写回短波带反照率，用于下一个物理步
  - 缓存与低频重算：可设 QD_ECO_ALBEDO_COUPLE_FREQ=subdaily/daily 控制写回频率
- 日级（慢路径）
  - 日末聚合的 A_b^surface 用于诊断/可视化；若 FEEDBACK_MODE=daily，则仅此时写回

5. 数值与实现要点（时级）

5.1 Plant.update_substep（见文档 13）
- 输入 WeatherInstant + competition inputs
- 计算即时 E_gain_inst ≈ Σ_b I_b(t)·A_b·Δλ_b，累积到 energy_buffer
- 水分胁迫：soil_water_index_inst 与阈值比较，累积 water_stress_hours
- 形态投资/状态转换/繁殖：推迟到日级应用（除非设置硬触发）
- 反射带 R_b 与叶面积变化在时级通常不变（除非叶片快速受损/落叶）；因此 A_b^surface 的时级变化主要来自短波 I_b(t) 与几何/云变化（见 5.2）

5.2 何时需要时级重算 A_b^surface
- 叶面积/结构显著变化（通常日级）
- 冠层排序/光衰减需要更新（超时或阈值变化）
- 缓解计算：默认每 QD_ECO_LIGHT_UPDATE_EVERY_HOURS 重算一次 canopy；其它步直接复用最近 A_b^surface（因为 R_b 与 W_i 几乎不变）

5.3 与能量闭合（docs/06）
- 即时耦合时，短波吸收/反射应随 I_b(t) 与云场等即时变化；生态仅调整地表反射（而非直接改动大气参数）
- 长期能量闭合仍以日/多年平均评估（TOA/SFC/ATM）

6. 环境变量（汇总至文档 04；新增标注“时级”）

主控与步长
- QD_ECO_ENABLE（默认 0）
- QD_ECO_DT_DAYS（默认 1.0）
- QD_ECO_SUBDAILY_ENABLE（默认 1）：开启时级接口
- QD_ECO_DT_HOURS（默认自动=物理步 dt_seconds 的小时；用于诊断）
- QD_ECO_SUBSTEP_EVERY_NPHYS（默认 1）：每 N 个物理步执行 1 次 step_subdaily
- QD_ECO_FEEDBACK_MODE（instant|daily，默认 instant）
- QD_ECO_ALBEDO_COUPLE_FREQ（subdaily|daily，默认 subdaily）

光竞争
- QD_ECO_LIGHT_K（默认 0.5）
- QD_ECO_LIGHT_UPDATE_EVERY_HOURS（默认 6）：时级重算冠层/反照率的最小间隔（小时）
- QD_ECO_LIGHT_RECOMPUTE_LAI_DELTA（默认 0.05）：若 LAI 相对变化超过阈值则强制重算

水竞争
- QD_ECO_WATER_PRIORITY（默认 root_mass）
- QD_ECO_SOIL_WATER_CAP（可选）

聚合与反馈
- QD_ECO_LAI_ALBEDO_WEIGHT（默认 1.0）
- QD_ECO_ALBEDO_COUPLE（默认 1）

演化与传播（慢路径）
- QD_ECO_SEED_BANK_MAX（默认 1000）
- QD_ECO_LONGDIST_FRAC（默认 0.05）
- QD_ECO_MUT_RATE（默认 1e‑3）

诊断
- QD_ECO_DIAG（默认 1）

7. 诊断与验收标准（更新，含时级）

功能性
- [ ] 时级与日级双路径均可单独/组合工作；即时回耦时短波响应与 I_b(t) 变化一致
- [ ] QD_ECO_SUBSTEP_EVERY_NPHYS、LIGHT_UPDATE_EVERY_HOURS 与 FEEDBACK_MODE 可控生效

涌现行为
- 与 v2 相同（干旱/湿润分化、光谱适应、群落演替）；新增：
- [ ] 晴天/多云日内：短波吸收变化可通过 A_b^surface + I_b(t) 体现，并影响 Ts 的日变化相位/幅度（趋势合理）

数值与耦合稳定
- [ ] 即时回耦开启后，多年平均 |TOA/SFC/ATM 净| 维持在阈值内（如 2 W·m⁻²）
- [ ] 性能可控：通过子采样与缓存策略，时级调用不会显著拖慢总步

8. 运行示例（占位）

即时回耦（推荐）
```bash
# 物理步：dt_seconds 按现有配置
export QD_ECO_ENABLE=1
export QD_ECO_SUBDAILY_ENABLE=1
export QD_ECO_SUBSTEP_EVERY_NPHYS=1
export QD_ECO_FEEDBACK_MODE=instant
export QD_ECO_ALBEDO_COUPLE_FREQ=subdaily
export QD_ECO_LIGHT_UPDATE_EVERY_HOURS=6
python3 -m scripts.run_simulation
```

仅日级（诊断/开发）
```bash
export QD_ECO_ENABLE=1
export QD_ECO_SUBDAILY_ENABLE=0
export QD_ECO_FEEDBACK_MODE=daily
python3 -m scripts.run_simulation
```

9. 与其它模块的边界与注意事项
- 能量（文档 06）：短波每物理步读取地表反照率；建议即时回耦下使用最新 A_b^surface（若无更新则沿用缓存）
- 湿度（文档 08）：蒸腾尚未显式加入；时级路径未来可将“植被覆盖/气孔行为”作为蒸发权重修正
- 水量（文档 09/14）：生态与路由的耦合可在后续通过“植被调制的蒸发/截留”体现；当前保持最小一致
- 数值稳定：日末清空个体缓冲；注意单位（E_gain 累积/小时计数）一致性

10. 变更记录（Changelog）
- 2025‑09‑25：v2 → v2.1：新增时级接口（PopulationManager.step_subdaily/WeatherInstant）、子采样与缓存策略、即时反照率回耦、双时序调度与环境变量；与 13/14/06/015 一致化
