# Project 015 (v2): 涌现生态学与光谱演化动力学（Emergent Ecology & Spectral Evolutionary Dynamics）

状态（2025‑09‑25）
- [x] v1 文档定稿
- [x] v2 方案修订与规范化（本文件）
- [ ] M1：气候—生态耦合接口（GCM Weather → EcologicalEnvironment）
- [ ] M2：虚拟生命内核（Genes / Plant FSM）
- [ ] M3：资源竞争模块（光/水）与 PopulationManager
- [ ] M4：光谱物理引擎（恒星光谱/植物吸收/反射聚合）
- [ ] M5：端到端生态演化运行与面板/诊断

关联文档
- docs/13-plant-model.md（虚拟生命内核：基因/个体/FSM）
- docs/14-adaptive-spectroscopy.md（适应性光谱物理：恒星光谱、吸收、反射）
- docs/15-ecology-framework.md（生态系统组织/竞争/传播/聚合反馈）
- docs/06-energy-framework.md（能量收支；地表反照率/短波/长波耦合）
- docs/04-runtime-config.md（运行参数总表；生态模块参数将在后续补入目录）
- scripts/run_simulation.py（主循环集成位点）

目标与范围
- 构建一个“自下而上”的涌现生态系统：从个体（Plant）—种群（Population）—群落与全球格局。
- 让植物光谱特征与双星动态光谱环境相互适应：颜色/吸收峰为演化结果而非先验常量。
- 双向耦合气候：生态聚合反射光谱 → 地表光谱反照率（短波路径）；水分/蒸腾影响湿度循环（演进阶段）。
- 提供可复现实验平台：形态—策略权衡、光谱—环境匹配、群落稳定性等问题的计算探究。

非目标（当前里程碑外）
- 动物/昆虫授粉与复杂食物网
- 土壤/微生物/地下水三维过程
- 植被—边界层湍流高级耦合（未来可与 P006/P008 扩展）

1. 架构与数据流（每日时间步）
1) GCM → Weather
- 来自主循环的每日/子日诊断：地表温度 Ts、近地气温 Ta、光照（光谱/或带平均）、云与降水、土壤水分代理（起步可用 P009 陆地桶 W_land 归一化）、风速与湿度指标等。
- 双星几何与大气调制产生地表光谱辐照度 Irradiance_Surface(λ, t) 或其离散波段 I_b(t)。

2) Weather → EcologicalEnvironment（每格点）
- PopulationManager 计算每个植物的“生态微环境”：
  - 可用光比例 light_availability_i（冠层遮挡/分层衰减）
  - 水分份额 water_share_i（根系权重/桶量约束）
  - 当日天气与光谱（或带平均 I_b），用于光合能量与胁迫更新

3) Plant.update_one_day(env_i)
- 依据 FSM（SEED/GROWING/MATURE/SENESCENT/DEAD）与内部记忆（积温/胁迫天数）更新
- 光谱积分（或带和）获得日能量收入 → 投资根/茎/叶/繁殖/储能
- 计算个体反射光谱 Reflect_i(λ) 或带反射 R_b,i 与叶面积/覆盖度

4) PopulationManager → Surface Spectral Albedo
- 聚合所有个体与土壤背景得到 Albedo_Surface(λ) 或带反照率 A_b
- 降维反馈给 GCM：短波模块使用 A_b → SW_atm/SW_sfc；TrueColor 使用全谱→RGB 可视化

5) 年度/季末 → 繁殖与传播
- 收集成熟个体的种子；大部留在本格“种子库”，少量按扩散模型远距传播至相邻格
- 发生微小突变 → 新基因型进入下一季竞争

2. 接口与模块（概览）
2.1 Weather / EcologicalEnvironment（GCM → Ecology）
- Weather（每日）：
  - Ts, Ta, wind10, RH/q 代理, precip_daily, soil_water_index（0..1）, spectral_irradiance 或 I_b[NB]
  - cloud_optical proxy（可选，用于感知光谱调制）
- EcologicalEnvironment（个体视角）：
  - effective_light_bands_i[NB] = I_b × light_availability_i（来自冠层模型）
  - water_supply_i（由水竞争权重与 soil_water_index/可用水决定）
  - day_length/seasonal_phase（可选，用于物候/状态转换）

2.2 Genes / Plant / PopulationManager（见 docs/13 / docs/15）
- Genes：身份/形态/光谱/生理/生命周期参数；光谱基因为高斯峰列表或样条系数
- Plant（FSM）：update_one_day(environment)；返回每日能量收支与反射特征
- PopulationManager：
  - step_daily(weather) → 更新内部所有 Plant
  - 竞争：分层光照衰减（类 Beer-Lambert）；根系权重水分分配
  - 聚合：反射光谱/带反照率；计算 LAI/总生物量/年龄结构等诊断
  - 繁殖/传播：本地播种 + 远距扩散 + 突变

2.3 Spectral Engine（见 docs/14）
- 输入：双星地表光谱 I(λ) 或带强度 I_b
- 植物吸收：Absorb_i(λ) 由基因生成（高斯峰叠加，截断至 ≤1）
- 能量收入（日）：E_gain_i = ∫ I(λ)·Absorb_i(λ)dλ ≈ Σ I_b·A_b,i·Δλ_b
- 反射：Reflect_i(λ) = 1 − Absorb_i(λ)
- 聚合：Albedo_Surface(λ) = 加权[Σ Reflect_i(λ)·leaf_area_i + soil(λ)·area_soil] / area_total
- 降维反馈：A_b = ∫ Albedo_Surface·W_b dλ（短波带权重）；RGB 可视化使用 CIE 变换

3. 环境变量（建议；将在 docs/04 增补目录）
主开关与步长
- QD_ECO_ENABLE（默认 0）：1 开启生态模块（每日步）
- QD_ECO_DT_DAYS（默认 1.0）：生态每日更新步长（天）

光谱与带宽（与 docs/14 对齐）
- QD_ECO_SPECTRAL_BANDS（默认 8）：短波离散波段数
- QD_ECO_SPECTRAL_RANGE_NM（默认 380,780）：可见范围
- QD_ECO_TOA_TO_SURF_MODE（默认 simple）：地表光谱调制模式（simple|rayleigh|custom）
- QD_ECO_SOIL_SPECTRUM（可选路径）：土壤背景光谱文件

竞争与形态
- QD_ECO_LIGHT_K（默认 0.5）：冠层光衰减系数
- QD_ECO_WATER_PRIORITY（默认 root_mass）：水分竞争权重（root_mass|depth_weighted）
- QD_ECO_STRESS_WATER_DAYS（默认 7）：累计水分胁迫进入衰老阈值（天）

演化与传播
- QD_ECO_SEED_BANK_MAX（默认 1000）：本地种子库容量
- QD_ECO_LONGDIST_FRAC（默认 0.05）：远距离传播比例
- QD_ECO_MUT_RATE（默认 1e-3）：突变概率（每种子）
- QD_ECO_INIT_SPECIES（默认 grass,tree）：初始化基因型集合（逗号分隔）

耦合与反馈
- QD_ECO_ALBEDO_COUPLE（默认 1）：将聚合反照率反馈到短波路径
- QD_ECO_LAI_ALBEDO_WEIGHT（默认 1.0）：叶面积在反照率聚合中的权重
- QD_ECO_DIAG（默认 1）：生态诊断打印与图层输出

4. 每日时间步顺序（与主循环的契合）
建议在 scripts/run_simulation.py 中（高层）按下列时序集成（与 docs/12 一致）：
1) 轨道/光谱入射：I(λ,t) 或 I_b（docs/02 + docs/14）
2) 动力步 + 湿度步 + 能量步（docs/06/07/08/10）
3) Ecology 步（每日触发；若当日跨越多气象步，聚合当日气象量后再调用）
   - 从 GCM 收集 Weather（当日平均或日中代表）
   - PopulationManager.step_daily(weather)
   - 产出 Albedo_Surface → 计算带反照率 A_b → 回写至能量模块（短波）
4) 海洋步 / 水文闭合 / 路由（docs/07/09/14）
5) 可视化与诊断：增加“生态”层（LAI/总生物量/颜色分布等）

4A. 时级接口与即时回耦（Sub‑daily / Hourly）
为满足“与物理步同步的快速耦合”，生态模块除日级接口外，新增时级接口与双时序方案（见 docs/15 v2.1 与 docs/13）：
- WeatherInstant（时级）：I_b(t)[NB]、Ts、RH、precip_rate、soil_moisture_proxy 等
- PopulationManager.step_subdaily(weather_inst, dt_seconds)：
  - 依据策略（子采样/缓存）选择是否重算冠层衰减/聚合反照率
  - 调用 Plant.update_substep 累积当日能量/胁迫缓冲（不做形态大跳转）
  - 返回 A_b^surface（可选），用于“即时短波回耦”（FEEDBACK_MODE=instant）
- 缓存与子采样：
  - QD_ECO_SUBSTEP_EVERY_NPHYS 控制每 N 个物理步调用一次
  - QD_ECO_LIGHT_UPDATE_EVERY_HOURS/RECOMPUTE_LAI_DELTA 控制冠层与反照率重算频率
- 回耦策略：
  - QD_ECO_FEEDBACK_MODE=instant：时级返回即写回短波带反照率，下一物理步生效
  - QD_ECO_FEEDBACK_MODE=daily：仅在日末写回（慢路径/诊断）
- 相关新增环境变量请见 docs/15 v2.1 §6 与后续将补充的 docs/04“生态模块（P015）”条目

5. 任务拆解与里程碑
M1：接口层（Weather ↔ Ecology）
- [ ] 在 scripts/run_simulation.py 新增 EcologyAdapter：负责从 GCM 汇总每日 Weather 并传递给生态模块；负责将生态聚合反照率回写短波带
- [ ] 提供带宽定义与 I→带强度 I_b 的降维器（与 docs/14 对齐）

M2：虚拟生命内核（docs/13）
- [ ] 实现 Genes/Plant（FSM + 内部记忆 + 投资策略 + 能量学）
- [ ] 提供 absorb_curve 生成器（高斯峰叠加）

M3：生态框架（docs/15）
- [ ] 实现 PopulationManager（光/水竞争、聚合反照率、繁殖/传播/突变）
- [ ] 定义本地种子库与邻域传播 API；与土壤水分/桶量代理的接口

M4：光谱物理（docs/14）
- [ ] 地表光谱与植物吸收/反射的波段降维；CIE RGB 可视化
- [ ] 与能量模块短波带 α_total(A_b, cloud) 的耦合路径

M5：端到端运行与面板
- [ ] 新增 EcologyDiag：LAI、总生物量、平均反照率（带/可见）、颜色分布、种群年龄结构
- [ ] 输出光谱生物地理图（颜色/LAI/物种主导分布）
- [ ] 长期运行与稳定性/耦合回归

6. 诊断与验收标准
- 功能性
  - [ ] Ecology 开关可控；每日步与 GCM 步协调一致
  - [ ] 反照率反馈生效：开启后短波吸收/反射变化方向合理
  - [ ] 光/水竞争对个体能量与存活有可解释的影响

- 涌现行为
  - [ ] 形态分化：干旱区根系投资高的基因型占优；湿润区叶/茎投资更高
  - [ ] 光谱适应：优势基因的吸收峰与当地入射光谱峰（G/K 星混合）具相关性
  - [ ] 群落演替：草→灌/乔木的时空序列在部分区域观察到

- 数值与耦合稳定
  - [ ] Ecology 打开后能量闭合（docs/06）长期无系统漂移（TOA/SFC/ATM |净| < 阈值）
  - [ ] 与 P008/P009 的湿度/水量闭合指标保持在阈值内（长期平均）
  - [ ] 性能与内存无异常增长（可分辨率/步长扫描）

7. 运行示例（占位；接口落地后有效）
启用生态模块（每日步），带反照率反馈
```bash
export QD_ECO_ENABLE=1
export QD_ECO_DT_DAYS=1
export QD_ECO_ALBEDO_COUPLE=1

# 光谱/带设置（可选）
export QD_ECO_SPECTRAL_BANDS=8
export QD_ECO_SPECTRAL_RANGE_NM=380,780
export QD_ECO_LIGHT_K=0.5
export QD_ECO_LONGDIST_FRAC=0.05
export QD_ECO_MUT_RATE=1e-3

python3 -m scripts.run_simulation
```

仅诊断（不反馈短波）
```bash
export QD_ECO_ENABLE=1
export QD_ECO_ALBEDO_COUPLE=0
python3 -m scripts.run_simulation
```

8. 与其它模块的交叉引用
- docs/13（Plant）：本项目调用其 FSM/能量学/光谱吸收生成器
- docs/14（光谱）：带划分、吸收/反射、CIE 颜色、I→A_b 降维与短波耦合
- docs/15（生态）：PopulationManager 的光/水竞争、播散/突变、聚合反照率接口
- docs/06（能量）：α_total(A_b, cloud) 进入短波，长波可能间接受云/温室影响（后续）
- docs/04（运行参数）：生态参数将在后续补入目录，保持与本文件一致的命名

9. 变更记录（Changelog）
- 2025‑09‑25：v2 文档规范化，与 docs/13/14/15 的 API/耦合关系对齐；新增环境变量建议、每日步序列、验收标准与运行示例
