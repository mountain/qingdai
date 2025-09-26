# 项目 P018：Plant v2 实施方案（虚拟植物模型落地与性能工程）

状态：进行中（M2/M2.5 已完成；M3（基因导出/日界集成）与 M4（扩展 autosave/周期保存）已完成核心项）  
版本：v2.3（2025‑09‑26）  
负责人：Cline（AI 软件工程师）  
关联：  
- 文档 13（个体/基因与 FSM）：[docs/13-plant-model.md](../docs/13-plant-model.md)  
- 文档 14（光谱带与 I(λ)→I_b）：[docs/14-adaptive-spectroscopy.md](../docs/14-adaptive-spectroscopy.md)  
- 文档 15（生态框架/PopulationManager/双时序）：[docs/15-ecology-framework.md](../docs/15-ecology-framework.md)  
- 运行参数目录（含 P015 条目）：[docs/04-runtime-config.md](../docs/04-runtime-config.md)  
- 现有代码（生态子包骨架）：`pygcm/ecology/`（adapter/spectral/plant/population/types 等）  
- 现有持久化产物（生态/水文/自动保存）：`data/eco_autosave.npz`、`data/species_autosave.json`、`data/phyto_autosave.npz`、`data/restart_autosave.nc`

### 实施状态更新（2025‑09‑26）
- 完成 M3（核心项）：基因导出 JSON schema v3（每个 gene 自包含 λ_* 带定义；顶层 bands 仅含 nbands/band_weights），在日界导出；新增冒烟脚本（`scripts/smoke_genes_export.py`）验证。
- 完成 M3+（个体池→种子库耦合）：IndividualPool 在日界将 per‑cell reproduction 能量写入 PopulationManager.seed_bank（受土壤门控与保留比例约束），随后在下一日 germination 转为 LAI；冒烟测试（`scripts/test_ecology_m3_plus.py`）通过。
- 完成 M4（核心项）：生态扩展 autosave（原子写 + 轮转备份）、周期保存（行星小时）、路径覆写与 Ctrl‑C 安全保存；新增轻量测试（`scripts/test_ecology_autosave.py`）验证（schema_version=1、roundtrip 通过、滚动备份生效）。
- 完成 M2：LAI 日级更新与带化回耦
  - PopulationManager 实现 E_day→growth 与 soil_water_index→senescence 的日级 LAI 更新（支持 K>1 分层捕获与分配）。
  - EcologyAdapter 在“日界”调用 get_surface_albedo_bands 并降维混入，减少逐步带化计算成本。
  - Ecology 面板新增 banded alpha 可视化。
- 完成 M2.5：空间扩张与生态可视化增强
  - 提供两种扩张机制：扩散扩张（vonNeumann/Moore 邻域）与种子传播（能量门控），参数可控，保守限幅。
  - 新增 canopy_height_map 与 species_density_maps 诊断；TrueColor 植被叠加使用带化反射与动态光谱权重。
- M1.5 辅助与持久化
  - 基因型吸收带缓存（Adapter 侧）降低谱线→带吸收计算开销。
  - 扩展 autosave（LAI、species_weights、bands、R_species_nb），向后兼容 legacy。
- 后续计划
  - 进入 M3（FSM/繁殖/传播联动个体池）与 M4（持久化 schema_version/RNG 等完善）阶段。

本文目标：将 Plant v2 的设计（文档 13/14/15）系统化落地到现有代码框架中，明确数据流、API 契约、性能策略、状态持久化方案与测试验收，确保“个体—种群—短波回耦”链路在小时级与日级均可稳定运行，且与能量/湿度/水文/海洋保持守恒一致与数值稳健；当作业中断时能无损续跑。


## 1. 背景与问题

当前仓库已具备：
- 生态模块框架（PopulationManager/Adapter/双时序接口）与光谱带协议（docs/15 v2.1、docs/14），运行开关与环境变量在 docs/04 已对齐。
- 代码骨架已存在（`pygcm/ecology/*.py`），但 Plant 个体的能量学仍存在“标量化能量权重”的遗留实现，不充分利用“基因→吸收带 A_b → 光谱带 I_b(t)”这一核心路径，导致：
  - 个体能量收入无法随“基因谱线 × 环境光谱”自洽变化；
  - 聚合反照率无法准确反映群落光学特征，削弱“生态→短波→气候”的反馈。
- 状态持久化尚未系统规划：中断后恢复的“最小必要状态”与“高性价比缓存”边界不清晰。

本提案修复与增强：
1) 将个体能量计算从“标量”升级为“带积分”（I_b × A_b × Δλ_b），用“基因型缓存 + 层级聚合”控制成本；  
2) 用 LAI/叶面积调制“有效吸收”，实现从个体形态到群落反照率的可解释路径；  
3) 将分层光竞争与水竞争收敛到 PopulationManager 层实现，并与“时级回耦/日级慢路径”统一在 Adapter 中调度；  
4) 新增“生态状态持久化与重启方案”，覆盖子步缓存/日级慢状态/RNG/物种库/缓存等关键持久量，定义自动保存节奏与原子写入策略。


## 2. 设计原则与数据流（与文档 13/14/15 一致）

核心原则（性能优先，物理一致）：
- 缓存（Cache）：与基因有关而与时刻/个体无关的计算（如吸收带 A_b）预计算一次，复用至多处，避免“个体×时间步”的高频重复。
- 聚合（Aggregate）：涉及群体交互的昂贵运算（冠层光衰减/分层竞争、地表带反照率聚合）在格点/层级统一计算，个体仅接收结果。
- 双时序（Sub‑daily + Daily）：  
  - 时级子步：只做“快速能量/光学回耦”（I_b(t) → A_b^surface → 短波）；不做个体状态大跳转；按缓存频率低频重算。  
  - 日级步：做“慢路径”（形态投资、生命周期转移、繁殖/传播/突变），产生日级诊断。

数据流（小时级回耦一圈）：
1) Forcing/光谱（docs/14）：计算地表短波带强度 I_b(t)（simple|rayleigh）。  
2) PopulationManager.step_subdaily(weather_inst, dt)：  
   - 选择性重算冠层衰减与带反照率聚合（由缓存与频率控制）；  
   - 返回 A_b^surface（可选）；  
3) Adapter 将 A_b^surface 立即降维为标量 α_surface_ecology（按 I_b 权重）并融合进能量模块的地表短波路径（docs/06）；  
4) 下一物理步使用更新后的 α_total（云/海冰/生态融合）。

数据流（日级慢路径）：
1) Adapter 聚合一日代表天气与 I_b_day；  
2) PopulationManager.step_daily(weather_day)：  
   - 先执行光/水竞争分配，再对每个 Plant 调用 update_one_day（吸收带 × I_eff 带积分获取能量 → 形态投资/FSM → 反射带/LAI 输出）；  
   - 聚合得到 A_b^surface_day 与生态诊断；  
3) 若 FEEDBACK_MODE=daily，则在此时写回带反照率。


## 3. 模块与 API 契约（与现有文件映射）

代码目录（现有 + 本提案不改变总体布局）：
```
pygcm/ecology/
  spectral.py     # 带定义/积分与 A_b 生成（文档 14）
  plant.py        # Plant（FSM + 每日步 + 子步缓冲）
  population.py   # PopulationManager（光/水竞争 + 聚合 + 缓存）
  adapter.py      # EcologyAdapter（双时序编排与短波回写）
  types.py        # 数据结构（WeatherInstant/WeatherDaily/Reports）
  genes.py        # Genes/Genotype 定义与加载
  diversity.py    # 物种集合/初始化/传播策略（与 docs/15 附录一致）
```

关键接口（新增/规范化）：

- spectral.py
  - make_bands(nbands, λ0, λ1, mode="uniform|adaptive") → SpectralBands
  - absorbance_from_genes(spectral_peaks, bands) → A_b[NB]  // 预计算用
  - toa_to_surface_bands(I_total, cloud, bands, mode) → I_b[NB]

- population.py
  - class PopulationManager(...):
    - genotype_absorbance_cache: Dict[genotype_id, A_b[NB]]
    - step_subdaily(weather_inst, dt_seconds) → Optional[A_b^surface[NB]]
    - step_daily(weather_day) → EcologyDailyReport
    - aggregate_surface_albedo(reports|state) → A_b^surface[NB]
    - summary_metrics() → dict

- plant.py
  - class Plant:
    - update_substep(env_inst, dt_seconds, *, A_b_genotype, bands) → SubstepReport
    - update_one_day(env_day, *, A_b_genotype, bands) → PlantDailyReport
    - effective_leaf_area()、get_reflectance_bands(bands)（反射带 R_b=1−A_b_effective）

- adapter.py
  - class EcologyAdapter:
    - step_subdaily_if_due(world_state, dt) → Optional[A_b^surface]
    - step_daily_if_boundary(world_state) → EcologyDailyReport
    - write_back_albedo_bands(A_b^surface) → None
    - 降维：α_surface_ecology = Σ_b A_b^surface · (I_b / Σ I_b)

约束与契约：
- Plant 不直接计算“谱线→A_b”，而是由 PopulationManager 从 cache 提供 `A_b_genotype`；  
- Plant 的“时级子步”只累积“当日能量/胁迫缓冲”，不改变形态/FSM；  
- PopulationManager 负责光/水竞争、聚合反照率与缓存策略；  
- Adapter 是唯一的“短波回写”入口，负责写回节奏（instant/daily）与降维。


## 4. 性能工程与数值稳定

- 基因型吸收带缓存（核心节流）：  
  - 生成时机：初始化 PopulationManager 时，对“初始种群/物种库”一次性计算；  
  - 扩展：新物种（突变/引入）出现时懒加载计算并缓存；  
  - 数量级：通常 Ns≈10–50，NB≈8–16；计算成本极低。

- 冠层光竞争缓存：  
  - canopy_cache 保存“上次层级排序+衰减权重+LAI 层统计”，按  
    - QD_ECO_LIGHT_UPDATE_EVERY_HOURS（默认 6 小时）  
    - 或 LAI 相对变化阈值（QD_ECO_LIGHT_RECOMPUTE_LAI_DELTA，默认 0.05）  
    触发重算；其余时刻直接复用，避免“每步排序”。

- 个体能量与反射：  
  - 每个体仅进行 NB 次乘加（I_eff × A_b × Δλ_b）与 1 次指数（LAI 调制，见下一节）→ 开销线性、可控。

- 数值稳定：  
  - 带积分采用带平均/权函数，避免逐 nm 采样的刚性；  
  - A_b/R_b 与 α_total 严格限幅 [0,1]；  
  - 时级回耦仅改变地表短波反照率分量，不直接改动大气参数；  
  - 与 P006/P008/P009 的能量/水量闭合由统一诊断守护（docs/06/08/09）。


## 5. 关键算法细节（v2 规范落地）

5.1 个体能量学（带积分）
- 预计算吸收带：A_b = absorbance_from_genes(genes.spectral_absorption_curve, bands)  
- 日能量收入（示意，文档 13 一致）：
  ```
  # light_availability, water_share 等由 PopulationManager 竞争分配
  I_eff = I_b * light_availability  # [NB]
  # LAI/叶面积调制（自遮蔽）：f_LAI = 1 − exp(−k_ext * LAI)
  A_b_effective = A_b * f_LAI(LAI)
  E_gain = Σ_b I_eff[b] * A_b_effective[b] * Δλ_b
  ```
- 反射带：R_b = 1 − A_b_effective（限幅 [0,1]）

5.2 LAI/叶面积调制（非线性、可解释）
- 推荐默认：`f_LAI = 1 − exp(−k_ext * LAI)`，k_ext≈0.3–0.6（可配置/物种差异）  
- 作用：小 LAI 时吸收近线性增长；LAI 大时渐饱和（自遮蔽）

5.3 冠层光竞争（PopulationManager）
- 对格点内个体按 height 降序分层；  
- 层间透过：`I_b,next = I_b,top * exp(−k_canopy * LAI_layer / A_cell)`（或离散化等效形式）；  
- 层内按叶面积或 LAI 权重分配 light_availability。

5.4 地表带反照率聚合
- 与文档 14/15 一致：  
  `A_b^surface = [Σ_i (R_b,i * W_i) + Albedo_soil_b * Area_soil] / Area_total`  
  - W_i 可用“leaf_area × QD_ECO_LAI_ALBEDO_WEIGHT”  
  - 土壤光谱可来自 `QD_ECO_SOIL_SPECTRUM` 或默认曲线

5.5 短波回耦（Adapter）
- 时级：若 FEEDBACK_MODE=instant 且本步重算，立即回写 A_b^surface；  
- 日级：若 FEEDBACK_MODE=daily，仅日末回写；  
- 降维策略：按当前 I_b 比例或固定 w_b^SW 将带反照率转为标量 α_surface_ecology 后与海冰/云/基底合成 α_total（docs/06）。

5.6 状态持久化与重启（新增）
为确保中断后可“无损续跑”，定义生态模块的持久化分层：最小必要状态（必须保存）+ 性价比缓存（可选保存）+ 诊断/元数据（推荐保存）。

A) 必须保存（Minimal Required State）
- 时间与随机性
  - t_seconds/t_step_index：当前模拟时刻（与 GCM 主时钟对齐）
  - rng_state：生态侧 RNG 状态（与 QD_ECO_RAND_SEED 一致化）
- 个体/种群状态（每格/每个体或聚合表）
  - Plant：state(FSM)、age_in_days、energy_storage、biomass{root,stem,leaf}、leaf_area、height、genes_id
  - 当日缓冲：last_day_energy_gain、water_stress_hours（用于日末慢路径）
- 群落/空间结构
  - 每格种群列表（轻量结构，引用 genes_id；数值字段用于重建）
  - 种子库 seed_bank（上限、物种类别、年龄/能量）
- 生态—短波接口
  - A_b^surface（上一回写值，若 FEEDBACK_MODE=instant，避免恢复步“突变”）
  - α_surface_ecology（标量降维值，便于落地到能量模块）
- Adapter 辅助
  - subdaily_accum_hours（已累计时长，用于 hodin 窗口与子采样频率）
  - last_canopy_recompute_ts（上次冠层重算时间戳）

B) 可选保存（High‑ROI Caches）
- genotype_absorbance_cache（基因型→A_b），可复算，但保存可提速
- canopy_cache（分层顺序与衰减权重），在 LAI 变化不大时重用
- soil_spectrum_bands（若读取外部土壤光谱）

C) 诊断与元数据（Recommended）
- LAI/biomass/主导物种/生态色彩统计的滚动窗口
- bands 定义（λ_edges/Δλ/权重）、FEEDBACK_MODE、子步频率
- 参数快照（关键 QD_ECO_*、QD_* 交叉模块参数），便于重现实验

D) 文件与格式
- 自动保存（高频、小体积）：`data/eco_autosave.npz`（Numpy NPZ，结构化数组/压缩）
  - 原子写入策略：先写临时文件 `*.tmp`，`fsync` 后 `os.replace` 原子替换
  - 滚动备份：`eco_autosave_YYYYMMDDTHHMMSSZ.npz` 最多保留 `QD_ECO_AUTOSAVE_KEEP` 份
- 日级快照（低频、全量）：`data/eco_daily_YYYYMMDDTHHMMSSZ.npz`
- 兼容主重启（选项）：在 `data/restart_equilibrium.nc` 中写入生态小节（NetCDF group），字段名与 NPZ 保持一致（可后续里程碑）

E) 触发与时机
- 子步：每累计 `QD_ECO_AUTOSAVE_EVERY_HOURS`（默认 6 h）触发一次 autosave
- 日界：完成 `step_daily` 后写日级快照（可选）
- 结束/中断：
  - 正常结束：`on_exit` 钩子写一次 autosave（若 `QD_ECO_AUTOSAVE_ON_EXIT=1`）
  - SIGINT/SIGTERM：在主循环外层 `try..except` 捕获，调用 Adapter.save_autosave()，保证生态侧与能量侧同步落盘

F) 恢复流程（Resume Semantics）
1. 若 `QD_AUTOSAVE_LOAD=1` 且 `data/eco_autosave.npz` 存在：  
   - 读取生态状态；校验时间戳与主重启文件（若存在）的 t_seconds 一致；不一致则以主重启为源，生态按照“最接近一致”的策略平滑恢复（A_b^surface/α_surface_ecology 用 autosave 值启动，下一子步内渐变）
2. 若 autosave 不存在：用物种库与默认初值初始化（docs/15）；记录提示
3. 缺省字段/版本迁移：通过 `schema_version` 映射做向后兼容；缺失字段用安全默认值填充

G) 关键环境变量（新增与复用）
- QD_AUTOSAVE_LOAD（默认 1）：是否自动加载 autosave
- QD_ECO_AUTOSAVE_EVERY_HOURS（默认 6）：子步自动保存间隔（小时）
- QD_ECO_AUTOSAVE_PATH（默认 data/eco_autosave.npz）
- QD_ECO_AUTOSAVE_KEEP（默认 4）：滚动备份数量
- QD_ECO_AUTOSAVE_ON_EXIT（默认 1）：退出时保存
- QD_ECO_RAND_SEED（已有）：确保恢复后 RNG 行为可复现

H) 与其它模块一致性
- 与能量（docs/06）：恢复后 α_surface_ecology 与上次写回保持连续，避免短波突跳
- 与湿度/水量（docs/08/09）：生态状态恢复不改变 E/P/R 的历史账目，仅影响未来时刻；诊断面板重置窗口时需标记
- 与水文路由/湖泊（docs/14/P014）：若湖面按水体处理，生态恢复不破坏路由缓存（由 routing 自己持久化）；生态侧仅记录“水体类型掩膜”版本号

I) 单元测试（持久化）
- save→load 等价：内存状态与恢复状态 key-by-key 比较
- 原子写中断模拟：半写入文件不应破坏上一个可用 autosave
- 不同 schema_version 迁移：老文件能被读取并补齐默认字段


## 6. 环境变量与默认（与 docs/04 P015 条目一致，新增持久化项）

主控与步长（小时级回耦推荐）
- QD_ECO_ENABLE=1  
- QD_ECO_SUBDAILY_ENABLE=1  
- QD_ECO_SUBSTEP_EVERY_NPHYS=1  
- QD_ECO_FEEDBACK_MODE=instant|daily（默认 instant）  
- QD_ECO_ALBEDO_COUPLE=1  
- QD_ECO_SPECTRAL_BANDS（默认 16，建议 8–16）  
- QD_ECO_TOA_TO_SURF_MODE=simple|rayleigh（默认 rayleigh）  
- QD_ECO_LIGHT_UPDATE_EVERY_HOURS=6  
- QD_ECO_LIGHT_RECOMPUTE_LAI_DELTA=0.05  
- QD_ECO_LAI_ALBEDO_WEIGHT=1.0  

持久化（新增）
- QD_AUTOSAVE_LOAD=1  
- QD_ECO_AUTOSAVE_EVERY_HOURS=6  
- QD_ECO_AUTOSAVE_PATH="data/eco_autosave.npz"  
- QD_ECO_AUTOSAVE_KEEP=4  
- QD_ECO_AUTOSAVE_ON_EXIT=1  

LAI/个体参数（建议）
- QD_ECO_LAI_K_EXT=0.4  
- 其余见 docs/13/15/04


## 7. 里程碑与任务拆解

M1 基线接通（谱带能量 + 缓存 + 回耦）
- [x] spectral.absorbance_from_genes：实现/完善高斯峰→带吸收 A_b  
- [x] PopulationManager：genotype_absorbance_cache；子步回耦返回 A_b^surface（带）  
- [x] Plant.update_one_day：改为使用 A_b_genotype 与带积分；移除旧“标量能量权重”路径  
- [x] Plant.update_substep_bands：带积分接口（A_b_genotype × f_LAI）已实现  
- [x] Adapter：instant/daily 回写与降维策略完备；日志标注 α_surface_ecology 来源  
- [x] 单元测试：A_b 生成一致性、带积分能量与限幅  

M2 LAI 调制与分层光竞争
- [x] Plant：新增 f_LAI（参数化 K_EXT）；R_b=1−A_b_effective 限幅  
- [x] PopulationManager：按高度分层、层透过/层内分配  
- [x] PopulationManager：canopy_cache 与重算策略  
- [x] 单元测试：层透过/分配守恒、LAI→能量饱和曲线合理  

M3 日级慢路径与 FSM/繁殖
- [x] Plant.update_one_day：形态投资、FSM 转换、繁殖/种子产量接口  
- [x] PopulationManager：日级 step_daily 聚合、种子库/传播（与 docs/15 附录策略一致）  
- [x] Adapter：日界驱动、日级诊断（LAI/生物量/主导物种/反照率）  
- [x] M3+：IndividualPool 每日繁殖能量→ seed_bank 耦合并冒烟测试（scripts/test_ecology_m3_plus.py）  
- [ ] 单元测试：生命周期转移与基本能量账目一致  

M4 持久化与重启（本次新增重点）
- [x] 定义 autosave schema（minimal/caches/diagnostics/schema_version=1）  
- [x] Adapter.save/load API + 原子写入/滚动备份/退出钩子  
- [ ] 与主重启文件 t_seconds 对齐与不一致处理（平滑恢复）  
- [x] 单元测试：save→load 等价、半写入容错、版本迁移（`scripts/test_ecology_autosave.py` 冒烟通过）  

M5 集成与守恒验收
- [ ] 与能量（docs/06）：开启生态后 TOA/SFC/ATM 多年平均 |净| < 2 W·m⁻²  
- [ ] 与湿度（水汽/潜热，docs/08）：⟨LH⟩≈⟨LH_release⟩ 不恶化  
- [ ] 与水文（docs/09/14）：E–P–R 闭合与河网/湖泊图兼容  
- [ ] 集成测试脚本：短程烟测 + 多年运行回归 + 关键诊断阈值  

M6 可视化与调参套件
- [ ] TrueColor 植被叠加（docs/15 参数）；颜色/带反照率导出  
- [ ] 点位调试脚本：`scripts/plot_ecology_point.py` 增强（能量/LAI/反照率时序）  
- [ ] 参数扫描：K_EXT、LIGHT_K、LAI_ALBEDO_WEIGHT 的敏感性与默认组固化  

M7 文档与示例
- [ ] 更新 docs/13/14/15 的“实现说明/运行示例”片段（与代码一致）  
- [ ] README：加入 P015 小节“启用生态小时级回耦”的最小命令组与持久化说明  
- [ ] 变更记录（Changelog）补记  

M8 性能基准与预算
- [ ] NB=8/16、网格 181×360、个体规模（每格 10/50/100）下的额外耗时占比  
- [ ] 子步重算频率与 LAI 阈值对成本/精度折中曲线  
- [ ] 目标：生态开启后的总步耗时增加 ≤ 10–15%（NB=16、默认缓存策略）


## 8. 验收标准与测试计划

功能性
- [ ] 生态开关/回耦模式（instant/daily）可控；小规模与默认网格均可稳定运行  
- [ ] 聚合反照率随群落结构（LAI/谱线）变化方向正确；TrueColor 渲染与带反照率一致  

数值与守恒
- [ ] 能量：多年平均 |⟨TOA_net⟩|、|⟨SFC_net⟩|、|⟨ATM_net⟩| < 2 W·m⁻²  
- [ ] 潜热：⟨LH⟩（SFC）≈ ⟨LH_release⟩（ATM）  
- [ ] 水量：E–P–R 闭合误差与无生态基线同阶  

性能与稳定
- [ ] NB=16、默认缓存策略下性能开销在预算内；无内存泄露  
- [ ] 时级回耦无数值振荡；极端光谱与阴雨连日场景无爆裂  

持久化
- [ ] 中断后恢复 run-to-run 等价（阵列字段数值匹配在容差内）  
- [ ] 半写入/损坏 autosave 的回退与日志提示  
- [ ] 跨版本 schema 迁移测试通过  

测试套件
- 单元：spectral/plant/population 的吸收/反射/分配/投资/持久化逻辑  
- 集成：Adapter 回写 + 能量闭合回归 + 诊断阈值报警  
- 端到端：短程烟测（小时级）、中期（年级）与基准回放比对


## 9. 风险与对策

- 性能瓶颈：大量个体 × 高频子步 → 用“基因型缓存 + 冠层缓存 + 子采样频率”化解；提供 NB/频率/阈值旋钮  
- 数值刚性：带积分与指数项叠加 → 严格限幅、平滑阈值、避免逐 nm 刚性  
- 反馈放大：生态→短波→温度→生态的正反馈 → 以 FEEDBACK_MODE=daily 与权重限幅做渐进试验  
- 复现性：生态传播/突变的随机性 → QD_ECO_RAND_SEED、autosave 策略与文档化  
- 接口漂移：Plant/PopulationManager/Adapter 的契约变更 → 本文档 API 契约为准，修改需同步 docs/12/15  
- 持久化一致性：生态 autosave 与 GCM 主重启文件时间戳不一致 → 以主时钟为准并平滑恢复 α_surface_ecology；打印黄色告警与差值

## 10. 代码改动清单（最小侵入）

- pygcm/ecology/spectral.py  
  - [新增/完善] make_bands/absorbance_from_genes/toa_to_surface_bands  
- pygcm/ecology/plant.py  
  - [重构] update_one_day/ update_substep：改为带积分路径，接收 A_b_genotype 与 bands  
  - [新增] f_LAI 调制、R_b 计算与限幅、子步能量缓冲  
- pygcm/ecology/population.py  
  - [新增] genotype_absorbance_cache、canopy_cache；分层光竞争与聚合  
  - [新增] step_subdaily/step_daily；aggregate_surface_albedo  
- pygcm/ecology/adapter.py  
  - [完善] 双时序调度；回写频率（instant/daily）；降维与日志  
  - [新增] save_autosave/load_autosave（原子写、滚动备份、schema_version）  
- pygcm/ecology/types.py / genes.py / diversity.py  
  - [补充] 数据结构字段（LAI、height、seed/repro、谱峰列表）与初始化策略  
- scripts/plot_ecology_point.py  
  - [增强] 输出 E_gain/LAI/反照率时序与带谱诊断  
- docs（13/14/15/04）  
  - [同步] API/变量与运行示例与持久化说明  
- scripts/run_simulation.py  
  - [核对] EcologyAdapter 调用时序与环境变量读取；能量模块短波读入 α_surface_ecology 融合无重复  
  - [新增] 主循环 try..except 与 on_exit 保存调用；SIGINT/SIGTERM 时保存生态 autosave


## 11. 迁移与兼容

- 旧“标量能量权重”路径将移除/兼容层标记弃用；若存在依赖脚本，需迁移至带积分接口  
- autosave（data/eco_autosave.npz）加载优先于初始化；如需忽略请设 `QD_AUTOSAVE_LOAD=0` 或清理文件  
- 运行脚本与 README 的“启用生态”与“持久化/恢复”示例将以本提案为准更新


## 12. 快速运行示例

小时级即时回耦（推荐）
```bash
export QD_ECO_ENABLE=1
export QD_ECO_SUBDAILY_ENABLE=1
export QD_ECO_SUBSTEP_EVERY_NPHYS=1
export QD_ECO_FEEDBACK_MODE=instant
export QD_ECO_ALBEDO_COUPLE=1
export QD_ECO_SPECTRAL_BANDS=16
export QD_ECO_TOA_TO_SURF_MODE=rayleigh
export QD_ECO_LIGHT_UPDATE_EVERY_HOURS=6
export QD_ECO_AUTOSAVE_EVERY_HOURS=6
python3 -m scripts.run_simulation
```

仅日级（慢路径/诊断）
```bash
export QD_ECO_ENABLE=1
export QD_ECO_SUBDAILY_ENABLE=0
export QD_ECO_FEEDBACK_MODE=daily
export QD_AUTOSAVE_LOAD=1
python3 -m scripts.run_simulation
```

性能保守模式（较少重算 + 持久化）
```bash
export QD_ECO_LIGHT_UPDATE_EVERY_HOURS=12
export QD_ECO_SPECTRAL_BANDS=8
export QD_ECO_AUTOSAVE_EVERY_HOURS=12
python3 -m scripts.run_simulation
```


## 13. 变更记录（Changelog）

- 2025‑09‑26：v2.2 实装 M2/M2.5：完成 LAI 日级更新与带化回耦、扩散/种子传播、冠层高度与物种密度可视化；扩展 autosave（bands/R_species），向后兼容 legacy；更新运行脚本与面板。  
- 2025‑09‑26：v2.1 增补“状态持久化与重启”方案（5.6 与 6/7/10/12/验收项对应更新）；对齐 autosave 文件与退出钩子；完善恢复语义与跨版本兼容。  
- 2025‑09‑26：v2.0 提案定稿：带积分能量学、基因型吸收缓存、冠层分层竞争、双时序回耦、LAI 调制、验收与性能基准、代码改动清单与运行示例。


## 14. 附录：数据结构要点（摘要）

- WeatherInstant：I_b[NB]、Ts/Ta、wind、RH 代理、precip_rate、soil_water_index_inst  
- WeatherDaily：I_b_day[NB]、Ts/Ta 日均、precip_daily、soil_water_index、day_length  
- PlantDailyReport：energy_gain、reflectance_bands R_b[NB]、leaf_area、state 变更、seed_count  
- EcologyDailyReport：A_b^surface[NB]、LAI、生物量、主导物种、颜色/反照率诊断、传播统计  
- Autosave（schema_version=1）：{
  t_seconds, rng_state, plants[], populations[], seed_bank, A_b_surface, alpha_surface_ecology,
  subdaily_accum_hours, last_canopy_recompute_ts,
  caches?: {genotype_absorbance_cache, canopy_cache, soil_spectrum_bands},
  params_snapshot, bands_def
}
