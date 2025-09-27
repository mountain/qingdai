# 项目 P020（修订版 v1.1）：架构重构（面向对象 + 可测试 + 渐进迁移）

状态：提案→可落地设计（2025‑09‑27）  
作者：Cline（AI 软件工程师）  
关联文档与代码：
- docs/12-code-architecture-and-apis.md（现有架构与模块职责）
- scripts/run_simulation.py（现行主循环）
- pygcm/*（核心物理模块）
- projects/016-jax-acceleration.md（JAX 兼容与加速路径）
- 所有 docs/04（运行时参数目录与环境变量），docs/06/07/08/09/10/11/14/15/16/18（物理与生态设计）

本文件目标：将“面向对象重构”的愿景落为可执行蓝图，明确对象边界、API 合约、迁移阶段、测试与验收标准、兼容策略与性能预算，确保以最小风险渐进完成重构。


## 0. 设计原则（与现有仓库一致）

- 渐进迁移：禁止“大爆炸”式重写。新增 OO 层首先作为 façade 包裹现有实现，逐步内收。
- 单一职责：对象封装“状态 + 行为”；物理公式与数值算子尽量做成纯函数模块（stateless）。
- 可测试：每个模块都有可独立构造的最小状态与可重复调用的 API；提供契约/回归测试。
- 向后兼容：环境变量、CLI 用法、输出产物默认不变；新 API 通过适配层注入。
- 守恒优先：任何阶段的修改都不得破坏 docs/06/08/09 规定的能量/水量闭合标准。
- 性能预算：重构后端到端步时开销不增加（±5% 范围内）；新增层次不可引入不必要复制。


## 1. 概念边界（Configuration / Parameter / State）

与 docs/12 对齐，但在工程层面明确“生命周期与存储”：

- Configuration（配置，运行前确定，运行期不可变）
  - 例：网格分辨率（n_lat, n_lon）、dt、数值方案/开关（QD_FILTER_TYPE, QD_USE_OCEAN）。
  - 表达：不可变 dataclass（`SimConfig`）；序列化到元数据（restart/header）。

- Parameter（参数，定义规律，单次实验通常恒定，可由实验更改）
  - 例：物理常数、温室参数、Bowen 比、生态基因库、光谱带定义。
  - 表达：不可变/可版本化 dataclass（`PhysicsParams`, `EcologyParams`, `SpectralBands` 等）。

- State（状态，时间快照）
  - 例：u, v, h, Ts, Ta, q, cloud, SST, uo, vo, η, h_ice, W_land, SWE, routing buffers, LAI/生态状态等。
  - 表达：可变 dataclass（`WorldState` 拆分子状态）；支持全量/分组持久化；含 schema_version。


## 2. 目标对象模型（类与模块）

顶层对象（面板）：
- QingdaiWorld（有状态对象）
  - config: SimConfig
  - params: ParamsRegistry（聚合 PhysicsParams/EcologyParams/SpectralBands 等）
  - grid: Grid
  - subsystems: Atmosphere, Ocean, Surface, Hydrology, Routing, Ecology, Forcing（见下）
  - state: WorldState（聚合各子状态）
  - methods:
    - step(): 推进一个物理步，按统一时序协调子系统
    - run(n_steps|duration): 运行循环
    - save_state(path): 持久化（NetCDF/NPZ）
    - load_state(path): 恢复（包含 schema 迁移）
    - diagnostics(): 收集 EnergyDiag/HumidityDiag/WaterDiag/OceanDiag/EcologyDiag
    - plotters(): 调用 plotting 模块输出图像

子系统（有状态对象）：
- Atmosphere
  - state: u, v, h, Ta（或 h 代理 Ta）、cloud, q, 辅助缓冲等
  - apply_forcing(), time_step(), apply_diffusion(), diagnostics()
- Ocean
  - state: uo, vo, η, SST
  - apply_wind_stress(), step(), advect_tracers(), polar_fix(), diagnostics()
- Surface
  - state: Ts, h_ice（厚度/冰盖）、C_s_map、α_base 等衍生图
  - surface_energy_balance(), seaice_thermo(), albedo_fusion() （与 Physics 组合器协作）
- Hydrology（P009 + P019 融合的陆面/雪被/桶）
  - state: W_land, SWE_mm, C_snow, snow_age（可选）
  - partition_precip_phase_smooth(), snowpack_step(), land_bucket_step(), diagnostics()
- Routing（P014）
  - state: internal buffers, flow_accum_kgps, lake_volume ...
  - step(runoff_flux, precip_flux?, evap_flux?), diagnostics()
- Ecology
  - state: 生态个体/群落/LAI/seed banks/缓存；海色 phyto（可选分为 MarineEcology）
  - step_subdaily(), step_daily(), aggregate_surface_albedo_bands(), diagnostics()
- Forcing（无状态组合器 + 适配）
  - methods: orbital geometry, insolation(I), toa→surface spectral modulation
  - 仅依赖 time, params, grid

纯函数模块（stateless）：
- physics（docs/06 组合器）：shortwave, longwave, boundary_layer_fluxes, calculate_dynamic_albedo 等
- numerics（docs/10 算子）：laplacian_sphere, hyperdiffusion, shapiro, spectral_filters
- humidity（docs/08 核心算子）：q_sat, evaporation_flux, condensation
- hydrology_core（P009/P019 算子）：sigmoid phase partition, snow/melt, runoff
- ocean_core（P011 数值核）
- ecology_core（docs/13/14/15/16 带/吸收/聚合）  
- jax_compat（与 projects/016 一致）：xp, map_coordinates 替换件

绘图与 I/O（stateless）：
- plotting：统一状态 → 面板/TrueColor/OceanColor
- io_netcdf / io_npz：读写、schema、版本迁移工具


## 3. 目录与文件布局（不破坏现有，新增 façade 层）

新增（建议）：
```
pygcm/
  world/
    __init__.py
    world.py           # QingdaiWorld（主编排）
    config.py          # SimConfig 及解析环境变量的加载器
    params.py          # ParamsRegistry + 各 Params dataclass
    state.py           # WorldState + 子状态定义、schema_version
    forcing_facade.py  # Forcing façade（组合 orbital/physics/spectral）
    adapters.py        # 旧脚本/模块的适配层（渐进迁移）
```

保留并渐进重构：
```
pygcm/
  dynamics.py          # 逐步内收进 Atmosphere（可先做轻薄代理）
  ocean.py             # 逐步内收进 Ocean
  hydrology.py         # 逐步内收进 Hydrology
  routing.py           # 保持，Routing 类化
  energy.py, physics.py, humidity.py, topography.py, forcing.py, numerics（新增）
  ecology/             # 维持结构；对接 Ecology 子系统
  jax_compat.py        # 继续保留
scripts/run_simulation.py  # 外观不变，逐步调用 world.QingdaiWorld
```

说明：第一阶段，只“新增 world.FAÇADE + dataclass”；子系统先以“代理到现有模块函数/类”的方式实现，确保 0 风险接入。


## 4. API 合约（签名草案）

关键 dataclass（示意）
```python
# pygcm/world/config.py
from dataclasses import dataclass
@dataclass(frozen=True)
class SimConfig:
    n_lat: int
    n_lon: int
    dt_seconds: float
    plot_every_days: float
    # feature flags
    use_ocean: bool
    use_seaice: bool
    use_routing: bool
    use_ecology: bool
    # numerics
    filter_type: str  # 'combo'|'hyper4'|'shapiro'|'spectral'
    # ...（从环境变量解析填充）

# pygcm/world/params.py
@dataclass(frozen=True)
class PhysicsParams: ...
@dataclass(frozen=True)
class SpectralBands: ...
@dataclass(frozen=True)
class EcologyParams: ...
@dataclass(frozen=True)
class ParamsRegistry:
    physics: PhysicsParams
    bands: SpectralBands
    ecology: EcologyParams
```

状态容器（示意）
```python
# pygcm/world/state.py
@dataclass
class AtmosState: u: np.ndarray; v: np.ndarray; h: np.ndarray; Ta: np.ndarray; q: np.ndarray; cloud: np.ndarray
@dataclass
class OceanState: uo: np.ndarray; vo: np.ndarray; eta: np.ndarray; sst: np.ndarray
@dataclass
class SurfaceState: Ts: np.ndarray; h_ice: np.ndarray; alpha_base: np.ndarray
@dataclass
class HydroState: W_land: np.ndarray; SWE_mm: np.ndarray; C_snow: np.ndarray
@dataclass
class RoutingState: flow_accum_kgps: np.ndarray; lake_volume_kg: np.ndarray; buffers: dict
@dataclass
class EcologyState: ...  # 见 docs/15/18/16
@dataclass
class WorldState:
    atmos: AtmosState
    ocean: OceanState
    surface: SurfaceState
    hydro: HydroState
    routing: RoutingState
    ecology: EcologyState
    t_seconds: float
    schema_version: int = 1
```

世界对象（示意）
```python
class QingdaiWorld:
    def __init__(self, config: SimConfig, params: ParamsRegistry, grid: Grid, state: Optional[WorldState]=None): ...
    def step(self) -> None: ...
    def run(self, n_steps: Optional[int]=None, duration_days: Optional[float]=None) -> None: ...
    def save_state(self, path: str) -> None: ...
    def load_state(self, path: str) -> None: ...
    def diagnostics(self) -> dict: ...
```

子系统（示意接口）
```python
class Atmosphere:
    def time_step(self, world: "QingdaiWorld") -> None: ...
class Ocean:
    def step(self, world: "QingdaiWorld") -> None: ...
class Hydrology:
    def step(self, world: "QingdaiWorld") -> None: ...
class Routing:
    def step(self, world: "QingdaiWorld") -> None: ...
class Ecology:
    def step_subdaily(self, world: "QingdaiWorld") -> None: ...
    def step_daily(self, world: "QingdaiWorld") -> None: ...
```

时序（统一顺序，兼容 docs/12）
1) Forcing/光谱/反照率预处理（必要的组合器，不改变状态）  
2) Atmosphere.time_step（含 P010 反噪）  
3) Humidity/E/LH/P_cond/LH_release 步（可以由 Atmosphere 内部调用）  
4) Surface.energy（shortwave/longwave/SH/LH/seaice）  
5) Ocean.step（风应力/平流/Q_net 注入/极点一致化）  
6) Hydrology.step（相态/雪/桶/径流）  
7) Routing.step（到水文步长时执行）  
8) Ecology（subdaily/daily）  
9) Diagnostics/Plotting（按频率）  


## 5. 现有模块映射（从 → 到）

| 现有 | 新结构中的归属 | 策略 |
|---|---|---|
| dynamics.py | Atmosphere + numerics | 先 façade 代理，再逐步内收 time_step/反噪到子类 |
| energy.py/physics.py | physics（纯函数模块） + Surface | 直接复用函数；Surface 仅封装状态与组合 |
| humidity.py | humidity（纯函数） + Atmosphere 集成 | 复用；在 Atmosphere 中组织调用与写回 |
| ocean.py | Ocean + numerics | 先 façade，保持现有接口，迁移 step/风应力/极点修正 |
| hydrology.py | Hydrology + hydrology_core | 将相态/雪/桶拆纯函数，Hydrology 负责读写状态 |
| routing.py | Routing（类化已基本具备） | 统一接口，缓存/诊断对齐 |
| ecology/* | Ecology（保持结构） | Adapter 对接 world；逐步替换直连 run_simulation |
| forcing.py/orbital.py | Forcing façade + physics.solar | 组合输出 I / T_eq / 带化 I_b（docs/14） |
| topography.py | Surface 初始化/属性图生成 | 初始化加载/插值逻辑不变 |


## 6. 向后兼容策略

- 运行入口不变：`python3 -m scripts.run_simulation`；环境变量目录（docs/04）保持有效。
- 第一阶段：`scripts/run_simulation.py` 仅新增“if USE_WORLD: world = QingdaiWorld(...); world.run()”，默认仍走旧路径；以环境变量门控迁移（例如 QD_USE_OO=1）。
- 全部图像/输出目录/命名保持；diagnostics 文本格式保留关键行（用于既有分析脚本）。
- Restart 文件：保持旧字段，同时新增 group 化结构（见 §8）；使用 schema_version 做迁移。


## 7. 迁移阶段（Milestones）

- 阶段 0（façade 注入，1–2 天）
  - 新增 `pygcm/world/*`（config/params/state/world）骨架；
  - `QingdaiWorld` 内部直接调用旧脚本逻辑（adapter 直连 scripts.run_simulation 的函数块或现有类）；
  - 新增 `QD_USE_OO=1` 开关，默认关闭；冒烟测试通过。

- 阶段 1（配置/参数/状态固化，2–3 天）
  - 实装 `SimConfig/ParamsRegistry/WorldState`，由环境变量解析；
  - 旧模块读 env → 改读 world.config/params（通过适配）；
  - 输出元数据与 schema_version 注入 restart。

- 阶段 2（Forcing/Physics 纯函数化，2–4 天）
  - 将短波/长波/BL/Teq/反照率组合器从过程调用改为纯函数模块（保留原逻辑）；
  - Atmosphere/Surface/Ocean 改为仅组织调用 + 写入 state。

- 阶段 3（Atmosphere/Ocean 子系统内收，4–7 天）
  - 将 dynamics/ocean 的 time_step/step 逻辑迁入类方法；numerics 保持独立；
  - 保持 façcade 兼容；步时与闭合逐步回归。

- 阶段 4（Hydrology/Routing/Ecology 对齐，5–8 天）
  - P009/P019 算子拆分；Routing 收束为类 API；Ecology 用 Adapter 对接 world 时序；
  - 日界/子步的调度固化到 world。

- 阶段 5（JAX 互操作 + 性能与基准，5–7 天）
  - 与 projects/016 的 `xp` 后端对齐；核心算子加 `@jit`；
  - 端到端对比基准、内存占用与步时日志。

每阶段均包含“止损点”：任何异常可回滚开关（QD_USE_OO=0）；保持产线可运行。


## 8. 状态持久化与 Schema

- 容器：NetCDF（主）、NPZ（轻量 autosave）；按 group 划分：/atmos, /ocean, /surface, /hydro, /routing, /ecology, /meta
- 元数据（/meta）：
  - schema_version（int）、git_hash、created_at_utc、grid dims、config snapshot、params snapshot
  - topo source/land fraction/albedo/friction stats 对齐 README/log
- 迁移工具：io_netcdf.migrate(path, from_version→to_version)，容忍少字段/新增字段默认化
- 安全写：tmp + fsync + atomic replace；滚动备份 N 份（按 docs/15 autosave 经验）  
- 与旧格式兼容：load() 优先尝试新 schema，否则读取旧字段 + 构造缺省 group；打印黄色兼容日志。


## 9. 测试矩阵与验收标准

- 单元（pytest/numba-free）：
  - numerics：laplacian, hyperdiffuse, shapiro（谱/空间一致性）
  - physics：SW/LW/BL 等通量维度与典型输入输出范围
  - hydrology_core：相态 Sigmoid、SWE/melt 守恒
  - ecology_core：带积分能量/反射与聚合限幅
  - io：schema 读写/迁移/原子写容错

- 契约（contract tests）：
  - 子系统 API：Atmosphere/Ocean/Hydrology/Routing/Ecology 的输入/输出字段与副作用范围
  - Forcing：给定 time/params，输出 I/I_b 的确定性（对随机种子固定）

- 回归（integration/regression）：
  - 选定基线运行片段（地形相同、随机种子固定），比较：
    - 能量闭合：|⟨TOA_net⟩|、|⟨SFC_net⟩|、|⟨ATM_net⟩| < 2 W/m²（docs/06）
    - 潜热一致：⟨LH⟩（SFC）≈ ⟨LH_release⟩（ATM）
    - 水量闭合：⟨E⟩ ≈ ⟨P⟩ + ⟨R⟩（docs/09）
  - 图像 Golden：状态图/TrueColor/OceanColor 结构差异（SSIM 或结构指标）不过阈

- 性能：
  - 步时统计：OO 开关前后 Δt_step 在 ±5% 内；内存峰值不高于 +10%
  - JAX 路径可选：CPU 下降 ≥ 30%（目标），GPU/TPU 更佳

- 通过标准（本阶段验收）：
  - OO 开关启用下端到端运行稳定（日尺度、年尺度）；
  - 守恒指标满足 docs/06/08/09，且对比旧路径不劣化；
  - 脚本与 README 中示例全部可复现。


## 10. 性能与内存策略

- 禁止重复分配：状态数组在 WorldState 内创建，子系统拿视图/引用；纯函数返回写入目标数组（out 参数）。
- IO 分批：大数组按需序列化；诊断分辨率/频率可降采样。
- JAX 与 Numpy 共存：`jax_compat.xp` 统一算子；绘图/NetCDF 前显式 `np.asarray`。
- 计算图粒度：将小函数合入大核（如 combined radiation）降低 Python 开销。


## 11. 最小代码骨架（可立即落地）

```python
# pygcm/world/world.py
class QingdaiWorld:
    def __init__(self, config, params, grid, state=None):
        self.config, self.params, self.grid = config, params, grid
        self.state = state or self._alloc_initial_state()
        # façades：旧模块的代理
        self.atmos = AtmosphereFacade(self)
        self.ocean  = OceanFacade(self)
        self.hydro  = HydrologyFacade(self)
        self.route  = RoutingFacade(self)
        self.eco    = EcologyFacade(self)
        self.forcing= ForcingFacade(self)

    def step(self):
        # 1) Forcing/precompute
        self.forcing.update()
        # 2) Atmos
        self.atmos.time_step()
        # 3) Surface energy（含 SW/LW/SH/LH 与海冰）
        #    可由 atmos/physics/surface 共同完成
        # 4) Ocean
        if self.config.use_ocean:
            self.ocean.step()
        # 5) Hydrology + Routing
        self.hydro.step()
        if self.config.use_routing:
            self.route.step()
        # 6) Ecology（subdaily/daily）
        if self.config.use_ecology:
            self.eco.step_subdaily()
            if self._hits_day_boundary(): self.eco.step_daily()
        # 7) Diag/Plot
        self._maybe_plot()
        self.state.t_seconds += self.config.dt_seconds
```

> 说明：Facade 初期直接调用旧模块函数，保证接入 0 风险；随后逐步替换为新类逻辑。


## 12. 风险与对策

- 风险：接口漂移导致产线中断  
  对策：QD_USE_OO 守门；facade 先代理旧实现；阶段化回归测试。

- 风险：守恒退化  
  对策：每阶段引入“守恒回归”并将失败视为阻断；只在通过后推进。

- 风险：性能下降  
  对策：性能基准纳入每阶段退出标准；profiling 针对热点做回退或 JAX 化。

- 风险：持久化不兼容  
  对策：引入 schema_version 与迁移工具；保持旧格式读取路径与默认填充。


## 13. 时间表（建议）

- Week 1：阶段 0–1（façade + 配置/参数/状态）；冒烟回归
- Week 2：阶段 2（Forcing/Physics 纯函数化）；基线回归
- Week 3：阶段 3（Atmosphere/Ocean 内收）；性能回归
- Week 4：阶段 4（Hydrology/Routing/Ecology 对齐）；端到端年尺度试跑
- Week 5：阶段 5（JAX 互操作 + 文档/示例/README 更新）

> 所有阶段可弹性并行（Ocean 与 Ecology 可并行），但合入主干前必须通过守恒与回归。


## 14. 交付与文档更新

- 本设计文档（P020）纳入仓库；每阶段结束补齐“变更记录与默认组变化”；
- 更新 docs/12 与 README 中“开发者指南/运行 GCM”一节，标注 OO 开关与迁移状态；
- 新增开发者示例：如何用 `QingdaiWorld` 以 10 行代码跑起一次短程仿真；如何在 test harness 中构建最小 WorldState 做单元/集成测试。


## 15. 变更记录（Changelog）

- 2025‑09‑27：v1.1 可落地蓝图：对象模型、API 合约、迁移阶段、测试矩阵、持久化 schema、性能策略与骨架代码；对齐 docs/12/04/06/07/08/09/10/11/14/15/16/18。
