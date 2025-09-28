# Project 021：双缓冲状态引擎（Double Buffering State Engine）

状态（2025‑09‑27）
- [x] 方案定稿（本文件）
- [ ] M1：核心抽象 DoubleBufferingArray 与单元测试
- [ ] M2：NumPy/JAX 后端验证与基准测试（前置）＋ WorldState 双缓冲化与子系统读/写契约改造
- [ ] M3：与 P020（OO 重构）集成，按步原子切换
- [ ] M4：回归与验收（守恒、不回归、性能）
- [ ] M5：文档与示例更新

交叉引用
- 架构与 API：docs/12-code-architecture-and-apis.md
- 面向对象重构：projects/020-oo-refactor.md
- 数值核心/稳定：docs/10-numerics-and-stability.md
- 脚本与主循环：scripts/run_simulation.py（集成位点）

---

## 0. 背景与动机

现有主循环在“读—改—写”过程中容易出现以下工程性问题：
- 状态原地修改导致跨模块“脏读”，难以追踪与调试；
- 为了保证原子性而创建“新状态对象”的策略带来大量数组分配与复制开销；
- JAX 路径（@jit）要求纯函数/无副作用，难以与“原地改写”混用。

目标：将“双缓冲（Double Buffering）”由“顶层流程约定”内化为“底层数组抽象”的自洽能力，让物理模块“自然”地在读/写分离的缓冲上计算，步末只需 O(1) 的指针翻转即可实现原子更新。

本项目提供一个轻量、后端无关（NumPy/JAX）的核心抽象 DoubleBufferingArray（简称 DBA），并将其无缝嵌入 P020 的 `WorldState`，实现“读缓冲只读、写缓冲只写、步末交换”的统一契约。

---

## 1. 目标与可交付

- 核心抽象（M1）
  - DoubleBufferingArray：封装两块同形状数组，提供 `.read`（只读视图）、`.write`（写视图）与 `swap()`（O(1) 交换）；
  - 兼容 NumPy/JAX 的数组后端（通过 `pygcm/jax_compat.py` 的 `xp` 抽象）；
  - 不破坏现有下游调用：支持基础索引、`__array__`（避免误用时崩溃）。

- 世界状态改造（M2）
  - 将 `WorldState` 中的各网格场（u/v/h/Ts/Ta/q/cloud/SST/uo/vo/η/...）替换为 DBA；
  - 提供 `WorldState.swap_all()` 原子翻转所有 DBA（步末调用）。

- 契约与集成（M3）
  - 明确所有子系统（Atmosphere/Ocean/Surface/Hydrology/Routing/Ecology）的读/写契约：
    - 只能从 `.read` 读取前一状态；
    - 只能往 `.write` 写入下一状态；
  - 主循环在所有子系统完成写入后再调用 `swap_all()`，实现步级原子性。

- 后端与性能（M4）
  - NumPy 与 JAX 双后端验证（`xp = numpy | jax.numpy`）；
  - 基准：相对“新建状态复制”方案，内存峰值下降、步时减少，目标：
    - 分配/复制次数 → ~0；
    - 端到端步时减少 ≥ 10–20%（视分辨率与模块混合而定）。

- 验收（M5）
  - 功能：结果与“无 DBA 参考实现”在容差内一致；
  - 守恒：能量/水量长期诊断不劣化（docs/06/08/09 标准）；
  - 性能：达到上文性能目标；JAX 路径可正常运行（见 §6 注意事项）。

---

## 2. 设计总览

双缓冲理念：每个状态场持有两块相同形状/类型的缓冲区（A 与 B）。在时间步 n：
- 模块只读 A，所有计算结果写入 B；
- 步末执行“指针翻转”：下一步读 B，写 A；
- 不复制数据、不分配新数组，原子性由翻转保证。

与 P020 的关系：P021 是 P020 的“状态演化策略实现”。P020 规定 World 对外暴露“值语义”的 step 行为；P021 在内部用双缓冲落地“伪不可变”策略（详见 P020 §1A）。

---

## 3. 核心抽象：DoubleBufferingArray（DBA）

文件：`pygcm/numerics/double_buffer.py`（建议）

关键接口（Python 伪代码）：
```python
from typing import Tuple, Optional
import numpy as np
from pygcm.jax_compat import xp  # xp = np or jnp

class DoubleBufferingArray:
    def __init__(self, shape: Tuple[int, ...], dtype=np.float64, initial_value=0.0):
        self._a = xp.full(shape, initial_value, dtype=dtype)
        self._b = xp.full(shape, initial_value, dtype=dtype)
        self._read_idx = 0  # 0=>_a read, _b write; 1=>_b read, _a write

    @property
    def read(self):
        return self._a if self._read_idx == 0 else self._b

    @property
    def write(self):
        return self._b if self._read_idx == 0 else self._a

    def swap(self):
        # O(1) 指针翻转
        self._read_idx ^= 1

    # 便利属性
    @property
    def shape(self): return self.read.shape
    @property
    def dtype(self): return self.read.dtype

    # 仅用于与少量 NumPy API 的互操作（谨慎使用）
    def __array__(self, dtype=None):
        arr = self.read
        return arr.astype(dtype) if dtype is not None else arr

    def __getitem__(self, key):
        return self.read[key]
```

实现要点
- `.write[:] = value` 赋值必须使用切片写入，避免替换底层数组对象；
- 不提供 `.write = ...` 属性替换，防止破坏内存占用优势；
- 不推荐将 DBA 直接传入 JAX @jit 编译的函数（见 §6），在 jitted 路径内应传入 `DBA.read`。

可选增强（后续）
- `from_array(a, like='read'|'write')`：用现有数组初始化 A/B；
- `zero_write()`：将写缓冲清零；
- `map_(fn_read_to_write)`：以函数封装 read→write 的映射（便于链式算子）。

---

## 4. WorldState 双缓冲化与子系统契约

WorldState 示例（概念）：
```python
# pygcm/world/state.py
from dataclasses import dataclass
from pygcm.numerics.double_buffer import DoubleBufferingArray as DBA

@dataclass
class AtmosState:
    u: DBA; v: DBA; h: DBA; Ta: DBA; q: DBA; cloud: DBA

@dataclass
class OceanState:
    uo: DBA; vo: DBA; eta: DBA; sst: DBA

@dataclass
class SurfaceState:
    Ts: DBA; h_ice: DBA; # ... 其它地表/冰雪/掩膜派生量可按需 DBA 化

@dataclass
class HydroState:
    W_land: DBA; SWE: DBA; # ... 等

@dataclass
class WorldState:
    atmos: AtmosState
    ocean: OceanState
    surface: SurfaceState
    hydro: HydroState
    t_seconds: float

    def swap_all(self):
        for sub in (self.atmos, self.ocean, self.surface, self.hydro):
            for _, value in vars(sub).items():
                if isinstance(value, DBA):
                    value.swap()
```

子系统读/写契约（强约束）
- 读取：仅从 `.read` 读取前一状态；
- 写入：仅向 `.write` 写入下一状态；
- 禁止在同一步中读取 `.write`、或写入 `.read`；
- 执行顺序：由世界对象统一编排；步末仅调用一次 `swap_all()`。

主循环（简化示例）：
```python
class QingdaiWorld:
    def step(self):
        # Forcing/诊断准备（只读）
        self.forcing.update(self.state)

        # 子系统更新（读 read，写 write）
        self.atmos.time_step(self.state, self.params, self.config)
        if self.config.use_ocean:
            self.ocean.step(self.state, self.params, self.config)
        self.hydrology.step(self.state, ...)  # SWE/桶/径流
        self.routing.step(self.state, ...)    # 到达水文步长再执行

        # 步末原子翻转
        self.state.swap_all()
        self.state.t_seconds += self.config.dt_seconds
        return self.state
```

---

## 5. 与 NumPy/JAX 的互操作策略

- 后端抽象：沿用 `pygcm/jax_compat.py` 提供的 `xp`（numpy | jax.numpy），DBA 内部仅调用 `xp.*`；
- NumPy：常规函数（如 `np.sin(dba)`）会通过 `__array__` 获取 `.read`，但会触发拷贝/回传，建议在热点路径显式用 `arr = dba.read`；
- JAX：
  - 约定将随 tests 逐步固化：推荐不要把 DBA 直接传给 @jit 函数，请传入 `dba.read`（`DeviceArray`/`Array`）；如存在必要特例，将通过测试用例与文档明确；
  - 所有 jitted 内核的输出应写回 `.write`（由外层 Python 调用完成）。

建议规范（贯穿代码评审；约定将随 tests 逐步固化）
- jitted 核心纯函数签名形如：`fn(out, *ins)` 或返回新数组，由外层 `.write[:] = fn(...)`；
- 子系统方法外层只组织 `.read` 与 `.write`，不在内部改变指针（pointer flip 始终在世界层统一执行）。

---

## 5A：Magic Methods Policy（非 JIT）

为在非 JIT 环境下简化 DBA 使用、减少显式 `.read/.write` 样板，同时保持读/写约束一致性，制定如下最小策略：

- 读取/写入基准
  - 读取一律来自 `read`；写入一律落到 `write`；步末 `swap()` 才“生效”到下一步。
- 建议实现的魔术方法
  - `__getitem__(key)` → `self.read[key]`（只读视图/切片）
  - `__setitem__(key, value)` → `self.write[key] = np.asarray(value)`（定向写 write，不触碰 read）
  - `__array__(dtype=None)` → 返回 `self.read`（必要时 astype），用于 NumPy 透明互操作
  - `__array_ufunc__(ufunc, method, inputs, kwargs)`（最小规则）
    - 无 `out` 时：将 DBA 输入替换为 `.read` 调用底层 ufunc，返回 ndarray（无副作用）
    - 有 `out` 且为 DBA 时：将该 `out` 替换为其 `.write`，调用底层 ufunc（定向写 write）
- 不建议实现
  - 原地运算魔术（`__iadd__` 等）——避免“立即原地修改”与 DBA 的“写 write + 需 swap 生效”语义冲突；若需写入请用 `np.add(a, b, out=dba)` 或 `dba.write[:] = ...`
- 边界与保护
  - 热点路径避免隐式 `__array__` 触发的拷贝，改用显式 `arr = dba.read`
  - 自别名写入：`dba[...] = dba` 建议禁止（raise），避免源/目标别名导致未定义行为
- 与 JIT 的关系
  - JIT 图内不使用自定义对象/魔术；jitted 核心仅接受底层数组（`dba.read`），输出由外层写回 `dba.write`；该约定将随 tests 逐步固化

测试覆盖（见 M1）
- `__getitem__/__setitem__`：写 write 不影响 read；`swap()` 后可见
- `__array__`：`np.asarray(dba)` / `np.sin(dba)` 仅从 read 取值
- `__array_ufunc__`：`np.add(dba, 1)` 返回 ndarray；`np.add(dba, 1, out=dba)` 写 write
- 自别名写入：`dba[...] = dba` 的期望行为（默认 raise）

---

## 6. 实现计划与里程碑

M1：核心抽象与测试（1–2 天）
- [ ] 新增 `pygcm/numerics/double_buffer.py`；
- [ ] 单元测试
  - 读/写分离与 swap 行为；
  - 切片写入保持对象不替换；
  - __getitem__/__setitem__：读 read / 写 write 语义与 swap 生效；
  - __array__：`np.asarray(dba)` / `np.sin(dba)` 仅从 read 取值；
  - __array_ufunc__：`np.add(dba, 1)` 返回 ndarray；`np.add(dba, 1, out=dba)` 写 write；
  - 别名写入保护：`dba[...] = dba`（明确行为，推荐 raise）；
  - 与 NumPy 基本 ufunc 互操作（`np.sum`, `np.asarray` 等）；
  - JAX 后端 smoke：`jnp.add(dba.read, ...)` 路径。

M2：后端验证与基准（前置，2–4 天）＋ 最小子系统接入（2–3 天，可并行）
- [ ] NumPy/JAX 双后端一致性回归（容差内）；
- [ ] 基准：与“每步新建状态复制”方案对比，记录内存峰值与步时；
- [ ] CI 增加 smoke/小基准，避免性能回退；
- [ ] 在世界对象 `WorldState` 中替换关键字段为 DBA；
- [ ] 改造 Atmosphere/Ocean/Hydrology 中最小热点路径以遵守读/写契约；
- [ ] 提供 `swap_all()`，主循环步末调用。

M3：与 P020 集成（2–3 天）
- [ ] 在 `QingdaiWorld` façade 中接入 DBA（`QD_USE_OO=1` 时启用）；
- [ ] 兼容旧路径（`QD_USE_OO=0`），不破坏现有运行；
- [ ] 按模块分阶段切换，保留止损开关。

M4：回归与验收（1–2 天）
- [ ] 守恒：TOA/SFC/ATM 与 E–P–R 长期平均在阈值内；
- [ ] 性能：达成 §1 的目标；

M5：文档与示例更新
- [ ] 更新本文件/README/开发者指南与示例；将 JIT 约定在单测中逐步固化并补充到文档。

---

## 7. 示例：在子系统中使用 DBA

以大气动量更新为例（简化）：
```python
def time_step(self, state, params, config):
    u = state.atmos.u.read
    v = state.atmos.v.read
    h = state.atmos.h.read

    # 计算倾向（示意）
    du, dv = compute_momentum_tendency(u, v, h, params, config)  # 纯函数

    # 写入下一步
    state.atmos.u.write[:] = u + du * config.dt_seconds
    state.atmos.v.write[:] = v + dv * config.dt_seconds
```

海表温度（SST）半拉氏平流（示意）：
```python
def advect_sst(self, state, dt):
    sst = state.ocean.sst.read
    uo = state.ocean.uo.read
    vo = state.ocean.vo.read
    sst_next = semi_lagrangian_advect(sst, uo, vo, dt)  # 纯函数（NumPy/JAX）
    state.ocean.sst.write[:] = sst_next
```

---

## 8. 测试与验收标准

功能一致性
- 与“无 DBA 参考路径”在相同随机种子/参数下，关键诊断（全球均值/能量/水量闭合）在容差内一致；
- 视觉回归：状态图/TrueColor/河网层无结构性差异（允许微小数值扰动）。

数值与稳定
- 禁止在一步内读取 `.write` 或写 `.read`（加入静态检查/审查清单）；
- JAX 路径：所有 jitted 内核均不接受 DBA 参数，仅处理底层数组。

性能
- 分配/复制次数接近 0（热点路径）；
- 端到端步时降低（基准脚本 `scripts/benchmark_jax.py` 可扩展比对 DBA/非 DBA）。

CI 与回归
- 新增 `tests/test_double_buffering.py`：
  - 读写隔离、swap、切片写入、NumPy/JAX smoke；
- 将一条“短程烟测 + 诊断阈值”纳入 CI（Ubuntu/macOS, Py3.11–3.13）。

---

## 9. 风险与对策

- 误把 DBA 传入 @jit 核心导致 host→device 复制/trace 失败  
  对策：代码规范 + 评审 + linters + tests；jitted 路径统一以 `arr = dba.read` 为入口。

- 在子系统中“提前 swap”或跨模块混合读写  
  对策：禁止在子系统内 swap；仅在世界步末 `swap_all()`；加入运行期断言（可选 DEBUG）。

- `__array__` 触发隐式拷贝带来性能异常  
  对策：在热点路径不使用隐式转换；以 `.read` 显式传递；基准守护发现回退。

- 与旧路径并存期间的维护负担  
  对策：P020 的 façade＝API 合约（见 020 §3.5）；按模块逐步切换，保留止损开关。

---

## 10. 运行与配置

- 默认随 P020 路径启用：`export QD_USE_OO=1`（OO façade 生效）；
- 若需对比/回退：`export QD_USE_OO=0`（遗留路径）；
- 无需新增环境变量；DBA 完全在内部抽象层生效。

---

## 11. 变更记录（Changelog）

- 2025‑09‑27：v1 方案定稿：定义 DBA 抽象、WorldState 双缓冲化、子系统读/写契约、NumPy/JAX 互操作策略、基准与验收标准。

---

## 附录 A：FAQ

Q：能否让 `np.sin(dba)` 直接在不拷贝的情况下工作？  
A：`__array__` 会产生一个数组视图/拷贝语义，具体取决于后端与调用；为避免隐藏成本，建议在热点路径显式使用 `arr = dba.read` 后交给数值函数处理。

Q：写入 `.write` 是否会与 `.read` 别名重叠？  
A：DBA 内部使用两块独立缓冲（A/B），只在指针翻转时切换角色，不存在别名重叠。

Q：JAX 下如何保证纯函数？  
A：jitted 内核仅接受原生 `DeviceArray/Array`；DBA 作为“Python 侧资源管理器”停留在图外；输出通过外层 `.write[:] = ...` 写回。

Q：是否支持并行/多线程？  
A：DBA 指针翻转不是线程安全原语；当前设计假定单线程步进（与现有主循环一致）。如需多线程/并行，将在更高层通过任务图/分域策略处理。

---

## 附录 B：所有权（Ownership）与耦合器写禁（No-Write Couplers）

本附录在 P021 双缓冲写入契约之上，引入“所有权（ownership）驱动的写入权限模型”，明确 Engine 与 Coupler 的责任边界，避免跨模块的隐式状态写入，提升系统的可理解性与可测性。

### B.1 核心原则（Single Owner, Write by Owner Only）

- 原则：WorldState 的每个“子状态（substate/fields）”都有且仅有一个 Owner（所有者）。只有 Owner 才被允许写入它所拥有的那部分状态。
- Owner 的写入必须遵循 P021 的双缓冲契约：只写入 next（write buffer），不可写入 current（read buffer）。

典型所有权映射（示意，实际映射在实现中集中定义）：
- AtmosphereEngine 拥有：大气动力与热力相关字段（u, v, h, Ta, q, cloud …）
- OceanEngine 拥有：海洋相关字段（uo, vo, η, SST …）
- HydrologyEngine/WaterCycle 拥有：陆面桶、SWE、相态等（W_land, SWE …）
- EcologyEngine 拥有：生态自身状态（LAI/biomass/genes/seed_bank 等）
- Routing（异步路由）拥有：其内部累加/诊断缓冲，对 WorldState 仅追加/回填约定字段

注：任何“跨域派生量”（如可视化缓存、诊断快照）不得反向写入核心物理状态。

### B.2 耦合器（Coupler）是纯计算者（只读→端口→返回）

- 定义：Coupler 只负责读 current_state，计算跨子系统的相互作用量（通量/应力/源汇等），并以 Ports（端口）形式返回结果。
- 禁止：Coupler 直接写入 WorldState（无论 read 还是 write buffer）。Coupler 不是状态所有者。

接口契约（纯函数）：
```
Coupler.compute(ports_in: Ports) -> ports_out: Ports
```

Ports 是对跨模块物理量的显式打包，具有类型安全与可组合性（示例见 B.5）。

### B.3 编排者（Orchestrator/QingdaiWorld）的职责

- 收集：从 current_state 组装 Coupler 所需的 ports_in（来自各 Owner 的只读视图）。
- 调用：按顺序或图（DAG）调用若干 Coupler，聚合其 ports_out。
- 分发：将聚合后的 forcings/ports 传入相应 Owner 的 Engine.step。
- 严禁：编排层绕过 Engine 直接写 WorldState 的 owned 区域；若确需写入非 Owned 辅助区，必须在“临时缓冲/诊断域”，不得污染核心子状态。

### B.4 与 P021 双缓冲的结合

- 写入面：只有 Owner 的 Engine.step 可以写 next（write buffer）中的 owned 字段；其余任何模块（包括 Coupler/Orchestrator）不得触碰 owned 字段的 write buffer。
- 读取面：Coupler 与 Engine 均从 current（read buffer）读取依赖量（自身/他域）。
- 步末：由世界对象统一执行一次 swap_all()，完成原子切换。

该模式将“写权限”与“缓冲角色”同时约束，最大限度避免脏写/脏读。

### B.5 API 契约（建议签名）

- Coupler（只读纯函数）：
```
class WindStressCoupler:
    def compute(self, ports_in: "PortsIn") -> "PortsOut":
        ...
```

- Engine（唯一写者）：
```
class OceanEngine:
    def step(
        self,
        read_state: WorldState,        # 只读（current buffer）
        write_state: WorldState,       # 只写（next buffer）——仅写 owned 字段
        forcings: "Ports"              # 来自一个或多个耦合器的聚合外强
    ) -> None:
        ...
```

- Ports（示例，具体以模块落地为准）：
```
@dataclass
class PortsIn:
    # 读取所需的最小集（只读视图/复制安全结构）
    u_atm: ArrayLike
    v_atm: ArrayLike
    uo: ArrayLike
    vo: ArrayLike
    land_mask: ArrayLike
    # … 可按耦合器扩展

@dataclass
class PortsOut:
    wind_stress: ArrayLike | None = None
    heat_flux: ArrayLike | None = None
    moisture_flux: ArrayLike | None = None
    # … 按需扩展（保持可空、可组合）
```

### B.6 示例：风应力耦合与海洋写入

- 读取：WindStressCoupler.compute() 从 current_state 读取 u_atm, v_atm 与海洋相对风信息（必要时含 uo, vo），形成 ports_out.wind_stress。
- 传递：Orchestrator 将 wind_stress 作为 forcings 传给 OceanEngine.step()。
- 写入：OceanEngine 作为 ocean 子状态的所有者，依据 wind_stress 仅写 write_state.ocean.*（如 uo_next, vo_next, η_next, SST_next）。

### B.7 约束与工程落地（建议）

- 所有权映射集中化：在单一源（例如 world/ownership.py）中定义“字段 → Owner 引擎”的映射表，并在 DEBUG 模式启用运行期断言：
  - 断言只有 Owner 在写入其 owned 字段的 write buffer；
  - 断言 Coupler 未触碰 WorldState（可通过接口不暴露 write_state 给 Coupler，辅以只读包装）。
- 静态/代码审查：
  - Coupler 目录禁止 import WorldState 写接口；仅允许使用只读型 DTO（Ports）。
  - 引入“写域门面（WriteFacade）”类型，在 Engine.step 之外不可见。
- 单元/契约测试：
  - 针对每个 Engine：验证只写 owned 字段与 write buffer；
  - 针对每个 Coupler：验证 compute 纯函数性（相同输入→相同输出；无副作用）。
- DBA（DoubleBufferingArray）层保护（可选）：
  - 在开发/测试构建中，给 read 视图包装“不可写 ndarray（writeable=False）”，在 NumPy 路径及时抛错；
  - 提供“写入探针”统计每步各字段的写入者，异常直接 fail-fast。

### B.8 与 P020/P021 的协同关系

- P020（OO 重构）定义对象边界与依赖注入；本附录给出“谁能写”的强约束，作为 façade 的 API 合约组成部分。
- P021（双缓冲）提供“写入时空域”的技术保障；与 Ownership 组合后，形成“字段级写者”与“缓冲级写面”的双重护栏。

### B.9 迁移与兼容

- 迁移期允许“临时适配层”将 legacy 写入改造为“Ports→Engine.step”路径；完成模块内收后移除适配。
- 若历史代码中存在 Coupler 写状态的路径，应立即改造为返回 Ports；任何直写行为必须删除或封装到对应 Engine 中。

---
