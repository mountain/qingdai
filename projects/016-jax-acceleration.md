# Project 016: JAX 加速——性能与可移植性提升

### 状态（2025-09-24）

  * [ ] 文档与方案定稿（本文件）
  * [ ] M1：核心计算函数（动力学/能量/湿度/海洋）JAX aot/jit 封装与性能基准
  * [ ] M2：Numpy-JAX 兼容层与类型/数组转换（`pygcm/jax_compat.py`）
  * [ ] M3：主循环（`run_simulation.py`）中的 JAX 路径切换（`QD_USE_JAX=1`）
  * [ ] M4：端到端性能与内存占用对比（JAX vs. Numpy/Scipy）
  * [ ] M5：可选：GPU/TPU 可移植性验证与文档

**关联模块与文档**

  * `docs/10-numerics-and-stability.md`（数值核心，JAX 优化的主要对象）
  * `docs/12-code-architecture-and-apis.md`（代码架构，JAX 路径的集成点）
  * `pygcm/dynamics.py`, `pygcm/energy.py`, `pygcm/humidity.py`, `pygcm/ocean.py`（将被 JAX 优化的核心模块）

-----

### 1\. 背景与动机

当前 PyGCM 完全基于 Numpy/Scipy，在 CPU 上运行良好但难以利用现代硬件加速器（GPU/TPU），且其中的循环（如 `_hyperdiffuse` 的子步）与 `scipy.ndimage` 依赖（如 `map_coordinates`）是性能瓶颈。

JAX 提供了“可组合的函数变换”（composable function transformations），包括自动微分（grad）、即时编译（jit）、自动向量化（vmap）与并行化（pmap），能以最小的代码改动显著加速科学计算，并实现跨设备（CPU/GPU/TPU）的可移植性。

本项目旨在引入 JAX 作为可选的计算后端，在不破坏现有 Numpy 接口与物理一致性的前提下，为 PyGCM 提供“即插即用”的性能加速与硬件可移植性。

-----

### 2\. 目标与范围

#### 核心目标

  * **M1 性能基准**：对 GCM 中计算最密集的函数（如拉普拉斯/超扩散、半拉氏平流、辐射通量计算等）进行 JAX jit/aot 封装，并建立性能基准。
  * **M2 兼容层**：创建 `pygcm/jax_compat.py`，提供 Numpy-JAX 数组的自动转换与设备管理，并对 `scipy.ndimage` 等 JAX 无直接对应的函数提供兼容实现（如 JAX-native 的 `map_coordinates`）。
  * **M3 JAX 路径切换**：在主循环 `run_simulation.py` 与核心模块中引入 `QD_USE_JAX=1` 开关，动态选择 JAX 或 Numpy 计算路径。
  * **M4 性能对比**：在相同配置下，端到端对比 JAX 与 Numpy 后端的“总运行时间/步”、“内存峰值占用”等指标。
  * **M5 可移植性验证**：在具备 GPU/TPU 的环境中（可选，如 Google Colab），验证 JAX 路径的可移植性与性能增益。

#### 不在本里程碑范围

  * 完全重写为 JAX-native GCM（保持 Numpy 接口与代码可读性）
  * 利用 JAX 进行自动微分或参数优化（可作为 P017 等后续项目）
  * 分布式并行（pmap）（首版聚焦单设备加速 jit/aot）

-----

### 3\. 架构与集成策略

#### “鸭子类型”与条件导入

  * 利用 JAX `numpy` (`jax.numpy` as `jnp`) 与 Numpy `numpy` (`numpy` as `np`) 接口的高度相似性，实现“鸭子类型”（duck typing）兼容。
  * 在 `pygcm/jax_compat.py` 中：
      * 根据 `QD_USE_JAX` 环境变量决定 `jnp` 或 `np` 作为计算后端 `xp`
      * 提供设备管理函数 `get_device()`, `to_device()`
      * 为 `scipy.ndimage.map_coordinates` 等提供 JAX 实现或包装器
  * 在核心模块（`dynamics.py`, `energy.py` 等）中：
      * `from .jax_compat import xp, device, jax_map_coordinates`
      * 数值计算函数（如 `_laplacian_sphere`）以 `xp` 为数组库，并用 `@jax.jit` 装饰
      * I/O 与 matplotlib 可视化仍使用原生 Numpy

#### JAX aot (Ahead-of-Time Compilation)

  * 对于不依赖于循环中动态形状的核心计算核（如 `_laplacian_sphere`），可使用 aot 预编译，降低首次调用的 jit 延迟。
  * `jax.scipy.ndimage.map_coordinates` 是 JAX-native 的 `map_coordinates` 实现，可直接替换。

-----

### 4\. 任务拆解与里程碑

  * **M1 JAX 封装与性能基准**
      * [ ] 在 `requirements.txt` / `pyproject.toml` 中添加 `jax` 与 `jaxlib`
      * [ ] 对 `_laplacian_sphere`, `_hyperdiffuse`, `_advect` 等核心函数创建 JAX 版本并用 `@jax.jit` 装饰
      * [ ] 编写孤立的性能基准脚本（`scripts/benchmark_jax.py`），对比 JAX vs Numpy 版本在不同分辨率下的耗时
  * **M2 Numpy-JAX 兼容层**
      * [ ] 创建 `pygcm/jax_compat.py`
      * [ ] 实现 `xp` 的条件导入与设备管理
      * [ ] 实现或包装 `jax.scipy.ndimage.map_coordinates` 等 Scipy 兼容函数
  * **M3 主循环 JAX 路径切换**
      * [ ] 在 `run_simulation.py` 中引入 `QD_USE_JAX` 开关
      * [ ] 在 `dynamics.py`, `energy.py` 等模块中，将 `np` 替换为 `xp`，`map_coordinates` 替换为 `jax_map_coordinates`
      * [ ] 确保 I/O（NetCDF）、matplotlib 可视化前，JAX 数组能正确转回 Numpy (`np.asarray(jax_array)`)
  * **M4 端到端性能对比**
      * [ ] 运行 spin-up.sh 或标准短程模拟，分别在 `QD_USE_JAX=0` 与 `QD_USE_JAX=1` 下记录总时长与内存占用
      * [ ] 输出对比报告图表
  * **M5 可移植性验证（可选）**
      * [ ] 在 Google Colab GPU/TPU 环境下运行 JAX 路径，记录性能增益
      * [ ] 在 README 与 `docs/10-numerics-and-stability.md` 中补充 JAX GPU/TPU 运行指南

-----

### 5\. 环境变量与接口（新增）

  * **`QD_USE_JAX`**（默认 0）：=1 时启用 JAX 作为计算后端
  * **`QD_JAX_PLATFORM`**（可选；`cpu|gpu|tpu`）：强制 JAX 使用特定设备平台
  * **`QD_JAX_AOT`**（默认 1）：启用 AOT 预编译（若实现）

-----

### 6\. 风险与注意事项

  * **“纯函数”约束**：JAX jit 要求被装饰函数为“纯函数”（pure function），即无副作用、不修改输入参数、不依赖全局状态。需对现有代码做少量重构以满足此约束（如将类方法改为静态方法，显式传递状态）。
  * **首次调用延迟（jit compilation lag）**：JAX 函数首次调用时会触发编译，导致耗时较长。可用 AOT 或在模拟开始前做一次“热身”调用（warm-up call）缓解。
  * **`scipy.ndimage` 兼容性**：`jax.scipy.ndimage.map_coordinates` 与 `scipy` 版本在 `mode` / `cval` 等参数上可能存在细微差异，需验证一致性。`convolve` 等其它函数也需找到 JAX 等价实现。
  * **内存管理**：JAX 数组默认在设备上分配，与 Numpy 的 CPU 内存交互（如 I/O、绘图）需显式转换，避免不必要的拷贝开销。

-----

### 7\. 验收标准

  * **功能性**：`QD_USE_JAX=1` 时，模型可稳定运行并产出与 Numpy 版本在数值容差内（\~1e-6）一致的结果。
  * **性能**：CPU 上 JAX jit 路径相比 Numpy/Scipy 平均每步耗时降低 ≥ 30-50%（目标值，视具体函数与分辨率而定）。
  * **可移植性**：GPU/TPU 环境下 JAX 路径性能显著优于 CPU（数量级提升）。
  * **代码质量**：JAX 路径通过条件导入与兼容层实现，对原有代码的侵入性最小化。

-----

### 8\. 运行示例

**A) 默认 Numpy 后端**

```bash
# (不设置 QD_USE_JAX 或设为 0)
python3 -m scripts.run_simulation
```

**B) 启用 JAX 加速（CPU）**

```bash
export QD_USE_JAX=1
python3 -m scripts.run_simulation
```

**C) 启用 JAX 加速（GPU，若环境支持）**

```bash
export QD_USE_JAX=1
export QD_JAX_PLATFORM=gpu
python3 -m scripts.run_simulation
```