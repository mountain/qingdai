# Project 005: 将地形（Topography）数据融入 GCM 模型

本项目在 P004 的基础上，把生成的地形与海陆属性数据（elevation / land_mask / base_albedo / friction）正式接入 GCM 主循环，作为动力与辐射的边界条件与参数化输入。该集成工作将替代旧的“矩形大陆”方案，并为后续的地形降水（orographic precipitation）与更真实的表面过程奠定基础。

## 0. 状态（2025-09-20）

- [x] P004 已完成里程碑 1：生成 L1+L3 地形、海平面自适应，导出标准化 NetCDF，提供基本可视化脚本。
- [x] P005 本文档（设计与任务拆解）
- [x] P005 代码实现：加载 NetCDF → 传入 GCM → 动态反照率使用 base_albedo → 运行示例
- [x] P005 扩展（可选）：地形降水（迎风/背风、上坡风抬升）参数化
- [x] P005 文档与示例更新（README、run_simulation 说明）

## 1. 目标

- 从 `data/*.nc` 读取由 P004 生成的地形与表面属性，并在 GCM 初始化时加载到模型：
  - `elevation`（m），`land_mask`（0/1）
  - `base_albedo`（无量纲），`friction`（s^-1）
- 在动力核心中使用 `friction_map`（已支持）。
- 在辐射/反照率计算中使用 `base_albedo` 替代过去固定的 `alpha_water` 标量。
- 提供配置开关与运行脚本，保证“读取外部 NetCDF”与“退回生成模式”均可工作。
- 预留 orographic precipitation 的参数化入口（基于地形坡度与迎风向）。

## 2. 数据与接口

### 2.1 NetCDF 数据字段（来自 P004）

- `lat`（degrees_north），`lon`（degrees_east）
- `elevation`（m）
- `land_mask`（0/1）
- `base_albedo`（1）
- `friction`（s^-1）
- 全局属性：`sea_level_m`、`target_land_fraction`、planet constants（用于可追溯性）

### 2.2 新增/调整 API（建议）

- 在 `pygcm/topography.py` 中新增：
  - `load_topography_from_netcdf(path, grid, *, regrid="auto") -> (elevation, land_mask, base_albedo, friction)`
    - 若文件分辨率与 `grid` 一致，直接读取。
    - 若不一致，提供简单插值（优先最近邻/双线性；经度按周期处理）。
    - 对维度/数值进行基本校验（NaN、范围）。

- 在 `pygcm/physics.py` 中扩展动态反照率接口（保持向后兼容）：
  - 现函数 `calculate_dynamic_albedo(cloud_cover, T_s, alpha_water, alpha_ice, alpha_cloud)`
    - 支持 `alpha_water` 既可为标量也可为 2D 数组（base_albedo_map）。当传入数组时，逐点使用。
  - 或新增函数：
    - `calculate_dynamic_albedo_with_base(base_albedo_map, T_s, alpha_ice, alpha_cloud)`

- `scripts/run_simulation.py` 启动参数/环境变量：
  - `QD_TOPO_NC`：如提供，则优先从该 NetCDF 读取 topography 与表面参数；否则退回旧生成流程（或直接使用 P004 生成器）。
  - `QD_USE_TOPO_ALBEDO`（默认 1）：使用 `base_albedo` 作为地表/海洋基反照率的底图并与云反照率混合。
  - `QD_OROG`（默认 0）：是否启用地形降水增强（可选，见 §4）。
  - `QD_OROG_K`（默认 5e-4 – 1e-3 量级，经验值）：地形抬升强度系数。

## 3. 集成点与主流程

### 3.1 run_simulation 初始化阶段（示意）

```python
# 1. Grid
grid = SphericalGrid(n_lat=..., n_lon=...)

# 2. Topography & Surface Properties
topo_nc = os.getenv("QD_TOPO_NC")
if topo_nc and os.path.exists(topo_nc):
    elevation, land_mask, base_albedo_map, friction_map = load_topography_from_netcdf(topo_nc, grid)
else:
    # Fallback (旧逻辑或直接 P004 生成器)
    land_mask = create_land_sea_mask(grid)  # 或：从生成器返回 elevation + mask
    base_albedo_map, friction_map = generate_base_properties(land_mask)
    elevation = None  # 若未生成 elevation，可置 None

# 3. Dynamics core
gcm = SpectralModel(grid, friction_map, H=..., tau_rad=..., greenhouse_factor=...)

# 4. Physics parameters
USE_TOPO_ALBEDO = int(os.getenv("QD_USE_TOPO_ALBEDO", "1")) == 1
alpha_ice = ...
alpha_cloud = ...

# 在主循环中计算 albedo：
if USE_TOPO_ALBEDO:
    # 允许 base_albedo_map 替代 alpha_water 标量
    albedo = calculate_dynamic_albedo(gcm.cloud_cover, gcm.T_s, base_albedo_map, alpha_ice, alpha_cloud)
else:
    albedo = calculate_dynamic_albedo(gcm.cloud_cover, gcm.T_s, alpha_water=0.1, alpha_ice=0.6, alpha_cloud=0.5)
```

要点：
- `friction_map` 已在 `SpectralModel` 的 `time_step` 中用于地表摩擦。
- `base_albedo_map` 通过 `calculate_dynamic_albedo` 参与辐射反照率。
- `elevation` 暂不强制使用（用于 §4 地形降水与可视化等）。

### 3.2 维度与插值（简要）

- NetCDF → Grid 对齐：若维度不一致，按经度周期性双线性插值；高纬处理需注意 cos(lat) 缩放。
- 插值后的 `land_mask` 建议用最近邻，以保持二值性。
- 插值函数建议统一收敛到工具函数中，避免分散实现。

## 4. （可选）地形降水参数化（里程碑 2 的先导）

目标：在现有的“动力-降水”基础上，引入地形抬升效应，增强迎风坡降水，减弱背风坡。

- 简化方案：在诊断出 `precip_base` 后，乘以地形因子
  ```
  n̂ = ∇H / |∇H|   （坡向单位向量，缺省时 ∇H≈0 则置因子=1）
  U = (u, v)       （近地面风）
  uplift = max(0, U · n̂)    # 迎风分量（上坡风为正）
  factor = 1 + k_orog * uplift
  precip = precip_base * factor
  ```
- 计算注意：
  - ∇H 用球面近似的经纬度梯度，按经纬度实际距离（a·cosφ·Δλ、a·Δφ）换算。
  - `k_orog` 规模化系数经验定标；避免过大以致数值爆裂。
  - 可对 `uplift` 做平滑，或对 factor 做上限截断（如 ≤ 2.0）。

- 可通过 `QD_OROG=1` 开启，`QD_OROG_K` 调整强度。

## 5. 任务分解

- T1. I/O：实现 `load_topography_from_netcdf(path, grid, regrid="auto")`
  - 最近邻/双线性插值（lon 周期），数据校验与日志打印。
- T2. 反照率：扩展 `calculate_dynamic_albedo` 接受 2D `alpha_water`（base_albedo）
  - 或新增 `_with_base` 版本。保留旧接口以确保兼容。
- T3. 主脚本：`scripts/run_simulation.py`
  - 支持 `QD_TOPO_NC`、`QD_USE_TOPO_ALBEDO` 环境变量。
  - 在可视化图中（如 plot_state / plot_true_color）可选叠加海岸线/地形等值线（提升辨识）。
- T4. 可选：`orographic_enhancement` 工具函数与接入点
  - 实现 `compute_orographic_factor(grid, elevation, u, v, k_orog) -> factor`
  - 在 `diagnose_precipitation` 之后应用（避免循环依赖），或封装成单独模块。
- T5. 文档与示例
  - 更新 README：新增 “使用外部 topography NetCDF 运行 GCM”。
  - 增加示例命令集：加载 data 下最新 nc，运行若干行星日并输出图像。

## 6. 接口与数据一致性约束

- `land_mask` 的加权面积分（cosφ 权重）应与全局属性 `target_land_fraction` 一致（±0.02）。
- `friction` 与 `base_albedo` 与 `land_mask` 的空间分布应物理合理（海洋低反照率/低摩擦，山地摩擦更大）。
- GCM 初始化后首次出图时，需在日志中打印：
  - “Topo source”（NetCDF 路径或生成模式）
  - “Land fraction (target / achieved)”
  - “Friction/Albedo stats（min/mean/max）”
  - “Elevation stats（min/mean/max）”（如可用）

## 7. 测试与验收标准

- 功能性：
  - [ ] 设置 `QD_TOPO_NC` 后能成功读取并运行 1–2 个大周期间隔出图；
  - [ ] 未设置时退回生成模式仍可运行。
- 一致性：
  - [ ] 模型网格与数据网格不一致时，插值结果无索引越界，海岸线无明显断裂。
- 物理合理性：
  - [ ] 风场、温度场稳定无爆裂；
  - [ ] 降水分布与海陆/地形具备基本一致性（陆地有更多降水中心，海上仍可有 ITCZ 风格带状特征）。
- 输出：
  - [ ] `output/` 中生成包含海岸线/陆架等可辨别细节的图像；
  - [ ] 日志中包含加载信息与地形统计。

## 8. 运行示例

- 使用 data 下最新 topography：
  ```
  export QD_TOPO_NC=$(ls -t data/*.nc | head -n1)
  export QD_USE_TOPO_ALBEDO=1
  python3 scripts/run_simulation.py
  ```
- 指定 orographic 增强（可选）：
  ```
  export QD_OROG=1
  export QD_OROG_K=7e-4
  python3 scripts/run_simulation.py
  ```
- 仍采用生成模式（不设置 QD_TOPO_NC）：
  ```
  python3 scripts/run_simulation.py
  ```

## 9. 后续扩展（与 P004/P006 的衔接）

- P004-L2：板块边界/侵蚀 → 更真实的山脉/海沟 → Orographic 参数化更具物理意义。
- 气候-地表反馈：利用 `elevation` 与 `land_mask` 引入地形/下垫面对边界层交换系数、粗糙度长度的影响。
- 可视化强化：地形阴影（hillshade）叠加、海岸线提取、多尺度剖面图、降水-地形相关性分析。

---

本项目聚焦“把 P004 生成的数据用起来”。完成后，GCM 将具有更合理的下垫面边界条件，为后续的地形降水与生态/文明模拟提供更可信的动力/热力背景。
