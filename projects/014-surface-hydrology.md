# Project 014: 地表水文与径流路由模型（Surface Hydrology & Runoff Routing）

状态（2025‑09‑24）
- [x] 文档与方案定稿（本文件）
- [x] M1：离线流向网络与湖泊盆地预计算工具（脚本：`scripts/generate_hydrology_maps.py`）
- [x] M2：在线径流演算与流量累积模块（`pygcm/routing.py`，RiverRouting）
- [x] M3：GCM 主循环集成与新诊断输出（含湖泊；`scripts/run_simulation.py`）
- [x] M4：可视化增强（状态图与 TrueColor 叠加河流网络/湖泊图层）
- [ ] M5：参数标定与验证（径流时标、汇流速度等）

关联模块与文档
- P004（地形生成）、P005（地形接入）→ elevation/land_mask
- P009（行星水循环闭合）→ 陆地“桶”模型与径流量 R 的来源
- docs/04-runtime-config.md（运行参数总表；本项目新增环境变量将补充至此）


## 1. 目标与范围

动机
- 现有 P009 将陆地径流量 R 视作“即时回海”，未体现“从高处流向低处→河网→湖泊→下游”的汇流过程，难以刻画局地水资源分布（河流/湖泊/湿地）。
- 我们引入地表径流路由模型，以极低计算成本在行星尺度上模拟“江河湖海”的形成，支撑生态/文明学应用（如农业、定居点、交通等）的基础地理层。

目标（最小但自洽）
- 离线一次性预计算：从地形 elevation 构造无洼地的水文连通表面、流向网络、演算顺序与湖泊盆地（含溢出口）。
- 在线异步路由：按更长“水文步长”将 P009 产出的陆地径流 R 自上游汇入下游，形成流量累积（flow accumulation）；在湖泊单元上做简化的水量平衡（入流+局地降水−蒸发），可溢出至下游。
- 与能量/湿度保持类型一致性：湖面在辐射、蒸发上按“水体”处理（默认与海面相同），并在可视化与诊断中新增河网/湖泊层。
- 保守性：质量守恒在可控误差内，周期经度/球面几何一致，单位清晰、数值稳定。

不在本里程碑范围（可选未来扩展）
- 河道内蓄滞洪/动力波、河道糙率/糙率随流量变化、显式地下水/土壤水三维模型、河网向海的淡水通量耦合到海洋盐度（P011后续扩展）等。


## 2. 架构总览

双层策略：离线预计算 + 在线异步更新

- 离线（一次性）预计算（M1）
  1) 洼地填充（Pit Filling）：在球面经纬网格上将局地洼地填至最小溢出高程，得到无洼地的“水文连通地形” H_filled。
  2) 流向（D8 或 D∞ 简化）与下游指针：对每个陆地格点选择最陡坡方向（考虑球面距离），记录其唯一“下游邻格”与线性索引（`flow_to_index`）。
  3) 演算顺序（拓扑序）：自上游到下游排列陆地格点索引（`flow_order`），保证在线演算时“先算上游再算下游”不回环。
  4) 湖泊识别与溢出口：在 Pit Filling 结果中识别湖泊盆地，记录湖泊掩膜 `lake_mask`、湖泊唯一溢出口格点（`lake_outlet`）、湖面参考高程区间与最大面积近似（用于水量—水位—出流近似）。

  产物写入 `data/hydrology_network.nc`，详见 §3。

- 在线（每个水文步长）路由（M2/M3）
  - 输入：上一个水文步长累计的陆地径流通量 `R_land_flux`（来自 P009，单位 kg m⁻² s⁻¹）、可选的降水/蒸发通量（用于湖面局地平衡）。
  - 演算：按 `flow_order` 将各陆地格点量（乘以网格面积）汇入其唯一下游格点，并在湖泊单元执行“入流+局地 P−E → 湖量更新与溢流”。
  - 输出：河道“流量累积”场（kg s⁻¹ 或 m³ s⁻¹）、湖泊水量/水位变化、下游入海通量统计。
  - 频率：每 `QD_HYDRO_DT_HOURS`（默认 6 小时）触发一次；其余 GCM 步长仅累计，不演算。

- 集成点（主循环）
  1) P009 计算得到陆地桶径流 R（kg m⁻² s⁻¹，定义在陆地格点）。
  2) 将 R 与 dt 累加到路由模块缓存；当累计时长 ≥ `dt_hydro`：执行一次路由/湖泊更新。
  3) 诊断输出（可按 `QD_PLOT_EVERY_DAYS`）叠加河网/湖泊；质量守恒统计打印。


## 3. 数据产品规范：`data/hydrology_network.nc`

维度与坐标
- `lat(n_lat)`（degrees_north），`lon(n_lon)`（degrees_east）：与 GCM 网格一致（经度周期）。
- `n_land`（int）：陆地格点数量（供压缩存储的拓扑序列）。
- `n_lakes`（int）：湖泊数量（识别的互不连通盆地个数）。

变量（建议类型/单位）
- `land_mask(n_lat, n_lon)`：uint8，1=陆地，0=海洋（源自外部地形 NetCDF 或 P005 接入）
- `elevation_filled(n_lat, n_lon)`：float32，m；Pit Filling 后的连通地形（可用于诊断）
- `flow_to_index(n_lat, n_lon)`：int32，线性索引（行主或列主，需在文件属性中声明；海洋/出海口取 -1）
- `flow_dir(n_lat, n_lon)`：int8，D8 方向编码（1..8，见附录 A），用于人为可读/调试
- `flow_order(n_land)`：int32，线性索引数组（拓扑序，自上游至下游，仅包含陆地格点）
- `lake_mask(n_lat, n_lon)`：uint8，1=湖泊格点，0=非湖（湖泊覆盖于陆地之上）
- `lake_id(n_lat, n_lon)`：int32，每个湖泊唯一 id（1..n_lakes；非湖为 0）
- 表格式湖泊属性（长度 `n_lakes`）：
  - `lake_ids(n_lakes)`：int32，1..n_lakes
  - `lake_outlet_i(n_lakes)`、`lake_outlet_j(n_lakes)`：int16，溢出口格点行/列
  - `lake_h_min_m(n_lakes)`、`lake_h_max_m(n_lakes)`：float32，湖面高程近似范围（m）
  - `lake_Amax_m2(n_lakes)`：float32，最大静态湖面面积（m²）（用于简化的面积—水位曲线）

全局属性（建议）
- `indexing`：字符串，例如 `"row-major (i=j*nx + i), i=x(lon), j=y(lat)"`
- `projection`：`"latlon"`
- `created_by`、`created_at_utc`、`source_topography_nc` 等元数据
- `notes`：Pit Filling 方法、D8 细节、球面度量说明

说明
- `flow_to_index` 优先作为在线路由的唯一权威；`flow_dir` 便于可视检查。
- 若某陆地格点之下游为“出海口”，其 `flow_to_index=-1`；若为湖泊格点，下游仍可能是湖内流向或湖溢出口。


## 4. 离线预计算脚本（M1）

新建 `scripts/generate_hydrology_maps.py`（一次性运行）
- 输入：P004/P005 产出的地形 NetCDF（`elevation`、`land_mask`，可通过 `QD_TOPO_NC` 指定）
- 步骤：
  1) Pit Filling：可使用改造版 Priority-Flood 算法，球面经纬网格、经度周期、极圈邻接特殊处理。
  2) D8 流向：对每个陆地格点，比较 8 邻的“高程差/球面距离”，取最陡坡方向；如全部更高，则落入湖盆（由 Pit Filling 收敛避免闭环）。
  3) 拓扑序生成：对 `flow_to_index` 构成的 DAG 做一次拓扑排序（Kahn/BFS），得到 `flow_order`。
  4) 湖泊识别：依据 Pit Filling 的填平路径与最低溢出口识别湖盆，填充 `lake_mask`/`lake_id`，并记录溢出口格点；粗略估计 `lake_h_*` 与 `lake_Amax_m2`（可按湖盆像元数×像元面积上界）。
- 输出：写入 §3 规范的 `data/hydrology_network.nc`。
- CLI 示例：
  ```bash
  uv run python3 -m scripts.generate_hydrology_maps \
    --topo data/topography_qingdai_*.nc \
    --out data/hydrology_network.nc
  ```


## 5. 在线路由模块 API（M2）

新增 `pygcm/routing.py`
```python
class RiverRouting:
    def __init__(self, grid, network_nc_path, dt_hydro_hours=6.0,
                 treat_lake_as_water=True,  # 湖面按水体处理（辐射/蒸发）
                 alpha_lake=None,           # 可选: 覆盖湖面基础反照率（默认与海相同）
                 diag=True):
        """
        - grid: SphericalGrid，与 GCM 共网格
        - network_nc_path: data/hydrology_network.nc
        - dt_hydro_hours: 水文步长（默认 6 小时）
        - treat_lake_as_water: 湖面在 E/SH/辐射路径上作为“水体”类型
        """
    def reset(self):
        """清空累计缓存（accumulators）与诊断量。"""

    def step(self, R_land_flux, dt_seconds, precip_flux=None, evap_flux=None):
        """
        - R_land_flux: kg m^-2 s^-1，陆地格点的径流通量（来自 P009）
        - dt_seconds: 当前 GCM 步长（秒）
        - precip_flux: 可选 kg m^-2 s^-1，全域降水通量（用于湖泊局地平衡）
        - evap_flux: 可选 kg m^-2 s^-1，全域蒸发通量（用于湖泊局地平衡）
        逻辑：
          - 逐步累积到 internal buffer（单位 kg），当累计时间 >= dt_hydro 触发一次路由：
              1) 将陆地像元 buffer_mass = sum(R*dt*Area) 转为“入河质量”，沿 flow_order 自上游汇至下游
              2) 对 lake_mask==1 的像元：汇总入流 + （P−E）·dt·Area → 更新 lake_volume；
                 若 lake_volume 超阈（由 lake_h_max/Amax 近似），溢流质量按溢出口指针汇出
              3) 刷新“河道流量累积”（kg s^-1）与“入海质量速率”诊断量
          - 触发后清空 buffer，保留残余累计时间
        """
    def diagnostics(self):
        """
        返回 dict:
          - flow_accum_kgps (n_lat, n_lon): 河道流量累积
          - lake_volume_kg (n_lakes) 或 (n_lat, n_lon): 湖泊水量
          - ocean_inflow_kgps: 总入海质量速率
          - mass_closure_error_kg: 本水文步长质量闭合误差
        """
```

实现要点
- 单位转换：像元面积 A_cell=R²·Δλ·(sinφ_{j+1/2}−sinφ_{j−1/2})，与 docs/12 中面积权重一致。
- 邻接与周期性：经度周期 wrap（lon±1）相连，纬向考虑极圈邻接；与离线时一致。
- 触发频率：内部 `t_accum += dt`；当 `t_accum >= dt_hydro` → 执行路由并 `t_accum %= dt_hydro`。
- 性能：`flow_order` 是线性数组，路由一步 O(n_land)；湖泊运算 O(n_lakes)；开销很小。


## 6. 主循环集成（M3）

在 `scripts/run_simulation.py` 中
1) 初始化
```python
from pygcm.routing import RiverRouting

hydro_nc = os.getenv("QD_HYDRO_NETCDF", "data/hydrology_network.nc")
routing = RiverRouting(grid, hydro_nc, dt_hydro_hours=float(os.getenv("QD_HYDRO_DT_HOURS", "6")))
```

2) 每步调用（在 P009 完成 R 更新之后）
```python
routing.step(R_land_flux=R_flux_land,   # kg m^-2 s^-1 on land（由 update_land_bucket 返回的 runoff）
             dt_seconds=dt,
             precip_flux=P_flux,        # 若可得，传入湖面 P，用于湖库局地平衡
             evap_flux=E_flux)          # 若可得，传入湖面 E，用于湖库局地平衡
```

3) 诊断/出图（按 `QD_PLOT_EVERY_DAYS`）
```python
diag = routing.diagnostics()
# 叠加河网：将 flow_accum_kgps 映射为线宽/颜色
# 叠加湖泊：lake_mask → 填色/等值区
# 打印 [HydroRouting] ocean_inflow_kgps / mass_closure_error_kg
```


## 7. 与能量/湿度/水量路径的一致性

- 类型与通量
  - treat_lake_as_water=True 时，湖面在能量与湿度模块上等同于“水体”（海洋）：α/ε/E/SH 参数取海值或 `alpha_lake` 覆盖；蒸发通量会参与湖泊局地水量平衡（如传入 evap_flux）。
  - 路由仅对“陆地桶的径流 R_land”做汇流；海面 E−P 不进入河网，湖面 E/P 可作为“湖库局地源汇”。

- 保守性
  - 全局水量守恒扩展式：
    d/dt [⟨CWV⟩ + ⟨M_ice/ρ_w⟩ + ⟨W_land⟩ + ⟨S_snow⟩ + ⟨V_lake/ρ_w⟩]
      ≈ ⟨E⟩ − ⟨P⟩ − ⟨R_land⟩ − ⟨R_lake_out⟩ → 0（长期平均）
  - 在线模块每个水文步长输出“路由质量闭合误差”，用于回归。


## 8. 环境变量（新增；最终以 docs/04 为准）

- QD_HYDRO_NETCDF（默认 data/hydrology_network.nc）：离线网络文件路径
- QD_HYDRO_DT_HOURS（默认 6）：水文步长（小时）
- QD_TREAT_LAKE_AS_WATER（默认 1）：湖面作为水体处理（辐射与蒸发）
- QD_ALPHA_LAKE（可选；不设则与海洋相同）
- QD_HYDRO_DIAG（默认 1）：打印路由诊断（入海通量、闭合误差、最大流量等）


## 9. 数值与性能建议

- D8 vs D∞：首版采用 D8，简洁稳健；未来可选 D∞ 改进分配精度（上游质量按坡度比例分配至多个下游）。
- 球面度量：计算坡度时使用球面距离（a·Δφ、a·cosφ·Δλ）；极区 cosφ 下限保护，避免病态比值。
- 演算稳定性：拓扑序保证上游→下游一次通过；湖泊内部先合并入流，再处理溢流。
- 频率与成本：`dt_hydro` 取 1–24 小时均可；O(n_land+n_lakes) 的演算成本远低于大气/海洋步。


## 10. 诊断与可视化

- 日志（建议每次水文步长后）：
  - `[HydroRouting] ocean_inflow=... kg/s, max_flow=... kg/s, n_lakes=..., mass_error=... kg`
- 图层（与状态图叠加；M4 已集成）：
  - 河网（flow_accum_kgps）：在状态图（Ts 面板、Ocean 面板）以二值等值线叠加主干河流；
    - 阈值：`QD_RIVER_MIN_KGPS`（默认 1e6 kg/s）；透明度：`QD_RIVER_ALPHA`（默认 0.35）
  - 湖泊（lake_mask）：在上述面板以等值线勾勒（透明度 `QD_LAKE_ALPHA`，默认 0.40）
  - TrueColor 图：将主干河流/湖泊以半透明方式混合进 RGB 底图（河网默认 alpha=0.45）
  - 总开关：`QD_PLOT_RIVERS=1`（默认开启；设 0 关闭叠加）
- 时间序列（可选）：全球总湖量、最大/平均流量、入海通量


## 11. 验收标准（建议）

- 质量守恒：单次水文步长质量闭合误差 |error|/（输入质量总量） < 1e‑6（数值阈可调）
- 几何合理：流向从高到低（相对 H_filled）无反坡；河网连通、出海口合理
- 湖泊行为：长期平均湖量稳定或具季节周期；溢流在降水期增强
- 性能：分辨率 181×360 下，水文步长演算耗时 ≪ 动力学步（通常 < 1%）


## 12. 运行示例

1) 生成水文网络（离线）
```bash
uv run python3 -m scripts.generate_hydrology_maps \
  --topo $(ls -t data/topography_*.nc | head -n1) \
  --out data/hydrology_network.nc
```

2) 启用路由并运行 GCM（部分参数示例）
```bash
export QD_TOPO_NC=$(ls -t data/*.nc | head -n1)
export QD_USE_TOPO_ALBEDO=1

export QD_ENERGY_W=1
export QD_USE_OCEAN=1
export QD_USE_SEAICE=1

export QD_HYDRO_NETCDF=data/hydrology_network.nc
export QD_HYDRO_DT_HOURS=6
export QD_TREAT_LAKE_AS_WATER=1

python3 -m scripts.run_simulation
```

3) 可视化（叠加河网/湖泊）
- 在现有 TrueColor/状态图函数中读取 `routing.diagnostics()` 并叠加：
  - `flow_accum_kgps` → 线宽/颜色；`lake_mask` → 填色


## 13. 风险与边界情况

- 大型湖盆/内流区：若湖盆无真实出海口，允许作为内流终端（`lake_outlet=-1`）并在诊断中统计“内流区蓄水”；可设置长期轻微蒸发导致平衡。
- 地形分辨率差异：若外部地形分辨率与模型网格差异较大，建议在 P005 插值后再进行 Pit Filling，避免海岸线锯齿导致伪河道。
- 极地与经度收敛：高纬 Δλ 空间尺度收敛，注意坡度与距离归一；必要时对极圈采用更保守权重或湖泊硬约束（参照 docs/07 极区处理思路）。
- 单位一致性：R/E/P 在模块间统一使用 kg m⁻² s⁻¹；路由内部转为 kg（乘以面积）再回到 kg s⁻¹（除以 dt_hydro）。


## 14. 里程碑与任务拆解

- M1 离线工具
  - [ ] Pit Filling（球面+周期）、D8 流向
  - [ ] 拓扑序与湖泊识别/溢出口
  - [ ] 写出 `hydrology_network.nc`（§3 规范）
- M2 在线路由模块
  - [ ] `pygcm/routing.py`（RiverRouting 类）
  - [ ] 单元测试：质量守恒、简单地形回归（阶梯坡）、湖泊溢流逻辑
- M3 主循环集成
  - [ ] 在 `scripts/run_simulation.py` 中初始化与调用
  - [ ] 诊断打印与出图叠加
- M4 可视化与参数
  - [ ] 河网阈值/对数色标；湖泊渲染；交互式图层（可选）
  - [ ] `docs/04-runtime-config.md` 增补 P014 变量条目
- M5 标定与验收
  - [ ] dt_hydro 扫描（1/3/6/12/24 h）
  - [ ] 质量闭合与河网稳定性回归；产出验收报告图表


## 附录 A：D8 方向编码（建议约定）

- 以 (i=x=lon, j=y=lat) 为索引基底，j 向北递增，i 向东递增，经度周期
```
  7  8  1
  6  .  2
  5  4  3
```
- 对每个陆地格点，计算 8 邻的“高程差/球面距离”最陡者；若并列，按 1→8 顺序取首个。


## 附录 B：湖面面积—水位近似（简化）

- 为避免引入复杂湖盆 hypsometry，首版采用线性近似：
  - A(h) ≈ Amax · clip((h − h_min)/(h_max − h_min), 0, 1)
  - V(h) ≈ ∫ A(h) dh（显式公式可预存或按步近似）
- 在线更新：以质量（kg）等价为体积（m³）用 ρ_w≈1000 kg m⁻³ 转换。
- 当 V 超过对应 h_max 所致的容量上限时，视为“满盆溢流”，多余出流沿 `lake_outlet` 汇出。


## 附录 C：单位与常量

- 径流/降水/蒸发通量：kg m⁻² s⁻¹（与 P008/P009 一致）
- 河道流量累积：kg s⁻¹（或 m³ s⁻¹，保持 ρ_w 常数换算）
- 像元面积：m²（球面几何；与 docs/12 度量一致）
- 水密度：ρ_w=1000 kg m⁻³（可从环境变量 QD_RHO_W 读取）
