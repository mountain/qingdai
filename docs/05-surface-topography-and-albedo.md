# 5. 地形、海陆掩膜与反照率接入（Topography, Land–Sea Mask & Albedo）

本章整合项目 P004“行星地形生成”与 P005“地形接入 GCM”的实施方案，形成可复现的数据生产与模型接入规范。完成后，GCM 将以真实的地形、海陆掩膜与基础反照率作为边界条件运行，并可选启用地形降水增强。

状态（2025-09-21）
- P004 里程碑 1 已完成：L1+L3 地形、海平面自适应、导出标准化 NetCDF，并提供可视化。
- P005 集成已完成：从 NetCDF 读取 elevation/land_mask/base_albedo/friction；在动力与辐射路径中使用；提供 orographic precipitation 可选增强。

关联文档
- docs/04-runtime-config.md：完整环境变量目录（强烈建议结合本章阅读）
- projects/004-topography-generation.md：算法与参数细节
- projects/005-topography-integration-into-gcm.md：接入路径与运行示例


## 5.1 目标与范围

- 生产标准化的地形与地表属性数据：
  - elevation（m）、land_mask（0/1）、base_albedo（1）、friction（s⁻¹）
- 将上述数据接入 GCM：
  - 在动力核心使用空间变率摩擦（friction_map）
  - 在辐射与反照率合成中使用 base_albedo 底图，并与云/冰反照率耦合
- 提供可选的地形降水增强模块（orographic precipitation）
- 形成从“数据 → 接入 → 诊断”的一体化 SOP


## 5.2 数据格式（NetCDF）与字段

输出文件位于 data/ 目录，命名示例：
- data/topography_qingdai_181x360_seed42_YYYYMMDDTHHMMSSZ.nc

字段
- lat（degrees_north），lon（degrees_east）
- elevation（m）：海拔高度，已完成海平面自适应与柔性海岸平滑
- land_mask（0/1）：海陆掩膜，1=陆地，0=海洋（最近邻保持二值）
- base_albedo（1）：基础反照率底图（海洋低、陆地与高山略高、高纬可适当升高）
- friction（s⁻¹）：地表摩擦参数（陆地更大、海洋更小；山地可增强）

推荐/可选全局属性
- sea_level_m、target_land_fraction（例如 0.29）
- 生成参数快照（种子、权重、hypsometry 配置）以保证可追溯性


## 5.3 生成管线（P004 摘要）

目标：生成具“自然感”的全球地形、海陆比例≈0.29，并派生 base_albedo / friction。

混合方法（详见 projects/004-topography-generation.md）
1) L1 大尺度大陆轮廓：球面随机场（低频增强）+ Poisson 盘采样的“大陆种子”高斯丘陵叠加  
2) L2 板块边界（计划的里程碑 2）：发散/汇聚/转换边界带抬升与海沟（后续加入）  
3) L3 多倍频分形粗糙度 + 轻度侵蚀平滑  
4) Hypsometry 标定与海平面自适应，目标海陆比≈0.29±0.02

命令
- 生成：python3 -m scripts.generate_topography
- 可视化：python3 -m scripts.plot_topography
- 输出：data/ 目录下 .nc 与对应 *_overview.png

里程碑（节选）
- [x] L1+L3+海平面自适应 → elevation/land_mask/base_albedo/friction 导出
- [ ] L2 边界提升/海沟/侵蚀
- [ ] Hypsometry 标定扫描与参数固化
- [ ] 更丰富的可视化（阴影、剖面、海岸线度量）


## 5.4 模型接入（P005 摘要）

初始化阶段（示意）
1) 构建网格：grid = SphericalGrid(n_lat, n_lon)  
2) 读取外部地形：
   - 若设置了 QD_TOPO_NC，则从 NetCDF 读取并（必要时）插值到模型网格
   - 否则回退到内置生成/旧逻辑
3) 将 friction_map 注入动力学核心；将 base_albedo_map 注入反照率合成

核心 API（pygcm/topography.py）
- load_topography_from_netcdf(path, grid, *, regrid="auto")
  - 返回 (elevation, land_mask, base_albedo_map, friction_map)
  - 维度不一致时：经度按周期双线性插值；land_mask 用最近邻
- generate_base_properties(mask, elevation)（作为回退/派生）
  - 依据海陆/海拔生成 base_albedo/friction 的合理底图

反照率合成（pygcm/physics.py，与 P006/007/008 一致）
- calculate_dynamic_albedo(cloud_cover, T_s, alpha_base, alpha_ice, alpha_cloud)
  - alpha_base 可为标量（回退）或 2D 数组（base_albedo_map）
  - 与云/海冰共同决定总反照率 α_total


## 5.5 动态反照率：海冰与云的耦合

总反照率采用混合形式（示意）
- α_total = α_surface(T_s, type)·(1−C) + α_cloud·C
- α_surface 由 base_albedo_map（海/陆/高山差异）、海冰状态（α_ice）与温度等决定
- C 来自云量融合（降水驱动 + 背景源，见 docs/08-humidity-and-clouds.md 与 docs/03-climate-model.md）

环境变量（详见 docs/04-runtime-config.md）
- QD_USE_TOPO_ALBEDO（默认 1）：启用 base_albedo 底图
- QD_ALPHA_WATER、QD_ALPHA_ICE、QD_TRUECOLOR_*：参与可视化/辐射的一致渲染


## 5.6 地形降水增强（可选）

在“动力-降水”诊断之后，依据上坡风抬升增强迎风侧降水、抑制背风侧（简化方案）：
- 记近地层风 U=(u,v)，地形坡向单位向量 n̂ = ∇H/|∇H|（H=elevation）  
- uplift = max(0, U·n̂)  
- precip = precip_base · factor, 其中 factor = 1 + k_orog · uplift

实现与注意事项
- 梯度 ∇H 按球面几何换算实际距离（a·Δφ、a·cosφ·Δλ）
- uplift 与 factor 建议平滑；可设 factor 上限（如 ≤2.0）
- 环境变量：QD_OROG（0/1 开关），QD_OROG_K（默认 7e-4，建议微调）

启用命令（示例）
- export QD_OROG=1  
- export QD_OROG_K=7e-4


## 5.7 插值与一致性约束

- NetCDF→Grid 对齐：纬向常规、经向周期；高纬考虑 cosφ 缩放  
- land_mask：最近邻以保持二值；其它字段双线性  
- 一致性校验（建议日志打印）：
  - Land fraction（目标/达成）
  - friction/base_albedo 统计（min/mean/max）
  - elevation 统计（min/mean/max）  
- 若分辨率差异较大，建议在海岸线附近增加轻度平滑以避免“锯齿”


## 5.8 日志与诊断（建议）

模型启动与首张图前打印
- Topo source（NetCDF 路径或“Generated”）
- Land fraction（target/achieved）
- Elevation/base_albedo/friction 的统计  
- 若启用 orographic：打印 k_orog 与增强阈值/上限等

可视化增强（可选）
- 在 TrueColor 图上叠加海岸线/地形等值线
- 输出降水-地形的空间相关图作为验证


## 5.9 与其它模块的关系

- P006（能量收支）：base_albedo 进入短波/长波路径；夜侧温度下限 T_floor 防止极端冷塌
- P007（平板海洋/海冰）：α_ice 与 base_albedo 融合；海冰显著提升 α_total 并抑制蒸发
- P008（湿度/云）：降水驱动云量；云反照率/长波发射率改变辐射收支
- P010（反噪）：风场/降水条纹通过选择性耗散得到抑制；地形雨影更物理
- P011（动态海洋）：真实海岸线/大陆布局影响风应力与海洋环流格局
- P012（极点处理）：与地形无直接接口，但更一致的极区场有助于诊断与可视化


## 5.10 运行示例

使用 data 下最新外部地形（推荐）
```bash
export QD_TOPO_NC=$(ls -t data/*.nc | head -n1)
export QD_USE_TOPO_ALBEDO=1
python3 -m scripts.run_simulation
```

开启地形降水增强（可选）
```bash
export QD_TOPO_NC=$(ls -t data/*.nc | head -n1)
export QD_USE_TOPO_ALBEDO=1
export QD_OROG=1
export QD_OROG_K=7e-4
python3 -m scripts.run_simulation
```

生成与检查地形
```bash
python3 -m scripts.generate_topography
python3 -m scripts.plot_topography
open data/*_overview.png   # macOS 可视化检查
```


## 5.11 参数速查（与 04-runtime-config.md 一致）

- QD_TOPO_NC：外部 NetCDF 地形路径  
- QD_USE_TOPO_ALBEDO（默认 1）  
- QD_OROG（默认 0）、QD_OROG_K（默认 7e-4）  
- 反照率/可视化：QD_ALPHA_WATER、QD_ALPHA_ICE、QD_TRUECOLOR_*  
- 其余能量/湿度/海洋参数参见 docs/04-runtime-config.md


## 5.12 验收标准（建议）

- 功能性：设置 QD_TOPO_NC 后可正常运行并出图；未设置时回退生成/旧逻辑  
- 对齐一致性：插值无越界、海岸线连续、land_mask 二值稳定  
- 物理合理性：风/降水与地形关系合理；陆地上降水中心与山地/迎风相关；海上仍可存在带状 ITCZ  
- 输出与日志：首段日志包含 topo 来源与统计；图像能辨识海岸线/大陆构型、反照率对比


## 5.13 变更记录（Changelog）

- 2025‑09‑20：完成 P004 里程碑 1（L1+L3、自适应海平面、标准化导出）；完成 P005 接入、orographic 可选增强  
- 2025‑09‑21：文档迁移与整合至 docs/05‑surface‑topography‑and‑albedo.md；与 04‑runtime‑config.md 参数对齐，补充运行示例与验收标准
