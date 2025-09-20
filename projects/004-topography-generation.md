# Project 004: 行星地形与海陆分布生成（Topography & Land-Sea Mask）

本文档提出一套可复现实用的“自然感”地形生成方案，用于替代当前矩形大陆，支持 2–3 大陆构型、接近地球的海陆比例（陆地约 29%），并提供海拔高度场（elevation），使其可派生出地表摩擦、基础反照率、地形降水等二次参数。

## 状态（2025-09-20）
- [x] 里程碑 1：L1 + L3 + 海平面自适应，生成并导出多字段 NetCDF（elevation、land_mask、base_albedo、friction）。
  - 生成：`python3 -m scripts.generate_topography`
  - 可视化：`python3 -m scripts.plot_topography`
  - 输出目录：`data/`（示例：`data/topography_qingdai_181x360_seed42_20250920T114018Z.nc` 与对应 `_overview.png`）
- [ ] 里程碑 2：引入 L2 板块边界提升/俯冲与侵蚀
- [ ] 里程碑 3：Hypsometry 标定与参数扫描
- [ ] 里程碑 4：集成到 GCM：导出 elevation/land_mask/friction/albedo，联动地形降水
- [ ] 里程碑 5：可视化工具完善（阴影、剖面、海岸线度量等）

## 1. 目标与约束

- 生成全球海拔高度场 `H(lat, lon)`（单位 m），并给出海平面高度 `H_sea`，以此得到二值海陆掩膜 `mask = 1[H >= H_sea]`。
- 大陆数量：可配置（目标：2–3 个大洲），形态自然，拥有大尺度轮廓与多尺度细节。
- 海陆面积比例：目标 ≈ 地球（陆地面积 ≈ 29%），允许 ±1–2% 浮动，通过自适应阈值收敛。
- 高度分布（类地 hypsometry）：海洋平均深度约 -3.7 km，陆地平均高度约 +0.75 km，分布呈双峰（海洋盆地与大陆块）。
- 可重复性：随机种子控制，参数化，可快速再现。
- 兼容 GCM：在 `pygcm/topography.py` 中可导出 `elevation`, `land_mask`, `base_albedo_map`, `friction_map` 等。

## 2. 总体思路（混合管线）

采用“低频大陆控制 + 板块风格边界 + 分形细节 + 侵蚀平滑 + 海平面自适应”的混合方法。核心包含三层频率：

1) 大尺度（L1）：控制大陆块数量与轮廓  
2) 中尺度（L2）：板块构造启发的造山、洋中脊、俯冲带  
3) 小尺度（L3）：分形粗糙度（多倍频噪声）与侵蚀平滑

最终高度为
H = w1·H_L1 + w2·H_L2 + w3·H_L3
经 hypsometry 映射与海平面阈值调整后输出。

## 3. 大尺度大陆轮廓（L1）

提供两种互补方案，可单独或叠加使用：

- 方案 A：球面随机场（谱合成）  
  - 在球面上生成高斯随机场，功率谱密度 P(k) ∝ k^(-β)，低频增强（β≈2–3），保留大尺度起伏。  
  - 可使用球谐随机合成（截断到较低波数 L≈8–12）或在经纬度网格上用带限噪声（对频谱做 Butterworth 低通）。

- 方案 B：Poisson 盘采样种子 + 高斯基函数叠加  
  - 在球面上用 Poisson-disk 采样 N_continents（2–3）个种子点 {c_i}，两点最小间距控制大陆分离度。  
  - 每个种子生成一个广义高斯丘陵：  
    H_L1_i(θ) = A_i · exp[-(d_gc(θ, c_i)/σ_i)^p]  
    其中 d_gc 为大圆距离，A_i>0（大陆），可少量放入海洋种子 A_i<0（形成大洋盆地）。  
  - 总和并归一化，得到 L1 轮廓；A_i、σ_i、p、N_continents 控制大陆数量与尺度。

建议：A+B 融合。A 提供自然的低频纹理，B 提供明确的大洲数与分离度控制。

## 4. 板块风格边界（L2）

- 随机生成 N_plate plate seeds（可与大陆种子相关或独立，典型 8–16）。  
- 对每个网格点分配最近的 plate seed → 球面 Voronoi（板块划分）。  
- 标注边界类型（概率或规则分配）：  
  - 发散（洋中脊）：边界附近给正向抬升带，幅值随距离衰减（带宽 300–600 km）。  
  - 汇聚（造山带/俯冲）：给非对称抬升与海沟负高，造山峰值靠陆侧，海沟在洋侧（双带结构）。  
  - 转换：轻微扰动即可。  
- 定义 L2 场：按边界类型计算提升/下沉，并在板块内部给少量缓起伏（模拟克拉通、中大陆高原）。

建议：距离函数 d_to_boundary 与法向侧（大陆/洋）联合控制：  
H_L2 = Σ_boundaries [ s_type · A_type · exp(-d^2/(2σ_b^2)) · side_factor ]

## 5. 分形细节与侵蚀（L3）

- 分形噪声：多倍频噪声（Perlin/Simplex/fBM）叠加 O(4–6) 个倍频，频率按 2^octave 递增，幅度按 Hurst 指数 H≈0.7–0.9 衰减。  
- 地形导向增益：在高海拔处增强粗糙度，在低海拔处稍弱，以模拟山地粗糙与海洋相对平缓。  
- 侵蚀/平滑：应用若干次迭代的“热力侵蚀（thermal erosion）”或各向异性扩散，限制坡度，去除伪锯齿。轻度即可（2–5 迭代）。

## 6. Hypsometry 与海平面自适应

- 目标分布：  
  - 海洋：均值 -3700 m（σ≈1200 m，长尾），  
  - 大陆：均值 +750 m（σ≈500 m，右偏），  
  - 双峰总体 hypsometry。  
- 方法：  
  1) 对 H 做单调映射 f(H) → 目标分布（可用分位数映射/直方图匹配，使分布近似目标曲线）。  
  2) 海平面 `H_sea` 通过二分搜索调整，使 `area(H >= H_sea) ≈ 0.29`。  
  3) 对海平面附近的高程施加轻柔 S 形过渡（例如 logistic），减少“海岸线台阶”。

## 7. 参数与默认值（建议）

- N_continents = 3（可 2–3）  
- Poisson 种子最小间距 ≈ 35–45°  
- L1：β=2.5，球谐截断 L≈10；或低通滤波截止波数对应 ~4000–6000 km 波长  
- L2：边界带宽 σ_b = 400 km；发散 A_div ≈ +1200 m（洋脊）；汇聚山脉 A_con ≈ +2500 m；海沟 A_trench ≈ -3000 m  
- L3：fBM octaves=5，Hurst H≈0.8，振幅总量 ≈ 600 m（山地增强、平原减弱）  
- 侵蚀迭代 3 次，热侵蚀角阈 35–40°  
- Hypsometry：目标均值（海洋 -3700 m，陆地 +750 m），直方图匹配强度 0.5–0.7（与原 H 保持一定相似度）  
- 最终海陆比目标 0.29，容差 0.01–0.02

## 8. 伪代码

```
Input: grid (lat, lon), seed

# L1: 大洲轮廓
S_low = band_limited_noise(grid, seed, beta=2.5, cutoff=L=10)
C = poisson_disk_samples_on_sphere(N_continents=3, min_dist=40°, seed)
H_l1 = normalize(S_low)
for each c in C:
    H_l1 += A_c * exp(-(dist_gc(grid, c)/sigma_c)^p)
H_l1 = normalize(H_l1)

# L2: 板块边界
P = poisson_disk_samples_on_sphere(N_plate=12, min_dist=20°, seed+1)
labels = spherical_voronoi_labels(grid, P)
B = extract_voronoi_boundaries(labels)
H_l2 = 0
for each boundary b in B:
    type = sample_type({div, conv, trans}, probs)
    d = distance_to_boundary(grid, b)
    side = continental_side_indicator(grid, b, labels)  # or random but consistent
    H_l2 += uplift_profile(type, d, side)

# L3: 分形与侵蚀
H_l3 = fractal_fbm(grid, seed+2, octaves=5, H=0.8)
H_l3 *= amplitude_modulate_by(H_l1 + H_l2)

# 融合与侵蚀
H_raw = w1*H_l1 + w2*H_l2 + w3*H_l3
H_eroded = thermal_erosion(H_raw, iterations=3, angle_thresh=38°)

# Hypsometry 映射与海平面调整
H_hist = histogram_match(H_eroded, target_hypsometry_curve)
H_sea = search_threshold(H_hist, target_land_fraction=0.29, tol=0.01)
land_mask = (H_hist >= H_sea)

# 柔性海岸线
H_final = coast_smooth(H_hist, H_sea, width=100 m)
Output: elevation=H_final, land_mask
```

## 9. 与 GCM 的集成

- 在 `pygcm/topography.py` 中新增：
  - `generate_elevation_map(grid, seed, params) -> elevation`
  - `create_land_sea_mask_from_elevation(elevation, target_land_frac=0.29) -> mask, H_sea`
  - `generate_base_properties(mask, elevation) -> base_albedo_map, friction_map`
    - `friction_map` 可随海拔增强（山地阻力更大），海洋更弱。
    - `base_albedo_map` 可在高纬／高海拔略增，低纬海洋维持低值。
- 未来扩展：或ographic precipitation（地形降水）参数化可直接使用 `elevation` 与坡向、迎风/背风指标。

## 10. 评估与校准

- 面积约束：确保 `area(mask==1)/area_total ≈ 0.29±0.02`。  
- 形态指标：大陆紧凑度、海岸线分形维数、盆地/山脉尺度分布。  
- 高程分布：对照目标 hypsometry 的 KL 散度或分位偏差。  
- 视觉检查：大陆数量、分布合理性、是否出现“格点伪迹”或“棋盘纹”。

## 11. 里程碑

1) 原型实现（L1+L3 + 海平面自适应），达到 2–3 大陆、地球海陆比。  
2) 引入 L2 板块边界提升/俯冲与侵蚀，形成山脉/海沟/洋脊特征。  
3) Hypsometry 标定与参数扫描，稳定再现目标分布。  
4) 集成到 GCM：导出 elevation/land_mask/friction/albedo，联动地形降水。  
5) 可视化工具：海拔、分布直方图、海岸线、山脉剖面、地形阴影渲染。

## 12. 参数清单（可通过配置/环境变量）

- `N_CONTINENTS`（默认 3）、`SEED`  
- `BETA`, `SPECTRAL_L_MAX`（L1 球谐/低通滤波）  
- `PLATE_COUNT`, `BOUNDARY_WIDTH_KM`, `A_div`, `A_con`, `A_trench`  
- `FBM_OCTAVES`, `HURST_H`, `FBM_AMP`（L3）  
- `EROSION_ITERS`, `EROSION_ANGLE_DEG`  
- `TARGET_LAND_FRAC`（默认 0.29）、`HYPSO_MATCH_STRENGTH`  
- `COAST_SMOOTH_WIDTH_M`

## 13. 参考实现建议

- 噪声：`numpy` + 简单 FBM；高级可用 `noise` 包（Perlin/Simplex）。  
- 球面几何：自行实现大圆距离或用 `pyproj` 等库（可选）。  
- 球谐：`pyshtools`（进阶）；原型阶段可在经纬网格用频域低通。  
- 侵蚀：实现简化版热力侵蚀（阈值坡度、物质转移），或各向异性扩散近似。

---

附注：上述方案能直接控制“大陆数量与海陆比例”，同时通过板块启发与分形细节实现自然外观。Hypsometry 标定保证海洋/陆地高度分布更“类地”，为后续地形-气候相互作用（如地形降水）打下基础。
