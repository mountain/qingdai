# Project 010: 更好的动力学——尺度选择性耗散与反噪（Hyperdiffusion）

本项目针对现有 GCM 输出中出现的“羽状云”（feather-like streaks）与对应的细丝状降水伪迹，给出标准、可控的数值稳定化方案：在动力核心中引入“尺度选择性耗散”（超扩散，∇⁴ 等），重点抑制网格尺度噪音，同时尽量不损伤大尺度环流结构与能量/水量闭合。

## 0. 状态（2025-09-20）

- [x] M1：问题诊断与方案选型（确认“羽状云”源自网格尺度噪音；选择 ∇⁴ 超扩散为主，Shapiro 滤波/谱掩膜为备选）
- [x] M2：实现 ∇⁴ 超扩散（u、v、h；可选 q、cloud 预留），显式稳定步进（σ₄ 自适应、极区稳定化、可选子步），参数由环境变量控制
- [ ] M3：参数扫描与默认值固化（σ₄、施加频率、场选择），并加入能量/水量不变性回归测试
- [x] M4：备选滤波器（Shapiro N=2/4；谱高波数带阻）与切换接口
- [ ] M5：验收与文档完善（图谱诊断、谱能量/散度高波数衰减报告）

关联实现基线：
- 当前 `pygcm/dynamics.py` 已有温和数值扩散（diffusion_factor≈0.995）、地转钳制与赤道正则化；未显式实现 ∇⁴ 超扩散或 Shapiro/谱滤波。

### 实现状态更新（2025-09-20）

- 代码落点（pygcm/dynamics.SpectralModel）：
  - 新增 `_laplacian_sphere` 与 `_hyperdiffuse`（两次 ∇² 复合得到 ∇⁴）
  - 极区稳定化：`cosφ` 下限 0.2；对输入场执行 NaN/Inf 清洗，避免度量发散
  - 在 `time_step` 中于摩擦后、云/降水诊断前施加超扩散（默认作用 `u,v,h`；`q/cloud` 可选）
  - 环境变量控制：
    - 开关与类型：`QD_DIFF_ENABLE`、`QD_FILTER_TYPE=hyper4`、`QD_DIFF_EVERY`
    - 强度与自适应：`QD_SIGMA4`（自适应 K₄=σ₄·Δx_min⁴/dt）、或显式 `QD_K4_U/V/H`（以及可选 `QD_K4_Q/CLOUD`）
    - 数值稳健：`QD_K4_NSUB`（子步数，默认 1）、`QD_DIFF_FACTOR`（全局温和扩散，默认 0.998）、`QD_DYN_DIAG`
    - Δx_min 取 `min(a·Δφ, a·cosφ·Δλ)`（最保守）
- 短程验证（dt=300 s，121×240）：
  - `QD_ENERGY_W=0` 与 `QD_ENERGY_W=1` 均稳定（12 步），未见数值爆裂；最大 |u|、|v|、|h| 在合理范围
- 建议默认参数：
  - `QD_SIGMA4=0.02`，`QD_K4_NSUB=1–2`，`QD_DIFF_EVERY=1`，`QD_DIFF_FACTOR=0.998`
  - 条纹严重时可增至 `QD_SIGMA4=0.03–0.04`；如需进一步稳健可设 `QD_K4_NSUB=2`
- 运行示例：
  ```
  export QD_ENERGY_W=1
  export QD_DIFF_ENABLE=1
  export QD_FILTER_TYPE=hyper4
  export QD_SIGMA4=0.02
  export QD_K4_NSUB=2
  export QD_DIFF_EVERY=1
  export QD_DIFF_FACTOR=0.998
  python3 -m scripts.run_simulation
  ```

## 1. 背景与问题

- 现象：云量图出现规则倾斜的高频条纹（“羽状云”），降水图与其一致。该结构与动力大尺度环流无关。
- 追因：降水基于散度/辐合诊断。风场中存在波长~2Δx 的网格尺度噪音（由平流/压梯/插值等数值误差激发），在求导时被放大，进入降水与云图，形成伪影。
- 结论：需在动力方程的数值离散中加入“强抑短波、弱抑长波”的耗散器，典型做法为高次拉普拉斯（∇⁴、∇⁶）或同等低通滤波（Shapiro、谱掩膜）。

参见原始分析草稿（保留于 `projects/10-better-dynamics.md`）。

## 2. 目标

- 抑制网格尺度与次网格尺度的噪音，显著减少云/降水中的条纹伪影。
- 保持行星波、急流等大尺度动力结构与能量/水量闭合性质。
- 以最小侵入方式集成：不重写核心框架，提供可选开关与参数，默认温和。
- 提供明确的稳定性约束与调参指导，含谱/图谱诊断与回归标准。

## 3. 设计与接口

### 3.1 核心方案：四阶超扩散（∇⁴）

- 动力方程附加项（以标量/矢量场 F 表示）：
  dF/dt = ... − K₄ ∇⁴F
  - ∇⁴F = ∇²(∇²F)
  - K₄（m⁴ s⁻¹）为超扩散系数，对最短波（~2Δx）抑制最强，对长波弱影响
- 作用对象（建议分场可控）：
  - 必选：u、v（风场）、h（高度/位势高度）
  - 可选：T_s（地表温度，通常不必）、q（湿度，温和）、cloud（云量，温和）
- 时间积分（显式，稳定性约束见 §5）：
  F ← F − K₄ ∇⁴F · dt
  - 实现上通过两次拉普拉斯：L = ∇²F，F ← F − K₄ ∇²L · dt

### 3.2 备选方案（可切换）

- Shapiro 滤波（N=2/4）：
  - 简易低通，每 N 步对选择字段执行一次，代价低、实现简单
- 谱掩膜（若谱变换管线就绪）：
  - 在 zonal FFT + 简化纬向投影后，对高波数带进行指数型阻尼
- 切换策略：
  - `QD_FILTER_TYPE=hyper4|shapiro|spectral`（默认 hyper4）
  - 允许组合：超扩散常开，Shapiro 低频率保底清噪（如每 6–12 步一次）

### 3.3 施加时机（time_step 内）

- 建议顺序（与现代码对齐的最小侵入改动）：
  1) 完成能量/湿度路径（计算 `u,v,h,T_s,q,cloud` 等更新）
  2) 进行 geostrophic 松弛、表面摩擦、半拉氏平流（现有）
  3) 在“云/降水等诊断进入下一步前”施加耗散：
     - Apply: ∇⁴ 超扩散于 u、v、h（必要时 q、cloud）
  4) 保留原有温和全局扩散（diffusion_factor≈0.995）或降低强度
- 目的：让降水/云的散度/辐合诊断“看见”已平滑的风场，切断“噪音→降水→云”的链条。

### 3.4 球面拉普拉斯（有限差分实现）

- 在半径 a 的球面上，纬度 φ、经度 λ 网格（均匀格距 Δφ、Δλ）：
  ∇²F = (1/a²) [ (1/ cosφ) ∂/∂φ (cosφ ∂F/∂φ) + (1/ cos²φ) ∂²F/∂λ² ]
- 数值实现要点：
  - 使用中心差分；经度方向周期边界；极区注意 cosφ→0 的数值稳健性（与现有 cos_lat_capped 一致）
  - 可重用 `self.dlat_rad`、`self.dlon_rad`、`self.grid.lat_mesh`、`self.a`
  - 计算 ∇²F 后再次对其套用 ∇² 以得 ∇⁴F

### 3.5 代码落点（建议）

在 `pygcm/dynamics.py` 中新增：

- 工具函数：
  - `_laplacian_sphere(F) -> ndarray`
  - `_hyperdiffuse(F, k4, dt, *, n_substeps=1) -> ndarray`
    - 若需稳定性裕度，可将 dt 拆分为 `n_substeps` 子步
- 主循环调用：
  - 在现有“温和扩散与 NaN 清理”段落前后，对选定字段执行：
    ```
    if diff_enabled and filter_type == "hyper4":
        u = _hyperdiffuse(u, k4_u, dt)
        v = _hyperdiffuse(v, k4_v, dt)
        h = _hyperdiffuse(h, k4_h, dt)
        if apply_k4_q:     q = _hyperdiffuse(q, k4_q, dt)
        if apply_k4_cloud: cloud_cover = _hyperdiffuse(cloud_cover, k4_c, dt)
    ```
- 诊断（可选）：
  - `QD_DYN_DIAG=1` 时每 ~200 步打印高波数方差比、散度谱峰值等（见 §6）

### 3.6 环境变量与默认值（建议）

- 开关与类型：
  - `QD_DIFF_ENABLE=1`（默认启用）
  - `QD_FILTER_TYPE=hyper4`（可选：`shapiro`、`spectral`）
  - `QD_DIFF_EVERY=1`（每步施加；或更大整数降低频率）
- 超扩散强度（两种设定方式，二选一）：
  1) 直接给系数（优先场分量可独立）：
     - `QD_K4_U=1.0e14`、`QD_K4_V=1.0e14`、`QD_K4_H=5.0e13`
     - 可选：`QD_K4_Q=2.0e13`、`QD_K4_CLOUD=1.0e13`（默认不启）
  2) 以无量纲 σ₄ 与网格/步长自适应（推荐）：
     - `QD_SIGMA4=0.02`（0.01–0.05 区间）
     - 计算：K₄ = σ₄ · (Δx_min)⁴ / dt，其中 Δx_min 取当前网格最小物理间距（高纬按 a·cosφ·Δλ）
- Shapiro 滤波（备选）：
  - `QD_SHAPIRO_N=2`、`QD_SHAPIRO_EVERY=6`
- 谱掩膜（备选）：
  - `QD_SPEC_CUTOFF=0.75`（相对 Nyquist 的截止），`QD_SPEC_DAMP=0.5`
- 诊断：
  - `QD_DYN_DIAG=1` 打印谱/方差指标

## 4. 参考实现（伪代码）

```python
def _laplacian_sphere(self, F):
    a = self.a
    dphi = self.dlat_rad
    dlmb = self.dlon_rad
    phi = np.deg2rad(self.grid.lat_mesh)
    cos = np.maximum(np.cos(phi), 1e-6)
    # 中心差分
    dF_dphi = (np.roll(F, -1, axis=0) - np.roll(F, 1, axis=0)) / (2*dphi)
    term_phi = (1.0/cos) * ( (np.roll(cos*dF_dphi, -1, axis=0) - np.roll(cos*dF_dphi, 1, axis=0)) / (2*dphi) )
    d2F_dlmb2 = (np.roll(F, -1, axis=1) - 2*F + np.roll(F, 1, axis=1)) / (dlmb**2)
    term_lmb = d2F_dlmb2 / (cos**2)
    return (term_phi + term_lmb) / (a**2)

def _hyperdiffuse(self, F, k4, dt, n_substeps=1):
    if k4 <= 0.0:
        return F
    sub_dt = dt / max(1, int(n_substeps))
    out = F
    for _ in range(int(n_substeps)):
        L = self._laplacian_sphere(out)
        L2 = self._laplacian_sphere(L)
        out = out - k4 * L2 * sub_dt
    return out
```

主循环调用（示意）：
```python
if int(os.getenv("QD_DIFF_ENABLE", "1")) == 1 and (step % int(os.getenv("QD_DIFF_EVERY", "1")) == 0):
    if os.getenv("QD_SIGMA4"):
        sigma4 = float(os.getenv("QD_SIGMA4"))
        dx_min = estimate_min_metric_length(self.grid, self.a)  # 取高纬 cosφ 缩放后的最小 Δx
        k4_base = sigma4 * dx_min**4 / dt
        k4_u = float(os.getenv("QD_K4_U", k4_base))
        k4_v = float(os.getenv("QD_K4_V", k4_base))
        k4_h = float(os.getenv("QD_K4_H", 0.5*k4_base))
    else:
        k4_u = float(os.getenv("QD_K4_U", "1.0e14"))
        k4_v = float(os.getenv("QD_K4_V", "1.0e14"))
        k4_h = float(os.getenv("QD_K4_H", "5.0e13"))
    self.u = self._hyperdiffuse(self.u, k4_u, dt)
    self.v = self._hyperdiffuse(self.v, k4_v, dt)
    self.h = self._hyperdiffuse(self.h, k4_h, dt)
    # 可选：q、cloud
```

## 5. 数值稳定性与调参指南

- 显式 ∇⁴ 稳定性（1D 标度）：约束为 (K₄·dt / Δx⁴) ≲ C_stab，经验常取 C_stab≈1/16–1/8。
- 在球面网格上，取最保守的 Δx_min（高纬方向受 cosφ 缩小），推荐以 `QD_SIGMA4=K₄·dt/Δx_min⁴` 表示，默认 0.02。
- 建议起步值：
  - u、v：σ₄≈0.02（羽状明显时 0.03–0.05）
  - h：σ₄≈0.01–0.02（避免过度平滑位势高度波）
  - q/云：仅在伪迹仍显著时施加，σ₄≤0.01
- 频率：每步（`QD_DIFF_EVERY=1`）；如需更弱耗散可设 2–4 步施加一次。
- 与现有扩散（0.995）关系：可保留但适当减弱其强度（例如 0.998），以把抑噪主责交给超扩散。

## 6. 诊断与验收标准

- 视觉验收：
  - [ ] 云/降水“羽状条纹”消失或显著减弱；无棋盘/格点印痕。
  - [ ] 大尺度环流（急流、行星波）保持结构与相位。
- 定量验收：
  - [ ] 高波数方差比例（上四分位数 k 带）较无超扩散时下降 ≥ 80%
  - [ ] 10–20 步时间窗内大尺度动能（截断于 k/k_Nyq<0.3）变化 |Δ| ≤ 5%
  - [ ] 全局能量诊断（P006）TOA/SFC/ATM 净差变化 |Δ| ≤ 2 W/m²（长期平均）
  - [ ] 水量闭合（P009）：⟨E⟩ 与 ⟨P⟩+⟨R⟩ 一致性变化 ≤ 3%
- 日志/谱诊断（可选实现 `QD_DYN_DIAG=1`）：
  - 打印高波数方差比、散度场的峰值/分位、谱能量在高波段的衰减系数

## 7. 任务拆解

- T1 实现算子
  - [ ] `_laplacian_sphere` 与 `_hyperdiffuse`，添加单元/数值测试（常数场→0、正弦波→特征值）
- T2 集成调用
  - [ ] 在 `SpectralModel.time_step` 中按 §3.3 施加于 u、v、h；加入环境变量读取
  - [ ] 调整现有 `diffusion_factor` 至较弱（如 0.998），避免重复过度耗散
- T3 参数扫描与默认组
  - [ ] 扫描 `σ₄`、`QD_DIFF_EVERY`，按 §6 指标选取默认
- T4 备选滤波器
  - [ ] 实现 Shapiro(N, every) 与谱掩膜（可只对 zonal FFT 维度）
- T5 诊断与文档
  - [ ] 输出高波数方差比/散度谱；更新 README 与本项目文档示例

## 8. 运行示例

- 基本：启用超扩散（自适应 σ₄）
```bash
export QD_DIFF_ENABLE=1
export QD_FILTER_TYPE=hyper4
export QD_SIGMA4=0.02
export QD_DIFF_EVERY=1
# 可选：弱化原先全局扩散（若暴露为环境变量，或在代码固定）
# export QD_ENERGY_W=1 推荐与能量收支框架共同运行
python3 -m scripts.run_simulation
```

- 强化抑噪（条纹严重场景）：
```bash
export QD_SIGMA4=0.04
export QD_DIFF_EVERY=1
python3 -m scripts.run_simulation
```

- 显式指定 K₄（避免自适应）：
```bash
export QD_DIFF_ENABLE=1
unset QD_SIGMA4
export QD_K4_U=1.0e14
export QD_K4_V=1.0e14
export QD_K4_H=5.0e13
python3 -m scripts.run_simulation
```

- 叠加 Shapiro（每 6 步一次）：
```bash
export QD_FILTER_TYPE=hyper4
export QD_SHAPIRO_N=2
export QD_SHAPIRO_EVERY=6
python3 -m scripts.run_simulation
```

## 9. 与其他项目的交叉引用

- P005（地形接入）：风场/降水与地形的相干性更清晰，抑制噪音带来的“假”地形雨影。
- P006（能量收支）：超扩散对能量收支仅以“数值耗散”形式出现；保持诊断与守恒核算（TOA/SFC/ATM）。
- P007（平板海洋/海冰）：不改变海冰热力学；平滑后的风场能提高 SH/LH 通量稳定性。
- P008（湿度 q）：可选对 q 施加较弱超扩散；减少 E−P 噪音耦合的回响。
- P009（水循环）：闭合诊断更稳定；E/P/R 时序更物理，极端像元抖动减少。

## 10. 实施备注与建议

- 先只对 u、v、h 施加 ∇⁴（推荐），观察“羽状云”是否消除。若仍有残留，再温和作用于 q/云。
- 高纬最严格：计算 K₄ 时使用 Δx_min（含 cosφ），必要时对子步（`n_substeps=2`）以增加稳定裕度。
- 默认不关闭现有温和扩散，但建议减弱，让“选择性耗散”成为主器。
- 验收以“高波数方差比下降 + 大尺度动能保持 + 能量/水量诊断稳健”为准。

## 11. 进一步改进建议（可选扩展）

- 非线性反混叠 / 谱去混叠（Dealiasing，2/3 规则）
  - 在经向 FFT 维度对高于 m_cut ≈ 2/3·m_Nyq 的波数做零化或指数阻尼，降低非线性混叠伪能量。
  - 建议配置：`QD_SPEC_DEALIAS=1`、`QD_SPEC_DEALIAS_RULE=two_thirds`。

- 发散选择性阻尼（div-damping）
  - 对散度分量施加更强的耗散以优先抑制“压缩型”噪音；对涡度分量使用较弱耗散以保留旋转流结构。
  - 网格近似路径：先诊断 `D=∇·V`、`Z=∇×V`，对 D 用较大 K₄；简化实现可对 `u,v` 直接添加与 `∇D` 成正比的校正项。
  - 参数建议：`QD_K4_DIV`、`QD_K4_VORT`（若未提供则回退至统一 `QD_K4_U/V`）。

- Smagorinsky 变系数涡粘（流场自适应）
  - 让耗散强度随形变率 |S| 自适应，K ∝ (C_smag·Δ)²·|S|；在剪切/锋区增强抑噪，平坦区弱化。
  - 参数：`QD_SMAG_C`（如 0.12–0.20）、`QD_SMAG_ORDER=2|4`（选择 K₂ 或 K₄ 形态）。

- 标量保形平流（monotone advection）
  - 为 T_s、q、cloud 引入通量限制器或 WENO5 等高阶、单调性保持的平流方案，降低插值振铃。
  - 切换接口（占位）：`QD_ADV_SCHEME=bilinear|monotone|weno5`，`QD_ADV_LIMITER=vanleer|minmod|superbee`。

- 准单调后处理（quasi-monotone filter）
  - 半拉氏插值后，将像元值裁剪到局部极值包络内，避免新生超/欠冲导致的高频伪迹。
  - 开关：`QD_QUASI_MONO=1`。

- 时间滤波器（RAW / Asselin）
  - 若未来改用 leapfrog，可引入 RAW 过滤器消除计算模；当前也可对 `u,v` 施加极小的时间平滑：
  - 参数：`QD_TIME_FILTER_EPS`（默认 0；建议 ≤ 0.01）。

- 自适应步长 / CFL 保护
  - 在运行脚本中按 `max(|u|,|v|)` 动态调整 `dt`，满足 `dt ≤ CFL·Δx_min/maxV`。
  - 参数：`QD_ADAPT_DT=1`、`QD_CFL=0.5`（默认建议 0.4–0.7）。

- 极区海绵层 / 极区强化耗散
  - 在 |φ| ≥ φ₀ 区域增大 K₄ 或提高 Shapiro 频率，缓解 cosφ 缩放导致的高纬最严格稳定性。
  - 参数：`QD_POLAR_SPONGE_LAT=70`、`QD_POLAR_SPONGE_GAIN=3`。

- 诊断增强（谱 / 能量 / 涡度-散度）
  - 输出 KE/Enstrophy 谱、Div/Vort 方差比、最大 CFL 数、谱能量高波段衰减因子。
  - 扩展诊断开关：`QD_DYN_DIAG_SPECTRA=1`（基于 `QD_DYN_DIAG=1`）。

- 网格与离散改进（路线图）
  - 远期可过渡到 Arakawa C-grid 或采用矢量不变形式与 Arakawa Jacobian，提高能量/涡量守恒；不在本里程碑实现。

- 临时缓解（仅限开发期）
  - 对降水/云图进行小窗口中值滤波或形态学开运算作为临时去噪（不建议长期启用，以免掩盖数值问题）。
  - 参数：`QD_PRC_MEDIAN=0|1`、`QD_PRC_WIN=3|5`。

---

附：原始草稿（projects/10-better-dynamics.md）指出问题源于“网格尺度噪音→散度放大→降水/云伪影”，本项目即以标准数值技术（超扩散/滤波）解决这一经典问题，使模型从“能跑”迈向“跑得好”。

---

### M4 备选滤波器实现与默认（2025-09-20）

已在 `pygcm/dynamics.py` 中实现并接入两类滤波器，支持与 ∇⁴ 超扩散组合（combo）使用：

- Shapiro 滤波（N=2/4，可配置）  
  - 实现：`_shapiro_filter`，通过分离的 1–2–1 核（经向 wrap、纬向 nearest）重复 N 次。  
  - 用途：低开销去除网格尺度锯齿，配合超扩散进一步抑制“羽状云”。

- 谱带阻（Zonal FFT 高波数阻尼）  
  - 实现：`_spectral_zonal_filter`，对相对 Nyquist 高于 `QD_SPEC_CUTOFF` 的波段按 `QD_SPEC_DAMP` 衰减。  
  - 用途：可选；在极端条纹场景与 combo 搭配，每若干步施加一次。

默认与推荐
- 默认 `QD_FILTER_TYPE=combo`（超扩散 + Shapiro）。  
- 默认参数：  
  - `QD_SIGMA4=0.02`（自适应 K₄=σ₄·Δx_min⁴/dt），`QD_DIFF_EVERY=1`，`QD_K4_NSUB=1`  
  - `QD_SHAPIRO_N=2`，`QD_SHAPIRO_EVERY=6`  
  - 谱阻尼默认关闭（需显式 `QD_SPEC_EVERY>0` 才启用）  
  - 温和全局扩散：`QD_DIFF_FACTOR=0.998`
- 若条纹仍然可见：  
  - 将 `QD_SIGMA4` 调至 0.03–0.04，或 `QD_K4_NSUB=2` 增稳定裕度；  
  - 保持 Shapiro N=2，每 6 步；必要时开启谱阻尼（如每 6–12 步，`QD_SPEC_CUTOFF≈0.7`，`QD_SPEC_DAMP≈0.5`）。

运行示例
```bash
# 默认温和（combo: hyper4 + Shapiro），已显著减弱羽状
export QD_FILTER_TYPE=combo
export QD_SIGMA4=0.02
export QD_SHAPIRO_N=2
export QD_SHAPIRO_EVERY=6
export QD_DIFF_EVERY=1
export QD_DIFF_FACTOR=0.998
python3 -m scripts.run_simulation
```

```bash
# 条纹较强时的增强方案（仍保持数值稳定）
export QD_FILTER_TYPE=combo
export QD_SIGMA4=0.03
export QD_K4_NSUB=2
export QD_SHAPIRO_N=2
export QD_SHAPIRO_EVERY=6
# 可选启用谱阻尼（低频率）
export QD_SPEC_EVERY=6
export QD_SPEC_CUTOFF=0.70
export QD_SPEC_DAMP=0.5
python3 -m scripts.run_simulation
```

接口（新增/更新）
- `QD_FILTER_TYPE=hyper4|shapiro|spectral|combo`（默认 combo）  
- Shapiro：`QD_SHAPIRO_N`（默认 2）、`QD_SHAPIRO_EVERY`（默认 6）  
- 谱阻尼：`QD_SPEC_EVERY`（默认 0=关闭）、`QD_SPEC_CUTOFF`（默认 0.75）、`QD_SPEC_DAMP`（默认 0.5）  
- 其余见 §3.6（`QD_SIGMA4`、`QD_K4_*`、`QD_DIFF_*` 等）
