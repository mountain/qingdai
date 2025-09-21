# 10. 数值稳定与反噪（Numerics & Stability: Hyperdiffusion, Filters, CFL）

本章迁移项目 P010“更好的动力学/反噪”的方案，规范 Qingdai GCM 的数值稳定化策略与调参方法，重点抑制网格尺度噪音（如“羽状云”条纹），同时尽量保留大尺度环流结构与守恒诊断稳定。

状态（2025-09-21）
- 已实现：四阶超扩散（∇⁴）用于 u、v、h（可选 q、cloud），自适应强度 σ₄；Shapiro 滤波与谱带阻作为备选/组合；默认“combo”已接入
- 待完善：参数扫描默认组固化、高波数方差/谱能量诊断面板

关联文档
- docs/06-energy-framework.md：能量诊断（TOA/SFC/ATM）
- docs/08-humidity-and-clouds.md：云–降水条纹与湿度一致性
- docs/07-ocean-and-sea-ice.md：海洋模块中的反噪参数（ocean_*）
- docs/04-runtime-config.md：运行时环境变量总表（含 P010 条目）
- projects/010-better-dynamics.md：详细设计与路线图


## 10.1 问题与目标

- 现象：云/降水图出现规则倾斜“羽状”条纹，源于风场的网格尺度噪声在求散度时被放大并传导到降水与云量。
- 目标：用“尺度选择性耗散”与低通滤波抑制短波噪音，保留行星波、急流等大尺度结构；保持能量/水量闭合诊断稳健。


## 10.2 方案概览

核心：四阶超扩散（∇⁴）
- 方程附加项：dF/dt = … − K₄ ∇⁴F（对 F∈{u, v, h}；可选 q、cloud）
- 性质：对最短波（~2Δx）抑制最强，对长波影响弱

备选/组合滤波
- Shapiro（N=2/4）低频率应用，低开销去除锯齿
- 谱带阻（zonal FFT）对相对 Nyquist 高波数带做阻尼
- 默认“combo”：超扩散 + Shapiro；谱带阻默认关闭

施加时机（time_step 内）
- 在地表摩擦、半拉氏平流等主更新后、云/降水诊断之前施加，切断“噪音→降水→云”的放大链条


## 10.3 数学与实现要点（简要）

球面拉普拉斯（有限差分）
- ∇²F = (1/a²)[ (1/cosφ)∂φ(cosφ ∂φF) + (1/cos²φ)∂²λF ]
- 数值：经度周期边界、极区 cosφ 下限保护、中心差分

超扩散算子
- 通过两次拉普拉斯获得 ∇⁴F：L = ∇²F；L2 = ∇²L
- 显式时间推进：F ← F − K₄ L2 · dt；支持子步（n_substeps）提升稳定裕度


## 10.4 环境变量与默认值（详见 docs/04-runtime-config.md）

主控
- QD_DIFF_ENABLE（默认 1）：启用数值抑噪
- QD_FILTER_TYPE（默认 combo，可选 hyper4|shapiro|spectral|combo）
- QD_DIFF_EVERY（默认 1）：施加频率（步）
- QD_DIFF_FACTOR（默认 0.998）：温和全局扩散（乘法因子）

超扩散（∇⁴）强度（两种设定）
- 自适应无量纲 σ₄（推荐）：QD_SIGMA4（默认 0.02）  
  K₄ = σ₄·(Δx_min)⁴/dt（Δx_min 包含高纬 cosφ 缩小）
- 显式 K₄（m⁴ s⁻¹）：QD_K4_U/V/H（字段可分开）；可选 QD_K4_Q/CLOUD

子步与诊断
- QD_K4_NSUB（默认 1）：超扩散子步
- QD_DYN_DIAG（默认 0）：动力学反噪诊断打印（高波数方差、散度峰值等，若实现）

Shapiro 滤波
- QD_SHAPIRO_N（默认 2）、QD_SHAPIRO_EVERY（默认 6）

谱带阻（默认关闭）
- QD_SPEC_EVERY（默认 0 关闭）、QD_SPEC_CUTOFF（默认 0.75）、QD_SPEC_DAMP（默认 0.5）


## 10.5 调参与稳定性建议

经验起点
- σ₄≈0.02（u、v、h），条纹严重可至 0.03–0.04；h 可稍弱（0.01–0.02）
- 施加频率每步（QD_DIFF_EVERY=1）；如需更弱可 2–4 步一次
- 保留温和全局扩散（0.998）但避免过强，主责交给 ∇⁴

显式稳定性约束（标度）
- (K₄·dt / Δx⁴) ≲ C_stab；经验 C_stab≈1/16–1/8；自适应 σ₄ 表示法可自动满足

何时对 q/云施加 ∇⁴
- 仅在“combo”后条纹仍残留时，用小 σ₄（≤0.01）温和作用于 q/云；优先处理风/位势高度

极区与高纬
- Δx_min 取包含 cosφ 的最严格尺度；可用子步（QD_K4_NSUB=2）提升稳定裕度
- 可与海洋极区 sponge/极点一致化（docs/07）配合


## 10.6 诊断与验收标准（建议）

视觉
- 云/降水“羽状”条纹消失或显著减弱，无棋盘印痕
- 大尺度环流（急流、行星波）结构/相位基本保持

定量
- 高波数方差比例（上四分位波段）较无超扩散时下降 ≥ 80%
- 10–20 步窗内大尺度动能（k/k_Nyq<0.3）变化 |Δ| ≤ 5%
- 能量闭合（docs/06）：TOA/SFC/ATM 长期平均净差变化 |Δ| ≤ 2 W/m²
- 水量闭合（docs/09）：⟨E⟩ 与 ⟨P⟩+⟨R⟩ 一致性变化 ≤ 3%


## 10.7 运行示例

默认温和（combo）
```bash
export QD_FILTER_TYPE=combo
export QD_SIGMA4=0.02
export QD_SHAPIRO_N=2
export QD_SHAPIRO_EVERY=6
export QD_DIFF_EVERY=1
export QD_DIFF_FACTOR=0.998
python3 -m scripts.run_simulation
```

增强抑噪（条纹显著）
```bash
export QD_FILTER_TYPE=combo
export QD_SIGMA4=0.03
export QD_K4_NSUB=2
export QD_SHAPIRO_N=2
export QD_SHAPIRO_EVERY=6
# 可选启用谱带阻（低频率）
export QD_SPEC_EVERY=6
export QD_SPEC_CUTOFF=0.70
export QD_SPEC_DAMP=0.5
python3 -m scripts.run_simulation
```

显式 K₄（不使用 σ₄）
```bash
export QD_DIFF_ENABLE=1
unset QD_SIGMA4
export QD_K4_U=1.0e14
export QD_K4_V=1.0e14
export QD_K4_H=5.0e13
python3 -m scripts.run_simulation
```


## 10.8 与其它模块的关系

- docs/06：能量诊断应随反噪参数变化保持长期近守恒
- docs/08：降水/云条纹由 P010 明显缓解；湿度—云一致性更易检验
- docs/07：海洋同样使用 ∇⁴/Shapiro（ocean_* 参数），并有极点一致化与 sponge
- docs/05：地形雨影在抑噪后更物理，避免噪音主导的“假”地形响应


## 10.9 变更记录（Changelog）

- 2025‑09‑20：实现 ∇⁴、自适应 σ₄、Shapiro/谱带阻接口；默认“combo”
- 2025‑09‑21：参数建议与诊断标准补充；与能量/湿度/海洋交叉引用
- 2025‑09‑21：文档迁移与整合至 docs/10‑numerics‑and‑stability.md；与 04/06/07/08/09 对齐
