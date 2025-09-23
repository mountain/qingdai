# Project 015: 生态适应与光谱生物地理学（Adaptive Spectral Biogeography）

状态（2025‑09‑24）
- [x] 文档与方案定稿（本文件，v1）
- [ ] M1：气候稳态“环境档案”生成器（Environment Archive）
- [ ] M2：生态适应性规则集（陆地/海洋）与 API
- [ ] M3：动态“生物反照率反馈”模块（与辐射路径耦合）
- [ ] M4：“盖亚假说”对比实验（生物耦合 vs. 参照行星）
- [ ] M5：全球光谱生物地理图谱与可视化面板

关联模块与文档
- P004–P005：地形、海陆掩膜与基础反照率
- P006：能量收支（SW/LW，反照率与云耦合）
- P007/P011：海洋与海冰、风驱动浅水海洋与 SST
- P008：湿度与 E–P–LH 闭环、云—辐射一致性
- P009：水循环闭合（E–P–R），陆地贮水与雪
- P010：数值稳定与反噪（确保诊断稳健）
- P012：极点处理（避免极区伪差污染诊断）
- P013：Spin‑up 与平衡态重启（本项目依赖稳态多年平均）
- P014：地表水文与径流路由（河网/湖泊作为生境要素）
- docs/04-runtime-config.md：运行参数目录（本项目新增环境变量将纳入后续修订）

---

## 1. 背景与动机

我们将 GCM 从“纯物理气候系统”扩展为“生物—地球物理双向耦合”的思想实验平台，探索生命对地表反照率（base_albedo）的适应性调节如何反过来影响气候。本项目关注“光谱适应”的一阶效应：不同能量—水分—光照条件下，陆地植被与海洋浮游生物可能演化出不同的颜色与反射特性，从而改变短波反照率，进而反馈到能量收支。

目标并非建立完整生物地球化学模型，而是在当前 GCM 框架内，以最小但自洽的方式检验三条核心假说（§2），并给出可重复的对比实验 SOP。

---

## 2. 核心假说

1) 适应性反照率假说  
- 生命根据局地水、光、热条件演化颜色（反射光谱），从而改变基础反照率 base_albedo。能量充沛区倾向“更亮”（抑制灼伤/过热），能量稀缺区倾向“更暗”（提升吸收）。

2) 生物—气候反馈假说  
- 生物调节反照率会直接影响 SW 收支并改变局地/全球温度、云/海冰等，从而塑造其栖息环境（正/负反馈并存）。

3) 盖亚稳定性假说  
- 由生命驱动的反照率反馈可降低温度/海冰等关键变量的振幅，提升系统对天文强迫（如“脉冲季”）与外部扰动的稳健性。

---

## 3. 总体目标与里程碑

- M1 环境档案：从长期稳态模拟中统计“生境指标”（光/热/水/云/冰/河湖/海洋上升流代理等），生成标准化 NetCDF（Environment Archive）。
- M2 生态规则：在 `pygcm/ecology.py` 定义陆地/海洋的“适应性反照率”规则族与平滑器，输入环境档案输出“目标生物反照率地图”。
- M3 动态反馈：在主循环按“行星年”节拍（或指定步长）更新 base_albedo_map（平滑过渡与数值保护），将生物反照率纳入 P006 短波路径。
- M4 实验对比：以相同初值分别运行“参照行星（无生物反馈）”与“生命行星（启用 M3）”，统计对比“温度/海冰振幅、扰动恢复速度、TOA/SFC 收支”等指标。
- M5 可视化：输出“全球光谱生物地理图谱”（最终/阶段性）、关键指标时间序列与差值图；可选叠加河网/湖泊与上升流代理。

---

## 4. 架构与数据流

文本架构图（自上而下）  
- Spin‑up（P013）→ 平衡态重启  
- M1：环境档案生成（多年窗口平均）  
- M2：生态规则 → 目标生物反照率 α_bio*（陆/海分开，掩膜冰/沙漠/湖泊）  
- M3：动态更新 base_albedo_map ← α_bio*（时序平滑/空间平滑/夹持）  
- P006：短波 SW 路径使用新 α_total = α_surface(Ts,type; α_base_bio)·(1−C)+α_cloud·C  
- 诊断与可视化：对比“生命行星 vs. 参照行星”的能量/水量闭合与气候态指标  

集成点（主循环）
- 建议每“行星年”的末尾触发一次生态更新（或按 QD_BIO_UPDATE_EVERY_YEARS），实现上用步数累计与年长换算（参见 docs/02 的年定义；代码层面按 `orbital.T_p_days` / `QD_TOTAL_YEARS` 或统计循环步数近似）。

---

## 5. 环境档案（M1）规范

文件
- 输出：`data/ecology/environment_archive.nc`
- 维度：`lat(n_lat)`, `lon(n_lon)`
- 时间窗：近 N 年多年平均与振幅（默认 N=20，配置见环境变量）
- 坐标属性：与 GCM 网格一致；经度周期

字段（建议，单位 / 说明）
- 光（短波入射/云）  
  - `I_mean` (W m⁻²)：到达大气顶/地表的年平均几何入射或地表吸收（二选一，默认地表吸收）  
  - `I_season_amp` (W m⁻²)：季节性振幅（例如 95th−5th 年内分位差或年谐波幅值）  
  - `cloud_mean` (1)：年平均云量（0..1）  
  - `cloud_var` (1²)：云量方差（可选）
- 热（表面温度）  
  - `Ts_mean` (K)、`Ts_ann_amp` (K)：年平均与年振幅  
  - `Ts_dn_amp` (K)：日夜温差的年平均（可选）
- 水（降水/蒸发/径流/冰）  
  - `P_mean` (kg m⁻² s⁻¹ 或 mm day⁻¹)  
  - `E_mean` (...)  
  - `R_mean` (...)（P009，陆地径流；海洋设 0）  
  - `ice_frac_mean` (1)：海冰年平均覆盖（海洋/湖泊）
- 地表与地形  
  - `land_mask` (0/1)、`lake_mask` (0/1)（来自 P005/P014）  
  - `elevation` (m)（用于沙漠/高山代理，可选）
- 海洋代理  
  - `tau_curl_mean` (N m⁻³)：风应力涡度（上升流代理，源自 P011）  
  - `SST_grad_mean` (K m⁻¹)：SST 水平梯度（锋区/营养盐代理）
- 可选光谱分量（若 diag_isr 输出两星分量）：  
  - `I_A_frac`, `I_B_frac`：两星分量比值（0..1），参与颜色映射的可视化（物理上不直接入规则，首版可跳过）

统计方法
- 运行期在线累加（推荐）：每步更新滚动和/平方和/极值统计，年末计算均值与振幅；或后处理读取输出场（成本更高）。
- 有效权重：面积权重 w=cosφ；但档案存储逐格点本地指标。

---

## 6. 生态适应性规则集（M2）

设计原则
- 连续、可微/可平滑，避免硬阈导致格点伪迹
- 与 P006/P007/P008/P009 物理路径一致（冰/雪与海冰优先，不被生物反照率覆盖）
- 保守夹持：α ∈ [α_min, α_max] 且接近 P005 给出的底图范围；海洋变化幅度小于陆地

术语
- α_base_topo：来自 P005 的基础反照率底图（既有）
- α_bio*：生态规则给出的“目标生物反照率”
- α_base_bio：用于 SW 路径的“（时空平滑后的）生物基础反照率”

### 6.1 陆地规则（示意）

输入：`Ts_mean, Ts_ann_amp, P_mean, cloud_mean, R_mean, lake_mask, land_mask, elevation`  
输出：`α_bio*_land`

- 干湿指数（简化）：`W = tanh((P_mean − E_mean)/W_ref)`，E_mean 可由档案或 `P−R` 近似  
- 能量丰度指标：`H = tanh((I_mean − I_ref)/I_scale)` 或用 `Ts_mean` 代理  
- 季节性/气候严酷度：`S = tanh((Ts_ann_amp − S_ref)/S_scale)`  
- 规则（平滑加权）：
  - 湿暖 & 高光：倾向“变亮”防过热，α↑：目标 `α_hi_land`（典型 0.25–0.35）
  - 寒冷/低光：倾向“变暗”以吸收，α↓：目标 `α_lo_land`（典型 0.10–0.15）
  - 干旱（W 低）/沙漠：目标 `α_desert`（典型 0.30–0.40），受 elevation/纬度修正
- 混合：  
  `α_bio*_land = clip( w_warm*α_hi_land + w_cold*α_lo_land + w_dry*α_desert, α_min_land, α_max_land )`  
  其中权重由 (W, H, S) 的连续函数给出（tanh/sigmoid 组合，权重归一）

特殊掩膜
- 冰雪覆盖/寒带高山：保留 α_ice 或 P005 的高反照率，不被生态调节覆盖
- 湖泊（land_mask==1 且 lake_mask==1）：按“水体”或单独的 `α_lake` 处理（见 §6.3）

### 6.2 海洋规则（示意）

输入：`tau_curl_mean, SST_grad_mean, I_mean, ice_frac_mean`  
输出：`α_bio*_ocean`

- 上升流/营养盐代理：  
  `U = tanh((tau_curl_mean − U_ref)/U_scale) + tanh((SST_grad_mean − G_ref)/G_scale)`  
- 光照强弱：`L = tanh((I_mean − I_ref)/I_scale)`
- 规则：  
  - 低营养/高光（副热带环流中心区）：浮游生物稀少 → 维持接近 `α_water`（0.06–0.10）  
  - 上升流强/锋区（U 高）：高叶绿素 → “生物变暗”或“色度改变”导致 α 略降（幅度 ≤ 0.02）  
- 目标：  
  `α_bio*_ocean = clip( α_water − k_u*U + k_l*(1−L), α_min_ocean, α_max_ocean )`  
  并在 `ice_frac_mean>阈值` 区域退回 α_ice

### 6.3 湖泊（可选）

- 参数 `QD_TREAT_LAKE_AS_WATER`（P014）开启时，湖面在能量/湿度路径按水体处理  
- 本项目：`α_bio*_lake = α_lake`（默认 = α_water 或由环境变量覆盖）；若将湖泊纳入生态色谱，复用海洋规则

### 6.4 数值平滑与夹持

- 空间平滑：Gaussian σ≈1–2 格点（`QD_BIO_SMOOTH_SIGMA`）  
- 时间平滑：指数加权移动平均（EWMA）  
  `α_base_bio(t+)= (1−γ)·α_base_bio(t) + γ·α_bio*`，`γ=QD_BIO_TEMPORAL_ALPHA`  
- 全域夹持：`α ∈ [QD_BIO_ALPHA_MIN_{LAND/OCEAN}, QD_BIO_ALPHA_MAX_{LAND/OCEAN}]`

---

## 7. 动态“生物反照率反馈”（M3）

功能
- 在主循环中按年节律/设定周期更新 `α_base_bio` 并传入 P006 短波路径
- 保证与冰/云/海冰耦合一致：最终 `α_total` 仍按 docs/06 定义

触发时机
- 每当累计模拟时长达到 `QD_BIO_UPDATE_EVERY_YEARS` 的整数倍时：  
  1) 从“在线档案累积器”或 `environment_archive.nc` 读取指标  
  2) 计算 `α_bio*`（陆/海/湖）  
  3) 空间/时间平滑、夹持与掩膜处理  
  4) 替换/融合 `base_albedo_map` → 进入后续步的 `calculate_dynamic_albedo`

保护
- 禁止在海冰格点改变 α（以 α_ice 优先）  
- 限制一次更新的最大跃迁幅度（`QD_BIO_MAX_STEP`，如 ≤0.02）  
- 统一日志与图像输出：打印全局统计（mean/min/max）与示意图

---

## 8. 新增 API 与脚本（建议）

模块：`pygcm/ecology.py`
```python
@dataclass
class BiogeoParams:
    # 门限/尺度
    I_ref: float; I_scale: float
    W_ref: float; S_ref: float
    U_ref: float; G_ref: float
    # 目标范围与步进
    alpha_hi_land: float
    alpha_lo_land: float
    alpha_desert: float
    alpha_water: float
    alpha_lake: Optional[float]
    alpha_min_land: float; alpha_max_land: float
    alpha_min_ocean: float; alpha_max_ocean: float
    temporal_alpha: float  # γ
    smooth_sigma: float
    max_step: float
    update_every_years: float

def build_environment_archive_online(acc, fields, dt, grid) -> None:
    """在线累加器：滚动均值/方差/极值等，用于年末出档案或直接供生态规则读取。"""

def adaptive_albedo_from_archive(archive, base_albedo, masks, params) -> np.ndarray:
    """根据环境档案与规则生成 α_bio* 目标场（陆/海/湖分开处理）。"""

def update_adaptive_albedo(alpha_prev, alpha_target, masks, params) -> np.ndarray:
    """时间平滑与夹持，返回下一期 α_base_bio。"""
```

脚本（可选）
- `scripts/generate_environment_archive.py`：从输出文件离线汇总生成 `environment_archive.nc`
- `scripts/preview_biogeo_map.py`：读档案 + 规则 → 输出生物反照率地图与对比图
- `scripts/run_simulation.py`：读取 `QD_BIO_*` 环境变量，初始化生态模块与在线累加器，并在年末触发更新

---

## 9. 环境变量（建议，后续纳入 docs/04）

开关与节拍
- `QD_BIO_ENABLE`（默认 0）：启用生态适应与生物反照率反馈
- `QD_BIO_UPDATE_EVERY_YEARS`（默认 1）：更新周期（行星年）
- `QD_BIO_ARCHIVE_YEARS`（默认 20）：档案的多年平均窗口（年）

范围与夹持
- `QD_BIO_ALPHA_MIN_LAND`（默认 0.08）、`QD_BIO_ALPHA_MAX_LAND`（默认 0.45）
- `QD_BIO_ALPHA_MIN_OCEAN`（默认 0.06）、`QD_BIO_ALPHA_MAX_OCEAN`（默认 0.12）
- `QD_BIO_MAX_STEP`（默认 0.02）：单次时间更新的最大跃迁

目标与规则参数（典型值，需标定）
- `QD_BIO_ALPHA_HI_LAND`（默认 0.30）
- `QD_BIO_ALPHA_LO_LAND`（默认 0.12）
- `QD_BIO_ALPHA_DESERT`（默认 0.35）
- `QD_BIO_ALPHA_WATER`（默认 0.08）
- `QD_BIO_ALPHA_LAKE`（可选，不设则与 WATER 同）

平滑
- `QD_BIO_TEMPORAL_ALPHA`（默认 0.3）：时间平滑 γ
- `QD_BIO_SMOOTH_SIGMA`（默认 1.0 格点）：空间平滑尺度

阈值与尺度（经验起点）
- `QD_BIO_I_REF`、`QD_BIO_I_SCALE`（光照）
- `QD_BIO_W_REF`（干湿）、`QD_BIO_S_REF`（季节性）
- `QD_BIO_U_REF`、`QD_BIO_G_REF`（上升流/锋代理）

说明：所有默认值“温和不爆裂”，首轮实验按本文件运行示例执行，再按 M4 标定与回归确定默认组。

---

## 10. 运行示例

前置：生成/选择稳态重启文件（P013 SOP），并已具备 P004/P005 外部地形与 P006–P012 默认参数。
```bash
# 选择最新地形与启用能量/诊断（建议）
export QD_TOPO_NC=$(ls -t data/*.nc | head -n1)
export QD_USE_TOPO_ALBEDO=1
export QD_ENERGY_W=1
export QD_ENERGY_DIAG=1
```

A) 参照行星（无生物反馈）
```bash
unset QD_BIO_ENABLE
python3 -m scripts.run_simulation
```

B) 生命行星（启用生态适应与年更新）
```bash
export QD_BIO_ENABLE=1
export QD_BIO_UPDATE_EVERY_YEARS=1
export QD_BIO_ARCHIVE_YEARS=20
export QD_BIO_TEMPORAL_ALPHA=0.3
export QD_BIO_SMOOTH_SIGMA=1.0
export QD_BIO_MAX_STEP=0.02

# 可选调参（温和默认）
export QD_BIO_ALPHA_HI_LAND=0.30
export QD_BIO_ALPHA_LO_LAND=0.12
export QD_BIO_ALPHA_DESERT=0.35
export QD_BIO_ALPHA_MIN_LAND=0.08
export QD_BIO_ALPHA_MAX_LAND=0.45
export QD_BIO_ALPHA_MIN_OCEAN=0.06
export QD_BIO_ALPHA_MAX_OCEAN=0.12

python3 -m scripts.run_simulation
```

C) 扰动恢复对比（可选）
```bash
# 对两套实验在年边界施加相同扰动（示例：降低短波 2%）
export QD_SW_A0=0.06   # 正常
# 扰动窗口内
# export QD_SW_A0=0.0612    # 等效 2% 变化（或在脚本中实现统一开关）
```
随后对比两者 TOA/SFC 净通量回归速度、Ts 恢复时间常数、海冰面积恢复曲线等。

---

## 11. 验收标准（建议）

- 功能性：  
  - [ ] QD_BIO_ENABLE=1 时每“生态更新周期”均有 α_bio 日志打印与示意图，且 SW 路径使用新 α_base_bio  
  - [ ] 档案生成（在线或离线）字段完备、无 NaN/Inf、维度与网格一致

- 守恒与一致性（长期平均）  
  - [ ] 能量：TOA/SFC/ATM 净通量 |净| < 2 W m⁻²（与 P006 标准一致）  
  - [ ] 潜热：⟨LH⟩ ≈ ⟨LH_release⟩（与 P008 标准一致）  
  - [ ] 水量：⟨E⟩ ≈ ⟨P⟩ + ⟨R⟩（与 P009 一致；生态更新不破坏闭合）

- 对比实验（生命 vs. 参照）  
  - [ ] 全球年温度振幅（Ts_ann_amp）的多年代均值降低（阈值/幅度按参数扫描报告）  
  - [ ] 大扰动后的恢复时间（e‑folding）缩短  
  - [ ] 海冰面积年振幅下降、极端事件概率（尾部分位）降低

- 视觉与诊断  
  - [ ] 生物反照率地图与地表类型/水系/上升流区域呈合理空间格局  
  - [ ] 无条纹/块状伪迹（P010 combo 默认），极区无经向伪差（P012 开启）

---

## 12. 风险与注意事项

- 规则的“物理直觉”与“光谱真实”仍有差距；首版聚焦“反照率幅度”的一阶效应验证，后续可引入更物理的光谱模型与云—生物相互作用
- 规则过强可能引发不稳定（正反馈过度降低/升高 α），务必使用夹持与最大步长限制，并从温和默认起步
- 档案统计窗口过短会引入年际噪音；建议 ≥10–20 年
- 海冰/雪优先权不可破坏；冰下海洋不做生物变暗；湖泊可先与海洋同参
- 运行时间成本：生态更新周期间隔较长（年级别），总体开销小；离线档案生成应避免大文件多次 IO

---

## 13. 与其它模块的关系

- P006：α_base_bio 仅改变 SW 地表反照率基底，LW/SH/LH 路径保持一致；能量诊断阈值不放宽  
- P007/P011：SST 影响 E 与云，反过来也受 α 变化影响；上升流代理来自风应力涡度/SST 梯度  
- P008：云—辐射一致性未改变；若未来引入“生态—云”联动（如蒸腾/气溶胶），另立项目  
- P009：E–P–R 闭合与库容记账不修改；生态色谱不直接引入独立水库  
- P010/P012：仍需维持 combo 抑噪与极点修正以保证生态地图平滑与极区一致  
- P014：河网/湖泊图层用于生态可视化与“湖泊按水体处理”的一致性

---

## 14. 变更记录（Changelog）

- 2025‑09‑24：v1 文档定稿（架构/API/参数/验收/运行示例）；新增在线/离线两路档案方案；提出温和默认参数组与对比实验 SOP
