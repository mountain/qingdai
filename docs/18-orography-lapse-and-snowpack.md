# 18. 地形递减率、雪线与雪被水库：从山地降温到稳定基流（Orographic Lapse Rate, Snowline and Snowpack Reservoir）

状态（2025‑09‑27）
- [x] 知识框架与参数化方案
- [x] 与现有 P004/P005/P006/P009/P014 的接口设计
- [x] 环境变量建议与运行示例
- [ ] 代码落地（见 projects/019-orography-lapse-and-snowpack.md）

关联
- 地形与或ographic：docs/05、projects/004/005
- 能量收支与反照率：docs/06
- 湿度/云与相态：docs/08
- 水量闭合与径流/桶：docs/09
- 在线路由与河网/湖泊：docs/14
- 运行参数目录：docs/04（本章新增变量将补充到 04）

---

## 0. 背景与动机

当前模型已具备：
- 可选地形抬升降水增强（`QD_OROG`；docs/05 §5.6）
- 降水相态阈值与固定融雪率（`QD_SNOW_THRESH`、`QD_SNOW_MELT_RATE`；docs/04 §6、docs/09）

但仍缺少对“海拔导致的气温递减（lapse rate）与雪线”这个一阶山地效应的显式考虑。其后果包括：
- 山地高程区的大气/地表温度偏高，固相降水/雪被低估；
- 冬季固相（雪）作为“季节性水库”未能有效蓄水/释水，春融基流不稳；
- 反照率反馈（积雪高反照率）与光照/能量收支的地形耦合不足。

本章提出一个轻量、数值稳健、与现有模块一致化的参数化方案：引入“海拔温度递减 → 雪线/固相分配 → 雪被水库（SWE）→ 融雪出流”的闭环，并与动态反照率和水文路由协同，产生稳定的季节性基流（baseflow）。

---

## 1. 海拔温度递减（Lapse Rate）与“修正温度”

经验与物理依据：
- 干绝热递减率 Γ_dry ≈ 9.8 K/km，湿绝热 Γ_moist ≈ 4–7 K/km（取决于 T/q）
- 气候学平均常用“环境递减率” Γ_env ≈ 6.5 K/km（对全球山地温度有较好一阶刻画）

参数化（每格点 i）：
- 令海拔高度 `H_i = elevation(i)`（m）
- 取参考海拔 `H_ref`（默认 0 m）与递减率 `Γ`（K/km）
- 定义“海拔修正空气温度”：
  ```
  T̂_a(i) = T_a(i) − Γ · (H_i − H_ref)/1000
  ```
- 对地表温度可选作同样修正以决定相态/雪反照率门控：
  ```
  T̂_s(i) = T_s(i) − Γ_s · (H_i − H_ref)/1000    # 默认为 Γ_s = Γ
  ```
- 夜侧/极端保护：在能量路径已有 `T_FLOOR` 保护，修正后可再限幅：
  ```
  T̂_a ← max(T̂_a, T_floor),  T̂_s ← max(T̂_s, T_floor)
  ```

建议：将 `T̂_a` 用于“降水相态与融雪”，`T̂_s` 用于“积雪反照率与可视化雪层”。

几何与边界约束（新增）
- 有效地形高度上限：为避免不现实的极端地形，采用 `H_eff = min(H_bedrock + h_ice_eff, QD_LAND_ELEV_MAX_M)`，默认 `QD_LAND_ELEV_MAX_M≈10000 m`。  
- 极区冰厚上限：在极区（`|φ| ≥ QD_POLAR_LAT_THRESH`，默认 60°）对雪/冰几何厚度采用 `h_ice_eff ≤ QD_POLAR_ICE_THICK_MAX_M`（默认 ≈4500 m），代表近地球的极冠厚度上限。  
- 几何换算：雪的几何厚度 `h_snow ≈ SWE / QD_RHO_SNOW`（默认 `QD_RHO_SNOW≈300 kg·m⁻³`），用于上式中的 `h_ice_eff` 计算。

物理上限与几何约束（新增）
- 陆地有效海拔上限：当考虑冰雪/冰盖叠加时，用于 lapse/snowline 的“有效海拔”定义为
  H_eff = min(H_bedrock + H_ice_eff, QD_LAND_ELEV_MAX_M)。
  其中 H_ice_eff 可由（i）雪被厚度（SWE 转换为几何厚度，近似 h_snow ≈ SWE_mm·1e-3·ρ_w/ρ_snow，ρ_snow 典型 300 kg/m³）与（ii）未来可能的“大陆冰盖厚度”字段（若存在）组成。
- 极地冰川厚度上限：在极区（|φ| ≥ QD_POLAR_LAT_THRESH）对冰盖/冰川几何厚度实施上限 H_ice_eff ≤ QD_POLAR_ICE_THICK_MAX_M（默认采用与地球同量级的上限，见 §6）。
- 以上约束确保“冰雪叠加后的陆地海拔”不超过 10 km，且极地冰川/冰盖厚度不超过地球量级上限，用于 lapse/雪线与相态判定的一致性几何基准。

---

## 2. 雪线与降水相态分配

定义“雪线高度” H_snowline 满足 T̂≈T_thresh（典型 273.15 K）：
```
H_snowline ≈ H_ref + 1000 · (T_a,ref − T_thresh)/Γ
```
在模型中不必显式构造高度场上的等值线，可直接逐格点按“修正温度”分配相态。

平滑的固/液分配函数（避免硬阈值齿状伪迹）：
- 设过渡半宽 ΔT（典型 1.5 K）与 Sigmoid：
  ```
  f_snow(i) = σ( (T_thresh − T̂_a(i)) / ΔT ) = 1 / (1 + exp((T̂_a−T_thresh)/ΔT))
  ```
- 则在格点 i 的降水（质量通量）分配为：
  ```
  P_snow(i) = f_snow(i) · P_total(i)
  P_rain(i) = (1 − f_snow(i)) · P_total(i)
  ```

这将自然地在山地高程区/寒区形成更高的固相比例，从而在 3 章进入雪被库。

---

## 3. 雪被水库（SWE）与融雪出流

引入每陆地格点的“雪当量水层” SWE(i)（单位：mm 水层或 kg/m²；二者仅差 ρ_w）：

质量守恒（每日/物理步离散）：
```
SWE_next = SWE + P_snow · dt − M_snow · dt − Subl_snow · dt
```
- 本版最小实现忽略 Sublimation（可后续加入温度/风依赖），默认 `Subl_snow≈0`
- 融雪 `M_snow` 采用度日法（degree‑day）或回退至常数融化率：

度日融雪（推荐）：
```
M_snow = DDF · max( T̂_a − T_melt, 0 )         # DDF 单位：mm / (K·day)，T_melt≈273.15 K
```
常数融化率（与现有 `QD_SNOW_MELT_RATE` 对齐；mm/day）：
```
M_snow = M0 · 𝟙[T̂_a ≥ T_melt]
```
出流去向：
- 融雪 `M_snow` → 进入陆地桶 `W_land`（docs/09 §9.4），再按桶/径流时标排出
- 可选“快速面融分量” φ_fast（0..1），将一部分直接作为快流进入路由（默认 0）

水库限幅与稳定性保护：
- `SWE ← clip(SWE, 0, SWE_MAX)`（可选上限），负值截为 0
- 连续时间步内，`M_snow · dt ≤ SWE`（安全地 min）

---

## 4. 积雪反照率与能量/可视化耦合

雪具有高反照率（0.6–0.9），对能量收支与可视化影响显著。建议引入“积雪覆盖率” C_snow(i)（0..1）调制地表反照率：

最小方案（基于 SWE 的光学厚）：
```
C_snow = 1 − exp( − SWE / SWE_ref )      # SWE_ref 典型 10–20 mm
α_surface_eff = α_base · (1 − C_snow) + α_snow · C_snow
```
- `α_snow`：新雪反照率（默认 0.70）；可选“雪龄衰减” τ_alb（天）：
  ```
  α_snow(t+Δt) = α_fresh − (α_fresh − α_old) · (1 − exp(−Δt/τ_alb))
  ```
- 该 `α_surface_eff` 在能量模块（docs/06）替换原有陆面基反照率通道，仅对陆面起作用；海冰反照率仍由 P007 路径处理

TrueColor 可视化：与 `QD_TRUECOLOR_SNOW_BY_TS` 一致，可切换为“按 SWE/C_snow 渲染”，避免仅靠地表温度判定时的误检。

---

## 5. 与现有模块的接口与顺序

建议主循环（简化，突出新增步骤）：
1) 动力/湿度/云：得到 `P_total`、`T_a`、`T_s`、云量等
2) Lapse 修正：`T̂_a = T_a − Γ·ΔH`；（可选）`T̂_s = T_s − Γ_s·ΔH`
3) 相态分配：`P_snow, P_rain`（Sigmoid 分配）
4) 雪被步（本章新增）：
   - `SWE ← SWE + P_snow·dt − M_snow·dt`（M_snow 由 DDF 或常数）
   - `W_land ← W_land + (1−φ_fast)·M_snow·dt`；快流分量 → 路由输入缓存
   - 诊断：`⟨SWE⟩`、`⟨P_snow⟩`、`⟨M_snow⟩`
   - 冰盖掩膜与分流（新增）：定义陆地冰盖掩膜  
     `glacier_mask = (land=1) ∧ (C_snow ≥ QD_GLACIER_FRAC ∨ SWE ≥ QD_GLACIER_SWE_MM)`；其中 `QD_GLACIER_FRAC` 默认 0.60、`QD_GLACIER_SWE_MM` 默认 50 mm。  
     • 冰盖像元上的“雨”视为冻结沉积：`SWE ← SWE + P_rain_glacier·dt`，不进入陆地桶。  
     • 冰盖融水不入桶，直接作为“下游源项”进入路由（可理解为冰下管网）；非冰盖像元仍采用“入桶→线性径流→路由”的路径。  
     • 生态掩膜：冰盖像元 `soil_idx=0`、`LAI=0`，并在个体抽样（IndividualPool）与日级聚合中跳过冰盖像元以节省计算。
5) 反照率合成：计算 `C_snow(SWE)` → `α_surface_eff` 并与云/海冰合成总反照率（docs/06）
6) 能量步：短/长波 + SH/LH；海洋步（docs/07/11）
7) 水文步：陆地桶/径流（docs/09）；到水文步长时路由（docs/14）
8) 可视化：状态图/TrueColor 叠加积雪/雪线诊断（可选）

与（现有）变量对齐：
- 相态阈值 `QD_SNOW_THRESH` 继续使用，但改以 `T̂_a` 判定
- 固定融雪率 `QD_SNOW_MELT_RATE` 在未设 DDF 时作为后备

---

## 6. 环境变量建议（汇总，最终以 docs/04 收录为准）

主控
- QD_LAPSE_ENABLE（默认 1）：开启海拔温度递减
- QD_LAPSE_K_KPM（默认 6.5）：环境递减率 Γ，单位 K/km
- QD_LAPSE_KS_KPM（默认 = QD_LAPSE_K_KPM）：地表温度修正的 Γ_s

相态与雪线
- QD_SNOW_THRESH（默认 273.15 K）：雨/雪温阈
- QD_SNOW_T_BAND（默认 1.5 K）：相态 Sigmoid 过渡半宽 ΔT

雪被与融雪
- QD_SWE_ENABLE（默认 1）：开启雪被水库
- QD_SWE_INIT_MM（默认 0）：初始 SWE（mm）
- QD_SWE_MAX_MM（可选，默认不设）：SWE 上限（mm）
- QD_SWE_REF_MM（默认 15）：C_snow 的参考厚度（mm）
- QD_SNOW_ALBEDO_FRESH（默认 0.70）
- QD_SNOW_ALBEDO_OLD（默认 0.45）
- QD_SNOW_ALBEDO_DECAY_DAYS（默认 10）
- QD_SNOW_MELT_MODE（degree_day|constant，默认 degree_day）
- QD_SNOW_DDF_MM_PER_K_DAY（默认 3.0）：DDF（mm/K/day）
- QD_SNOW_MELT_RATE（默认 5.0 mm/day）：常数融雪率（后备）
- QD_SNOW_MELT_TREF（默认 273.15 K）：融雪起始温度
- QD_SNOW_FASTFLOW_FRAC（默认 0.0）：融雪快流比例 φ_fast

可视化与诊断
- QD_PLOT_SNOWLINE（默认 1）：绘制雪线诊断（见 7 节）
- QD_TRUECOLOR_SNOW_BY_SWE（默认 1）：TrueColor 按 SWE 渲染雪

几何上限与冰盖厚度（新增）
- QD_LAND_ELEV_MAX_M（默认 10000）：冰雪/冰盖叠加后的陆地“有效海拔”上限（m），用于 lapse/snowline/相态计算。
- QD_POLAR_ICE_THICK_MAX_M（默认 4500）：极地（|φ| ≥ 阈值）冰川/冰盖厚度上限（m），用于计算 H_ice_eff。
- QD_POLAR_LAT_THRESH（默认 60）：极地纬度阈值（deg），用于应用极地厚度上限。

---

## 7. 诊断与验收标准

必备诊断（建议每 N 步打印/输出）：
- [SnowDiag] `⟨SWE⟩`、`⟨P_snow⟩`、`⟨M_snow⟩`、`雪面覆盖率 ⟨C_snow⟩`
- 雪线统计：按纬带输出“`H_snowline(φ)` 或 `f_snow(φ)`”代表值（或以等值线叠加到状态图）
- 能量闭合（docs/06）：开启雪后多年平均 |⟨TOA/SFC/ATM_net⟩| < 2 W/m²
- 水量闭合（docs/09）：长期平均 `⟨E⟩ ≈ ⟨P⟩ + ⟨R⟩`；雪库趋势 `d⟨SWE⟩/dt` ~ 0
- 基流与季节性：中高纬/山地流域的“春融季节基流”显现且平滑，无非物理锯齿

可视化建议：
- 状态图叠加“`C_snow≥0.2` 等值域”为雪盖，或直接绘 SWE 等值线
- TrueColor：启用 `QD_TRUECOLOR_SNOW_BY_SWE`，对 C_snow≥阈值增强白度
- 河网/湖泊叠加（docs/14）：春融期间主干流量明显增强

---

## 8. 数值与单位一致性

单位
- 降水/融雪通量：kg m⁻² s⁻¹ ↔ mm/day（换算：1 mm/day ≈ 1.1574e-5 kg m⁻² s⁻¹）
- SWE（库存）：kg m⁻² ↔ mm（ρ_w=1000 kg/m³，1 mm = 1 kg/m²）

稳定性
- 相态分配使用 Sigmoid，避免硬阈值导致“棋盘雪线”
- 融雪以 `min(SWE/Δt, M_snow)` 限幅，避免负库存
- 反照率限幅 [0,1]；与云/海冰/生态叠加时保持凸组合

---

## 9. 与相关模块的耦合要点

- P006 能量：通过 `α_surface_eff` 与 `T_floor` 同步保护；雪盖提升 α → 地表冷却 → 可能延迟融雪（合理反馈）
- P008 湿度：积雪覆盖区域蒸发减弱（现已通过陆面蒸发缩放）；可后续在雪盖区进一步衰减 E
- P009 水量：新增雪库 `S_snow` 进入闭合核算；`M_snow` 进入陆桶 R 链
- P014 路由：春融期入河流量增强；快流分量 φ_fast（若>0）直接加入路由输入缓存
- P005 orographic：地形降水增强 + lapse 相态将共同塑形更现实的山地雪被/河网

---

## 10. 运行示例

基础：开启 lapse + 雪被 + 度日融雪 + 诊断
```bash
# 海拔递减
export QD_LAPSE_ENABLE=1
export QD_LAPSE_K_KPM=6.5

# 相态与雪被
export QD_SWE_ENABLE=1
export QD_SNOW_THRESH=273.15
export QD_SNOW_T_BAND=1.5
export QD_SWE_REF_MM=15

# 融雪：度日法
export QD_SNOW_MELT_MODE=degree_day
export QD_SNOW_DDF_MM_PER_K_DAY=3.0
export QD_SNOW_MELT_TREF=273.15

# 雪反照率与可视化
export QD_SNOW_ALBEDO_FRESH=0.70
export QD_SNOW_ALBEDO_OLD=0.45
export QD_SNOW_ALBEDO_DECAY_DAYS=10
export QD_TRUECOLOR_SNOW_BY_SWE=1
export QD_PLOT_SNOWLINE=1

# 运行
python3 -m scripts.run_simulation
```

强化山地效应：叠加地形降水增强
```bash
export QD_OROG=1
export QD_OROG_K=7e-4
python3 -m scripts.run_simulation
```

---

## 11. 未来扩展与研究方向

- 湿绝热/温湿度一致的 lapse（随 q、T 动态调整 Γ）
- 能量平衡融雪（弥散短波、长波、湍流通量、雨致融雪/冷却）
- 雪密度与压实、积雪热传导（多层雪）
- 亚网格尺度 hypsometry：将格点内地形分布用于精细雪线/雪盖分数估计
- 与生态/植被（docs/15）耦合：积雪遮盖期改变 LAI、推迟物候；反照率/光谱反馈联动

---

## 12. 验收清单（建议）

- [ ] lapse 修正生效：中高山地区冬季 T̂ 明显降低，雨/雪分配合理
- [ ] 雪被库守恒：多年平均 d⟨SWE⟩/dt ≈ 0（非单调漂移）
- [ ] 融雪期基流：河网主干春季平稳增大，无数值锯齿
- [ ] 能量/水量闭合：阈值内（docs/06/09 标准）
- [ ] 反照率/可视化：雪盖区域 α_total 与 TrueColor 表现一致
- [ ] 诊断面板：SnowDiag/Snowline 输出齐备，可复现实验
