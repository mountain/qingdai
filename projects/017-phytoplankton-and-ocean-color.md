# Project 017: 浮游生物与海色（Phytoplankton & Ocean Color）

状态（2025‑09‑25）
- [x] 方案定稿（本文件）
- [ ] M1：最小浮游生物（Phyto）表层混合层模型 + 海色/带反照率耦合
- [ ] M2：营养盐与简化生长/呼吸/分裂（与 docs/16 Genes 参数一致）
- [ ] M3：海流平流与垂向浮力/混合（与 P011 动态海洋对接）
- [ ] M4：诊断与可视化（Ocean Color/TrueColor/Kd(490)/Chl‑a），参数标定
- [ ] M5：湖泊/内陆水体一致化（与 P014 路由/湖泊）

交叉引用
- docs/16‑phytoplankton‑model.md：水生个体/菌落级 Phyto 模型（Genes/浮力/吸收带/营养盐）
- docs/14‑adaptive‑spectroscopy.md：带离散、I(λ)→I_b、带反照率 A_b 与 TrueColor
- docs/07‑ocean‑and‑sea‑ice.md：混合层海洋/海冰、SST 平流与 Q_net 注入
- docs/06‑energy‑framework.md：短波/长波、地表反照率 α_total 的耦合
- docs/15‑ecology‑framework.md：Adapter/双时序接口；与陆生生态（P015）一致的回写策略
- projects/011‑ocean‑model.md：风驱动浅水海洋（uo/vo/η）与 SST 平流
- projects/014‑surface‑hydrology.md：湖泊掩膜与水体类型一致化

---

## 1. 目标与范围

- 在 Qingdai GCM 中引入“水生自养浮游生物（phyto）”的最小但自洽模型，使其：
  1) 作为表层混合层内的生物量场（Chl‑a 代理）随光照/温度/营养在日尺度演化；
  2) 通过“水体光学—海色—反照率”的路径改变短波吸收与地表能量收支（海—气耦合）；
  3) 受海流（P011）平流/扩散与简化垂向过程（浮力/混合）控制空间分布；
  4) 产出观测相容的海色产品（TrueColor/带反照率、Kd(490)、Chl‑a），用于可视化与验证。

- 最小闭环（M1–M2）：
  I_b(surf) → Kd_b(Chl) → Ī_b(MLD) → μ(I,T,N) → C_phyto → R_rs(λ)/A_b^water(Chl) → α_total → SW_sfc

- 非目标（本里程碑外）：
  - 完整三维生物地球化学（多营养盐、浮游动物摄食等）；
  - 逐纳米/多重散射辐射传输；
  - 海冰下藻类与复杂季节学。

---

## 2. 概念与变量（混合层平均方案）

### 2.1 状态与单位
- C_phyto（mg Chl m⁻³）表层混合层体积浓度（预报量）
- B_phyto（mg Chl m⁻²）面密度：B = C_phyto · H_mld（诊断/用于光学）
- N（mmol m⁻³，可选）营养盐浓度（M2 起），单一限制（如 DIN）
- T_w（K）表层水温（来自 P007/011，用于温度生理因子）
- H_mld（m）混合层深度（复用 `QD_MLD_M` 或 `QD_OCEAN_H_M`）

说明：M1 采用“单层均匀混合”的 0D‑in‑MLD 近似；M3 引入水平平流/扩散与简化垂向项。

### 2.2 光谱—海色近似（带离散）
- 地表带强度 I_b(surf)（docs/14）：双星 + 大气调制后按 NB 个短波带离散。
- 水体辐照度在混合层的带平均（指数衰减）：
  Ī_b = (1/H_mld) ∫_0^{H_mld} I_b(surf) · e^{−Kd_b z} dz
      = I_b(surf) · (1 − e^{−Kd_b H_mld}) / (Kd_b H_mld)
- 带消光系数 Kd_b（单参 Chl 近似）：
  Kd_b(Chl) = Kd0_b + k_chl_b · Chl^m
  其中 Kd0_b 为本征水体/溶解有色物背景，k_chl_b、m 为经验参数（典型 m≈0.5）。

---

## 3. 生长—损耗模型（M1–M2）

### 3.1 净增长率 μ（d⁻¹）
- 光限制（带整合；与 docs/13/14 接口一致）：
  - 简化缓饱和：f_I = tanh( (α_P · Ī_PAR) / μ_max )
  - 或吸光曲线：f_I = 1 − exp( − α_P · Σ_b Ī_b · Δλ_b / (μ_max · I_ref) )
- 温度因子（Q10 近似）：f_T = Q10^{(T_w − T_ref)/10}
- 营养限制（M2）：f_N = N / (K_N + N)

净增长：μ = μ_max · f_I · f_T · f_N − m_0
- μ_max：最适生长率；α_P：光利用系数；m_0：背景损耗（呼吸/死亡/微食）

### 3.2 预报方程（混合层平均）
- M1（无平流/扩散）：
  dC_phyto/dt = μ · C_phyto − λ_sink · C_phyto / H_mld
  其中 λ_sink 为等效“下沉/流失”率（可设 0–小值）
- M2（含营养）：
  dN/dt = − (Y^{-1} · μ_grow · C_phyto) + R_remin
  其中 μ_grow = μ_max · f_I · f_T · f_N，Y 为同化产量，R_remin 简化再矿化（可设 0）

数值：与生态/物理日步对齐（默认 1 天）；若与物理步耦合，采用子步缓存/日末聚合（见 §8）。

---

## 4. 海色—反照率耦合（与 docs/06/14 对齐）

### 4.1 水体带反照率 A_b^water(Chl)（经验近似）
- 目标：将 Chl（或 C_phyto）映射为短波带反射（绿增强、蓝减弱）。
- 最小形式（单参幂律 + 谱形）：
  A_b^water(Chl) = A_b^pure + c_reflect · S_b · Chl^{p_reflect}
  - A_b^pure：本征水体的带反射（常数或查表）
  - S_b：谱形系数（绿/黄带为正、蓝带为负或 0；由 docs/14 带定义给出）
  - c_reflect, p_reflect：经验参数（默认 c_reflect≈0.02、p_reflect≈0.5；可标定）
- 限幅：clip A_b^water ∈ [A_min_b, A_max_b]

可选更物理化（后续）：由 a(λ), b_b(λ) → R_rs(λ) → 带平均 A_b^water。

### 4.2 与总反照率 α_total 合成（docs/06）
- 海洋网格处，将短波带级“水体反照率”设为 A_b^water(Chl)：
  α_water_eff(b) ← A_b^water(Chl)
- 若能量模块仍使用标量 α_water，则降维：
  α_water_eff = Σ_b A_b^water · w_b^SW（w_b^SW 可用当前 I_b 比例或固定权重）
- 最终 α_total = α_surface(type, Ts, ice) · (1 − C_cloud) + α_cloud · C_cloud
  其中 type = ocean 时使用 α_water_eff。

说明：仅改变地表短波反射，不改动大气参数；能量闭合由 docs/06 统一诊断。

---

## 5. 传输与垂向（M3）

### 5.1 水平平流与扩散（与 P011）
- 在 ocean.step 之后，以海流（uo, vo）平流：
  ∂C/∂t = − (uo · ∇)C + K_h_phy ∇² C
- K_h_phy（m² s⁻¹）：生物量水平混合系数（默认与 `QD_KH_OCEAN` 同阶，或独立配置）

### 5.2 垂向（最小）
- M1 忽略显式垂向；λ_sink 表示等效下沉/流失
- M3 可加入“浮力因子”（docs/16 `genes.buoyancy_factor`）+ 混合项对 C_phyto 的缓慢调制（仍保持混合层均匀假设）

### 5.3 湖泊（P014）
- `QD_TREAT_LAKE_AS_WATER=1` 时，湖面采用与海面相同的 A_b^water 公式，但允许独立参数（`QD_PHYTO_*_LAKE`）；Chl 动力学同海洋但默认无平流（或极弱）。

---

## 6. 架构与 API（建议）

### 6.1 模块
- 新增 `pygcm/ecology/phyto.py`：
  - class PhytoManager(grid, nbands, band_defs, params)
    - `step_subdaily(weather_inst, dt_seconds)`（可选：收集 I_b(t)、即时诊断）
    - `step_daily(weather_day)` → 更新 C_phyto（与 N）、计算 A_b^water；返回 `PhytoDailyReport`
    - `advect_diffuse(uo, vo, dt_seconds)`（M3）
    - `optics_bands(chl)` → (A_b^water[NB], Kd_b[NB])
    - `diagnostics()` → {`Chl_mean`, `Kd490`, `alpha_water_eff`, …}

- 报告/数据结构（与 docs/15/13 风格一致）：
  - `PhytoDailyReport`：
    - `chl_mean`, `chl_map`（可选降采样）、`alpha_water_bands`, `alpha_water_scalar`
    - `kd_bands`, `kd_490`, `productivity_day`（可选）

### 6.2 主循环集成（scripts/run_simulation.py）
- 初始化：
  - 读取 `QD_PHYTO_ENABLE` 与光谱带定义（与 docs/14 一致）
  - 构造 `PhytoManager`（NB、带边界/权重、默认参数）
- 每物理步：
  1) 计算/缓存 I_b(surf)、云场
  2) （可选）若到达子采样频率，调用 `phyto.step_subdaily` 收集诊断
  3) 海洋步 `ocean.step` 之后（若 M3 启用），调用 `phyto.advect_diffuse`
- 日界：
  1) 由 Adapter 聚合 `WeatherDaily`（Ī_b、T_w 等）
  2) `phyto.step_daily(weather_day)` → 生成 `alpha_water_bands`
  3) 以 `QD_PHYTO_FEEDBACK_MODE` 决定“即时/日末”回写短波带 α_water；计算 `α_total`
  4) 记录/出图：Chl、Kd(490)、OceanColor/TrueColor 叠加

---

## 7. 环境变量（建议；与 docs/04 对齐）

主控
- `QD_PHYTO_ENABLE`（默认 0）：开启/关闭浮游生物模块
- `QD_PHYTO_DT_DAYS`（默认 1.0）：日级更新步长
- `QD_PHYTO_SUBDAILY_ENABLE`（默认 0）：启用小时级接口（收集诊断/可选即时回耦）
- `QD_PHYTO_SUBSTEP_EVERY_NPHYS`（默认 6）：每 N 个物理步一次子步
- `QD_PHYTO_FEEDBACK_MODE`（instant|daily，默认 daily）：带 α_water 回写策略
- `QD_PHYTO_ALBEDO_COUPLE`（默认 1）：是否将 A_b^water 回写到短波

光谱/带
- `QD_PHYTO_SPECTRAL_BANDS`（默认与 `QD_ECO_SPECTRAL_BANDS` 一致，建议 8–16）
- `QD_PHYTO_SPECTRAL_RANGE_NM`（默认 380,780）
- `QD_PHYTO_TOA_TO_SURF_MODE`（simple|rayleigh|custom，默认沿用 docs/14 配置）

光学参数（带数组或缩放系数）
- `QD_PHYTO_KD0`（默认按带/经验表）
- `QD_PHYTO_KD_CHL`（默认按带/经验表）
- `QD_PHYTO_KD_EXP_M`（默认 0.5）
- `QD_PHYTO_REFLECT_C`（默认 0.02）
- `QD_PHYTO_REFLECT_P`（默认 0.5）
- `QD_PHYTO_SPECTRAL_SHAPE`（路径或预设：blue/green/yellow 权重表）

生理与损耗
- `QD_PHYTO_MU_MAX`（d⁻¹，默认 1.5）
- `QD_PHYTO_ALPHA_P`（默认 0.04）
- `QD_PHYTO_Q10`（默认 2.0）
- `QD_PHYTO_T_REF`（K，默认 293.15）
- `QD_PHYTO_M_LOSS`（d⁻¹，默认 0.05）
- `QD_PHYTO_LAMBDA_SINK`（m d⁻¹，默认 0–1）
- `QD_PHYTO_I_REF`（W m⁻²，默认 50；若使用吸光式 f_I）

营养盐（M2）
- `QD_PHYTO_YIELD`（mg Chl / mmol N，默认 1）
- `QD_PHYTO_KN`（mmol m⁻³，默认 0.5）
- `QD_PHYTO_REMIN`（mmol m⁻³ d⁻¹，默认 0）

传输（M3）
- `QD_PHYTO_KH`（m² s⁻¹，默认与 `QD_KH_OCEAN` 同阶）
- `QD_PHYTO_ADVECTION`（0/1，默认 1）
- `QD_PHYTO_BUOYANCY_SCALE`（默认 0，启用后按 genes 浮力微调垂向损耗项）

诊断与可视化
- `QD_PHYTO_DIAG`（默认 1）：打印诊断
- `QD_PLOT_OCEANCOLOR`（默认 1）：输出 Ocean Color/TrueColor 图（含海色）
- `QD_OC_KD_BAND_REF_NM`（默认 490）：Kd(490) 参考波段 nm 值（近似选最近带）

湖泊（P014）
- `QD_PHYTO_LAKE_ENABLE`（默认 1）
- `QD_PHYTO_REFLECT_C_LAKE`、`QD_PHYTO_REFLECT_P_LAKE`、`QD_PHYTO_SPECTRAL_SHAPE_LAKE`（可选覆盖）

---

## 8. 运行时序与数值注意

- 推荐时序（与 docs/12 一致）：
  1) 轨道/短波入射 → I_b(surf)（docs/14）
  2) 动力/湿度/能量 → 海洋步（docs/06/07/08/10/11）
  3) Phyto 传输（M3）：`advect_diffuse`
  4) 日界：`step_daily` → A_b^water → α_total 回写
  5) 水文/路由（P009/P014）→ 诊断与出图

- 与能量闭合：生态回写仅改变 α_water（短波）；TOA/SFC/ATM 诊断（docs/06）在长期平均应维持近守恒。

- 数值稳定：
  - f_I 采用缓饱和（tanh）形式更稳健；
  - μ、Kd 随 Chl 的幂律关系建议限幅；
  - 传输步用 CFL‑安全的半拉氏或稳定插值/扩散；
  - 带/权重在初始化时固定，降低日常重算成本。

---

## 9. 诊断与验收标准

功能性
- [ ] `QD_PHYTO_ENABLE=1` 后可稳定运行，产出 Chl、Kd(490)、Ocean Color/TrueColor、`alpha_water_eff`
- [ ] `QD_PHYTO_FEEDBACK_MODE=instant/daily` 生效；能量诊断无异常漂移

物理合理
- [ ] 高光区（清澈热带）Chl 低、蓝带强反射；富营养/冷水上升区 Chl 高、绿带增强
- [ ] Kd(490) 与 Chl 呈正相关；海色随季节/环流展现合理时空变化

守恒与一致
- [ ] 多年平均 |⟨TOA_net⟩|、|⟨SFC_net⟩|、|⟨ATM_net⟩| < 2 W m⁻²（docs/06）
- [ ] 与 P011 的传输耦合稳定，无网格噪音或爆裂

性能
- [ ] NB=8–16 时，额外开销可控（<10–15%）

---

## 10. 运行示例

仅诊断（不回写短波）
```bash
export QD_PHYTO_ENABLE=1
export QD_PHYTO_ALBEDO_COUPLE=0
export QD_PHYTO_DT_DAYS=1
export QD_PHYTO_SPECTRAL_BANDS=8
export QD_PHYTO_DIAG=1
python3 -m scripts.run_simulation
```

开启海色回耦（日末回写）
```bash
export QD_PHYTO_ENABLE=1
export QD_PHYTO_ALBEDO_COUPLE=1
export QD_PHYTO_FEEDBACK_MODE=daily
export QD_PHYTO_SPECTRAL_BANDS=8
python3 -m scripts.run_simulation
```

即时回耦（时级子步，每 N 个物理步重算一次）
```bash
export QD_PHYTO_ENABLE=1
export QD_PHYTO_ALBEDO_COUPLE=1
export QD_PHYTO_FEEDBACK_MODE=instant
export QD_PHYTO_SUBDAILY_ENABLE=1
export QD_PHYTO_SUBSTEP_EVERY_NPHYS=6
export QD_PHYTO_SPECTRAL_BANDS=16
python3 -m scripts.run_simulation
```

启用传输（与海流耦合），并输出 Ocean Color/TrueColor
```bash
export QD_PHYTO_ENABLE=1
export QD_PHYTO_ALBEDO_COUPLE=1
export QD_PHYTO_FEEDBACK_MODE=daily
export QD_PHYTO_ADVECTION=1
export QD_PHYTO_KH=5000
export QD_PLOT_OCEANCOLOR=1
python3 -m scripts.run_simulation
```

湖泊支持（按水体处理）
```bash
export QD_PHYTO_ENABLE=1
export QD_TREAT_LAKE_AS_WATER=1
export QD_PHYTO_LAKE_ENABLE=1
python3 -m scripts.run_simulation
```

---

## 11. 里程碑与实施要点

- M1：混合层 0D 生物量 + 海色回耦
  - 完成 PhytoManager；f_I（tanh）+ Q10；A_b^water 经验映射；标量/带回写链路
- M2：营养盐闭环
  - 引入 N 与 μ_grow、Y、K_N、R_remin；确保长期稳定不爆裂
- M3：平流/扩散与简化垂向
  - 半拉氏/稳定插值；与 P011 的 CFL/反噪参数对齐；λ_sink/浮力因子
- M4：Ocean Color 可视化与标定
