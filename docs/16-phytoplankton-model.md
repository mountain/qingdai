# 文档 16：浮游生物模型（Phytoplankton Model, P017）

状态（2025‑09‑27）
- [x] M1：混合层 0D 生物量 + 光谱—海色—短波回耦（多物种）  
- [x] M2：单营养盐池（N）竞争（Monod），再矿化，持久化/导入导出  
- [x] M3：海流平流/水平扩散（与 P011 对接）  
- [x] M4：诊断与可视化（Ocean Color/TrueColor/Kd(490)/Chl‑a）  
- [ ] M5：湖泊一致化细节与参数标定（延续）

交叉引用
- 文档 14：适应性光谱（带离散、I(λ)→I_b、权重）  
- 文档 06/07/11：能量收支、混合层海洋与风驱动海洋  
- 项目 017：浮游生物与海色（路线图、运行示例）

---

## 0. 目标与范围

在 Qingdai GCM 中引入“多物种浮游植物”的最小但自洽模型，使其：
- 作为混合层（MLD）平均的生物量场 C_i（mg Chl m⁻³）随“光谱光照/温度/营养盐”在日尺度演化；
- 通过“水体光学—海色—反照率”的路径改变短波吸收与地表能量收支（海—气耦合）；
- 受海流（P011）平流/扩散控制空间分布；
- 提供观测一致的海色诊断（TrueColor/带反照率、Kd(490)、Chl‑a）。

本版本已实现：
- 多物种（默认 Ns=10，可配置），物种具有独立的“带级光谱形状（Gaussian 权重）与生理参数”；  
- 光竞争（谱带整合）+ 温度响应（Q10）；  
- 单营养盐池 N（Monod，K_N/Yield/再矿化）；  
- 海流平流/扩散；  
- Ocean Color 面板/TrueColor 叠加、Kd(490) 诊断；  
- JSON/NetCDF 导入、导出与 autosave。

---

## 1. 模型变量与方程

### 1.1 状态与单位
- C_i（mg Chl m⁻³）：第 i 种浮游生物混合层体积浓度（i=1..Ns）。  
- N（mmol m⁻³）：单营养盐池（可启用/关闭）。  
- H_mld（m）：混合层深度（默认读取 QD_OCEAN_H_M 或 QD_MLD_M）。  
- 诊断：Chl_total=Σ_i C_i；Kd_b（m⁻¹）；Kd(490)；A_b^water（带）；α_water（标量）。

### 1.2 光场与带平均
- 地表带强度 I_b(surf)：由双星短波离散到 NB 个波段（docs/14）。  
- 混合层带平均：
  Ī_b = I_b(surf) · (1 − exp(−Kd_b H_mld)) / (Kd_b H_mld)

- 带消光（总色素的一参式，性能友好）：  
  Kd_b = Kd0_b + k_chl_b · (Chl_total)^m，默认 m=0.5（QD_PHYTO_KD_EXP_M）。

### 1.3 光谱—生长—竞争
- 每物种定义“带级光谱形状”shape_s[b]（Gaussian 权重，∑_b shape_s[b]=1）。  
- 物种光照代理（谱带整合）：
  E_s = Σ_b (Ī_b · shape_s[b] · Δλ_b)

- 光限制（tanh）：
  muL_s = tanh( α_P · E_s / μ_max_s )

- 温度因子（Q10）：
  f_T = Q10^((T_w − T_ref)/10)

- 营养限制（启用 N 时）：
  f_N,s = N / (K_N,s + N)

- 净增长率（d⁻¹）：
  μ_grow,s = μ_max_s · muL_s · f_T · f_N,s
  μ_s = μ_grow,s − (m0_s + λ_sink/H_mld)

- 物种演化（mg Chl m⁻³ d⁻¹）：
  dC_s/dt = μ_s · C_s

- 营养池演化（mmol m⁻³ d⁻¹）：
  Uptake_s = (μ_grow,s · C_s) / Y_s  
  dN/dt = − Σ_s Uptake_s + R_remin

注：此处 Y_s 单位（mg Chl per mmol N），R_remin 为再矿化速率（常数或后续温度依赖）。

---

## 2. 水体反照率与海色

- 水体带反照率（经验近似，含纯水+浮游生物反射项）：
  A_b^water = A_pure_b + Σ_s [ c_reflect_s · shape_s[b] · (Chl_s)^p_s ]  
  其中 c_reflect_s, p_s 为物种反照率映射系数（默认 c≈0.02, p≈0.5）。

- 标量 α_water（用于能量收支短波）：
  α_water = Σ_b (A_b^water · w_b)，w_b 为短波带权重（Rayleigh/simple）。

- Kd(490)：从 Kd_b 中取最接近 490 nm 的带值作为诊断。

- TrueColor/Ocean Color：用带反照率与动态 per-band 辐照权重合成 RGB（与植被可视化一致，见 scripts.run_simulation 的 plot_true_color）。

---

## 3. 传输与边界

- 平流/扩散：对每物种 C_i 用半拉氏平流与球面拉普拉斯扩散（K_h=QD_PHYTO_KH），与 P011 一致。  
- 极圈一致化：对极圈行做标量平均回填，避免“单物理点的多经度样本”不一致。  
- 湖泊：QD_TREAT_LAKE_AS_WATER=1 时，湖面按水体处理（辐射/蒸发/海色），后续可独立参数化（QD_PHYTO_*_LAKE）。

---

## 4. 环境变量（运行配置）

主控与步长
- QD_PHYTO_ENABLE（默认 1）：开启浮游生物模块  
- QD_PHYTO_NSPECIES（默认 10）：物种数 Ns  
- QD_PHYTO_DT_DAYS（默认 1.0）：日级更新步长  
- QD_PHYTO_ALBEDO_COUPLE（默认 1）：将 α_water 写回短波  
- QD_PHYTO_FEEDBACK_MODE（daily|instant，默认 daily）  
- QD_PHYTO_ADVECTION（默认 1）：启用海流平流/扩散  
- QD_PHYTO_KH（默认=QD_KH_OCEAN 或 5e3 m²/s）

光学带设置
- QD_PHYTO_KD0 / QD_PHYTO_KD_CHL / QD_PHYTO_APURE：每带数组（长度 NB）；或默认  
- QD_PHYTO_KD0_DEFAULT（默认 0.04）  
- QD_PHYTO_KD_CHL_DEFAULT（默认 0.02）  
- QD_PHYTO_APURE_DEFAULT（默认 0.06）  
- QD_PHYTO_KD_EXP_M（默认 0.5）  
- QD_OC_KD_BAND_REF_NM（默认 490）

物种光谱与反照率映射
- QD_PHYTO_SPEC_MU_NM：每物种光谱中心（nm）数组  
- QD_PHYTO_SPEC_SIGMA_NM：每物种宽度（nm）数组（默认 70）  
- QD_PHYTO_SPEC_C_REFLECT：每物种反照率系数（默认 0.02）  
- QD_PHYTO_SPEC_P_REFLECT：每物种幂指数（默认 0.5）

生理参数
- QD_PHYTO_SPEC_MU_MAX：每物种 μ_max（d⁻¹）  
- QD_PHYTO_SPEC_M0：每物种 m0（d⁻¹）  
- 共享默认：QD_PHYTO_ALPHA_P（默认 0.04）、QD_PHYTO_Q10（2.0）、QD_PHYTO_T_REF（293.15 K）、QD_PHYTO_M_LOSS（0.05 d⁻¹）、QD_PHYTO_LAMBDA_SINK（m d⁻¹）

营养池（可选）
- QD_PHYTO_ENABLE_N（默认 0）  
- QD_PHYTO_KN：每物种半饱和 K_N（mmol m⁻³，默认 0.5）  
- QD_PHYTO_YIELD：每物种 Y（mg Chl per mmol N，默认 1.0）  
- QD_PHYTO_REMIN（mmol m⁻³ d⁻¹，默认 0.0）  
- QD_PHYTO_N_INIT（初始 N，默认 1.0）

初始化与持久化
- QD_PHYTO_CHL0（总 Chl 初值，默认 0.05 mg/m³）  
- QD_PHYTO_INIT_FRAC（物种初始权重列表，归一化后生效）  
- QD_PHYTO_INIT_RANDOM（1 随机噪声初始化；0 确定性）  
- autosave：data/phyto_autosave.npz（含 C_i 与 N），data/plankton.json（bio/optics），data/plankton.nc（分布）

---

## 5. 接口与实现（简要）

模块：`pygcm/ecology/phyto.py`
- `class PhytoManager(grid, land_mask, bands=None, H_mld_m=None, diag=True)`  
- `step_daily(insA, insB, T_w, dt_days=1.0) -> (alpha_bands, alpha_scalar)`  
  - 计算 I_b(surf) → Kd_b → Ī_b → E_s → μ_grow,s/μ_s → 更新 C_s；  
  - 若启用 N，按 Uptake/Y 更新 N；  
  - 生成带反照率 A_b^water 与标量 α_water；更新 Kd(490)。  
- `advect_diffuse(uo, vo, dt_seconds)`：物种 C_i 平流/扩散。  
- I/O：`save_bio_json/load_bio_json`、`save_distribution_nc/load_distribution_nc`、`save_autosave/load_autosave`。

驱动脚本：`scripts/run_simulation.py`
- 每个行星日边界调用 `phyto.step_daily(...)`；依据 QD_PHYTO_ALBEDO_COUPLE 更新海面基反照率；  
- 与 P011：每物理步调用 `phyto.advect_diffuse`（若启用）；  
- 可视化：Ocean Color/TrueColor 叠加，Kd(490)/Chl 时间序列打印。

---

## 6. 运行示例

仅诊断（不反馈短波）：
```bash
export QD_PHYTO_ENABLE=1
export QD_PHYTO_NSPECIES=10
export QD_PHYTO_ALBEDO_COUPLE=0
export QD_PLOT_OCEANCOLOR=1
python3 -m scripts.run_simulation
```

回写短波（每日）：
```bash
export QD_PHYTO_ENABLE=1
export QD_PHYTO_NSPECIES=10
export QD_PHYTO_ALBEDO_COUPLE=1
export QD_PHYTO_FEEDBACK_MODE=daily
python3 -m scripts.run_simulation
```

启用营养池竞争 + 传输：
```bash
export QD_PHYTO_ENABLE=1
export QD_PHYTO_ENABLE_N=1
export QD_PHYTO_N_INIT=1.0
export QD_PHYTO_REMIN=0.0
export QD_PHYTO_ADVECTION=1
export QD_PHYTO_KH=5000
python3 -m scripts.run_simulation
```

---

## 7. 验收与数值建议

- 功能性：多物种（Ns=10）可稳定运行；海色与 α_water 随 Chl 增长而增加（绿带增强、蓝带减弱）。  
- 物理合理：高光区 Chl 低、蓝偏；富营养/上升流区 Chl 高、绿偏；温度分带形成纬向差异。  
- 能量守恒：多年平均 |⟨TOA_net⟩|、|⟨SFC_net⟩|、|⟨ATM_net⟩| < 2 W m⁻²（docs/06）。  
- 参数建议：先用 daily 回耦；m≈0.5；μ_max≈1–1.5 d⁻¹；K_N≈0.3–0.8；Y≈1.0；K_h≈5e3 m²/s。  
- 性能：NB=8–16；Ns=10；默认配置下额外成本可控（<10–15% 量级）。

---

## 8. 变更记录（Changelog）

- 2025‑09‑27：实现 Ns=10、谱带光竞争、温度响应与单池营养竞争；I/O 与海流平流/扩散；Ocean Color 可视化与 Kd(490) 诊断；更新运行示例与参数目录。  
- 2025‑09‑25：初稿（个体级思想草案，未实现）。
