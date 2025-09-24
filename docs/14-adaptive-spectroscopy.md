# 文档 14（v2）：适应性光谱物理学（Adaptive Spectroscopy）

状态（2025‑09‑25）
- [x] v1：理论与双星混合光谱、植物吸收曲线与反射、反馈闭环概念
- [x] v2：光谱带离散方案、接口与数据流、与 06/13/15/015 的耦合规范、环境变量与运行示例
- [ ] M2：参考实现（带划分器、带权函数、I(λ)→I_b、吸收峰生成器、CIE 转换）
- [ ] M3：端到端验证与默认参数固化

交叉引用
- 项目 015：涌现生态学与光谱演化（projects/015-adaptative-spectral-biogeography.md）
- 文档 13：虚拟生命内核（个体吸收/反射带接口，docs/13-plant-model.md）
- 文档 15：生态框架/聚合与反馈（PopulationManager → A_b^surface，docs/15-ecology-framework.md）
- 文档 06：能量收支（短波/反照率耦合 α_total，docs/06-energy-framework.md）
- 文档 02：轨道与入射（双星几何与层顶通量，docs/02-orbital-dynamics.md）

0. 目标与范围
- 将“连续光谱”降维为“有限波段（bands）”，使双星—大气—生态—能量路径可计算、可耦合、可诊断。
- 定义接口：层顶光谱 Irradiance_TOA(λ,t) → 地表光谱 Irradiance_Surface(λ,t) → 波段强度 I_b(t)；个体吸收 A_b → 反射 R_b → 地表带反照率 A_b^surface → 能量模块短波。
- 提供 CIE 可视化路径：Reflected_Light_Spectrum(λ) → RGB。

1. 光谱—波段离散方案

1.1 波段划分
- 范围：可见近可见 λ ∈ [λ_min, λ_max] = [λ0, λ1]（默认 380–780 nm；可扩展至 NIR）
- 波段数 NB：QD_ECO_SPECTRAL_BANDS（默认 8）
- 均匀或自适应划分：
  - 均匀：Δλ_b = (λ1−λ0)/NB
  - 自适应（可选）：按双星混合光谱或大气透过窗口确定可变带宽，保留能源集中的波段更高分辨率
- 带中心 λ_b、带宽 Δλ_b、带权函数 W_b(λ)（矩形/三角或平顶高斯）

1.2 带强度定义
- 带平均地表辐照度：
  I_b(t) = (1/Δλ_b) ∫_{λ∈band b} Irradiance_Surface(λ,t) dλ
- 若使用权函数（非矩形）：
  I_b(t) = [∫ W_b(λ)·Irradiance_Surface(λ,t) dλ] / [∫ W_b(λ) dλ]

2. 层顶 → 地表光谱（GCM 输入端）

2.1 双星混合（文档 02）
- 层顶连续光谱：
  Irr_TOA(λ,t) = L_A_spectrum(λ)/[4π d_A(t)^2] + L_B_spectrum(λ)/[4π d_B(t)^2]
- 或以表观常数 S_A,S_B 与归一化谱形 Φ_A,Φ_B 表示（便于参数化）：
  Irr_TOA(λ,t) = S_A(t)·Φ_A(λ) + S_B(t)·Φ_B(λ)

2.2 大气调制（简化参数化）
- 瑞利散射/分子吸收/云反射/吸收的波长依赖综合成有效透过率 T_atm(λ,t) 与云反照 C_ref(λ,t)：
  Irr_Surf(λ,t) = T_atm(λ,t) · Irr_TOA(λ,t)
- 简化模式（与 QD_ECO_TOA_TO_SURF_MODE 对应）：
  - simple：T_atm(λ) = const; 云仅作为带不相关的减弱因子
  - rayleigh：T_atm(λ) = T0 · (λ/λ_ref)^η（η≈4 的近似衰减指数）
  - custom：可插入表格/经验曲线（按气溶胶/臭氧/云厚等参数）

2.3 带降维
- I_b(t) 由 1.2 的定义从 Irr_Surf(λ,t) 集成得到（可在脚本或 physics/forcing 适当位置实现）

3. 生命吸收与反射（个体层，文档 13 接口）

3.1 基因吸收曲线生成
- 高斯峰叠加模型：
  Absorb(λ) = clip[ Σ_i h_i · exp{− ((λ − c_i)/w_i)^2 }, 0, 1 ]
  - 峰参数（center c_i[nm], width w_i[nm], height h_i≤1）由 genes.spectral_absorption_curve 给定或演化
- 带吸收率 A_b：
  A_b = (1/Δλ_b) ∫_{band b} Absorb(λ) dλ（或 W_b 加权）

3.2 个体能量获取（与文档 13 一致）
- E_gain_i(day) ≈ Σ_b I_b · A_b,i · Δλ_b
- 反射带：
  R_b,i = 1 − A_b,i（限幅 [0,1]）

4. 地表聚合与反馈（文档 15 接口）

4.1 地表带反照率聚合
- PopulationManager 聚合得到：
  A_b^surface = [ Σ_i (R_b,i · LeafArea_i^eff) + Albedo_soil_b · Area_soil ] / Area_total
  - LeafArea_i^eff 可取 = leaf_area_i × QD_ECO_LAI_ALBEDO_WEIGHT
  - Albedo_soil_b 可由 QD_ECO_SOIL_SPECTRUM 提供或采用默认土壤谱

4.2 向能量模块（文档 06）反馈
- 短波路径输入“地表带反照率” A_b^surface；若能量模块仍需要标量 α_total，可由带加权降维：
  α_total = Σ_b A_b^surface · w_b^SW
  - w_b^SW：辐射短波计算使用的带权重（可按太阳常数或当前 I_b 归一化）

4.3 TrueColor / 颜色可视化（可选）
- 反射光谱：
  Reflected(λ) = Irr_Surf(λ) · Albedo_Surface(λ)
- CIE 1931 XYZ → RGB：
  - X = ∫ Reflected(λ)·x̄(λ) dλ，Y = ∫ …ȳ(λ) dλ，Z = ∫ …z̄(λ) dλ
  - 转换到 sRGB 并 gamma 校正；或用带近似表（精度次优）

5. 数值与实现要点

5.1 选择合适的 NB 与带边界
- NB=8–12 对大多数任务足够；若需颜色/物种精细辨识，可升至 16–24
- 建议在 G 与 K 星峰值周围使用更细带（自适应分段），两端（深蓝/近红外）可用宽带

5.2 时间聚合与日尺度
- I_b(t) 可按日平均或日间代表时刻（正午）近似；需与 Ecology 的每日步一致
- 云与气溶胶日内变化可先简化为日均/中位

5.3 单位一致性
- I(λ)：W·m⁻²·nm⁻¹；I_b：W·m⁻²·nm⁻¹（带平均）；E_gain 使用 Σ I_b·A_b·Δλ_b（W·m⁻²）
- 反照率 A_b、R_b：无量纲 [0,1]；注意限幅与数值稳定

6. 接口与数据结构（建议）

6.1 SpectralBands
- fields：
  - nbands: int
  - lambda_edges: np.ndarray[NB+1]
  - lambda_centers: np.ndarray[NB]
  - delta_lambda: np.ndarray[NB]
  - weight_func: Optional[callable] 或枚举 {"rect","tri","gauss"}
- methods：
  - integrate_to_bands(I_lambda: callable or (λ, I) table) → I_b[NB]
  - from_range(nbands, λ0, λ1, mode="uniform"|"adaptive") → SpectralBands

6.2 EcologyAdapter（见项目 015/M1）
- gather_daily_weather → {Ts, Ta, wind, RH, precip, soil_index, I_b[NB], …}
- push_surface_albedo_bands(A_b^surface) → 回写给能量模块或全局状态

6.3 Plant / PopulationManager（见文档 13/15）
- Plant.update_one_day(env) 返回 R_b,i 与 leaf_area
- PopulationManager.aggregate_reflectance(R_b,i, leaf_area_i, soil_bands) → A_b^surface

7. 环境变量（收录至文档 04）

主控
- QD_ECO_SPECTRAL_BANDS（默认 8）
- QD_ECO_SPECTRAL_RANGE_NM（默认 380,780）
- QD_ECO_TOA_TO_SURF_MODE（simple|rayleigh|custom，默认 simple）

大气调制（示例）
- QD_ECO_RAYLEIGH_T0（默认 0.9）
- QD_ECO_RAYLEIGH_LREF_NM（默认 550）
- QD_ECO_RAYLEIGH_ETA（默认 4.0）

土壤与可视化
- QD_ECO_SOIL_SPECTRUM（文件路径；未设则用内置默认）
- QD_ECO_TRUECOLOR_ENABLE（默认 1）
- QD_ECO_TRUECOLOR_GAMMA（默认 2.2）

8. 运行示例（占位）

仅生成带强度与个体/地表带反照率（不回耦能量）
```bash
export QD_ECO_ENABLE=1
export QD_ECO_SPECTRAL_BANDS=8
export QD_ECO_SPECTRAL_RANGE_NM=380,780
export QD_ECO_TOA_TO_SURF_MODE=rayleigh
export QD_ECO_TRUECOLOR_ENABLE=1
python3 -m scripts.run_simulation
```

启用生态与短波耦合（带反照率回写）
```bash
export QD_ECO_ENABLE=1
export QD_ECO_ALBEDO_COUPLE=1
export QD_ECO_SPECTRAL_BANDS=8
python3 -m scripts.run_simulation
```

9. 变更记录（Changelog）
- 2025‑09‑25：v2 规范化：带离散方案、层顶→地表→带降维、个体吸收/反射与聚合、短波耦合路径、CIE 可视化、环境变量与运行示例；与 06/13/15/015 对齐
