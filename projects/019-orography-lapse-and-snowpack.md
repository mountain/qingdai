# Project 019：地形递减率 + 雪线 + 雪被水库 + 稳定基流（Orographic Lapse, Snowline, Snowpack and Baseflow）

状态（2025‑09‑27）
- [x] 设计文档（本文件）
- [x] 知识文档（docs/18-orography-lapse-and-snowpack.md）
- [ ] 代码落地（M1–M3）
- [ ] 运行参数目录更新（docs/04）
- [ ] README 引用与示例完善

交叉引用
- 知识/参数化说明：docs/18-orography-lapse-and-snowpack.md
- 地形与 orographic 增强：docs/05、projects/004/005
- 能量/反照率：docs/06
- 湿度与相态：docs/08
- 水量闭合/陆地桶：docs/09
- 路由与湖泊（基流展示）：docs/14、projects/014
- 主循环与代码架构：docs/12、scripts/run_simulation.py

概要
- 在现有 GCM 中引入“海拔温度递减（lapse rate）→ 雪线/相态分配 → 雪被水库（SWE）→ 度日融雪 → 反照率耦合 → 基流/路由”的闭环，使山地/高寒区形成合理的固相蓄水与季节性平滑基流，同时改善积雪反照率反馈与可视化一致性。

---

## 1. 目标与可交付

业务目标
- 用极少参数（Γ、ΔT、DDF/融雪律）补足山地一阶效应：气温随海拔降低、雪线出现、雪季蓄水、春融释水。
- 维持能量/水量长期近守恒，不破坏现有模块稳定性。
- 使河网在春季呈现稳定的“基流高原”，减少尖峰/锯齿。

技术可交付
1) 新的运行时开关与参数（见 §4），在 docs/04 增补条目。
2) 新增状态量与 I/O：
   - SWE_mm（雪当量水层，mm 或 kg/m²）
   - C_snow（积雪覆盖率，0..1；派生）
   - 可选 snow_albedo_age（雪龄诊断，日）
   - 保存/重启到 restart（NetCDF）与 autosave（与现有方案一致）
3) 代码改动（模块分工见 §2）：
   - hydrology：雪被水库与融雪出流、与陆地桶耦合
   - physics/energy：积雪反照率与总 albedo 合成
   - scripts/run_simulation：主循环插入点与诊断/绘图开关
4) 诊断与可视化：
   - SnowDiag（⟨SWE⟩、⟨P_snow⟩、⟨M_snow⟩、⟨C_snow⟩）
   - 雪线等值线/覆盖叠加、TrueColor 雪渲染（按 SWE）

验收标准（概述）
- 多年平均 |⟨TOA/SFC/ATM_net⟩| < 2 W·m⁻²（能量）
- ⟨E⟩ ≈ ⟨P⟩ + ⟨R⟩，d⟨SWE⟩/dt ≈ 0（长期）（水量）
- 山地/中高纬春季出现稳定基流；雪盖与 TrueColor/反照率一致

---

## 2. 设计与集成（代码层）

2.1 Lapse Rate（海拔温度递减）
- 位置：scripts/run_simulation（或 pygcm/physics 统一函数）
- 输入：T_a、T_s、elevation、Γ（K/km）、Γ_s（可=Γ）
- 输出：T_hat_a = T_a − Γ·ΔH_km；T_hat_s 类似（用于雪反照率/可视化）
- 保护：与 QD_T_FLOOR 一致的温度地板

2.2 相态分配（雪/雨混合）
- 位置：pygcm/hydrology.py（与现有 partition_precip_phase 接口对齐/新增 v2）
- 新接口：partition_precip_phase_smooth(P_total, T_hat_a, T_thresh, ΔT) → (P_rain, P_snow, f_snow)
- 说明：用 Sigmoid 平滑过渡，避免硬阈值棋盘伪迹

2.3 雪被水库（SWE）与融雪
- 位置：pygcm/hydrology.py
- 新增函数：
  - snowpack_step(SWE, P_snow, T_hat_a, dt, params) → (SWE_next, M_snow, C_snow, snow_albedo)
    - 融雪模式：degree_day（推荐）或 constant（回退到 QD_SNOW_MELT_RATE）
    - C_snow = 1 − exp(−SWE/SWE_ref)
    - snow_albedo：可选雪龄衰减（α_fresh→α_old，τ_alb）
- 桶与径流耦合：
  - W_land ← W_land + (1−φ_fast)·M_snow·dt
  - 快流 φ_fast·M_snow·dt → 路由输入缓存（可选，默认 0）

2.4 反照率耦合与能量
- 位置：pygcm/physics.py 或 energy.py 中组装层
- 修改 calculate_dynamic_albedo：
  - 对陆地：α_surface_eff = α_base·(1−C_snow) + α_snow·C_snow（限幅）
  - 与云/海冰合成总 α_total（原逻辑不变）
- TrueColor：支持按 SWE 渲染（QD_TRUECOLOR_SNOW_BY_SWE）

2.5 主循环时序（简述）
1) 动力/湿度/云 → 得到 P_total、T_a、T_s
2) 计算 T_hat_a/(T_hat_s)
3) 相态分配（P_rain/P_snow）
4) 雪被步：SWE/M_snow/C_snow → 写回 W_land 与快流缓存
5) 反照率合成（含 C_snow）
6) 能量步 → 海洋步
7) 水文步（桶/径流）→ 路由（到水文步长）
8) 诊断/绘图（SnowDiag/雪线/TrueColor）

2.6 状态/重启
- restart NetCDF（推荐字段）：
  - SWE_mm（float32, land-only；海洋/湖泊为 0）
  - snow_albedo_age_days（可选）
- autosave：与现有 data/restart_autosave.nc 或 eco_autosave 分离；按主脚本策略写入/加载

---

## 3. 数据结构与单位

- P_total、P_rain、P_snow、M_snow：kg m⁻² s⁻¹（可打印 mm/day）
- SWE_mm：mm（内部等价 kg m⁻²，1 mm = 1 kg m⁻²）
- C_snow：0..1（派生，不必持久化亦可）
- elevation：m；ΔH_km = (H − H_ref)/1000
- 反照率 α：无量纲 0..1

数值保护
- SWE_next = max(0, SWE + P_snow·dt − M_snow·dt)
- M_snow ≤ SWE/dt（步内不透支）
- α_total ∈ [0,1]；C_snow ∈ [0,1]

---

## 4. 运行参数（新增/沿用，最终以 docs/04 为准）

主控
- QD_LAPSE_ENABLE（1）/ QD_LAPSE_K_KPM（6.5）/ QD_LAPSE_KS_KPM（=K_KPM）
- QD_SWE_ENABLE（1）/ QD_SWE_INIT_MM（0）/ QD_SWE_MAX_MM（可选 None）

相态
- QD_SNOW_THRESH（273.15 K）
- QD_SNOW_T_BAND（1.5 K）

融雪/雪反照率
- QD_SNOW_MELT_MODE（degree_day|constant，默认 degree_day）
- QD_SNOW_DDF_MM_PER_K_DAY（3.0）
- QD_SNOW_MELT_RATE（5.0 mm/day；常数后备）
- QD_SNOW_MELT_TREF（273.15 K）
- QD_SWE_REF_MM（15）
- QD_SNOW_ALBEDO_FRESH（0.70）
- QD_SNOW_ALBEDO_OLD（0.45）
- QD_SNOW_ALBEDO_DECAY_DAYS（10）
- QD_SNOW_FASTFLOW_FRAC（0.0）

可视化/诊断
- QD_PLOT_SNOWLINE（1）
- QD_TRUECOLOR_SNOW_BY_SWE（1）
- 打印间隔（沿用 QD_PLOT_EVERY_DAYS 或独立 SnowDiag_EVERY）

---

## 5. 任务拆解与里程碑

M0 准备（已完成）
- [x] 需求/知识文档：docs/18
- [x] 本项目设计文档

M1 核心实现（lapse + 相态 + SWE + 融雪）
- [ ] hydrology.partition_precip_phase_smooth（Sigmoid）
- [ ] hydrology.snowpack_step（SWE/M_snow/C_snow/可选雪龄）
- [ ] run_simulation：插入 lapse/相态/SWE 调用，写回 W_land 与快流缓存
- [ ] 状态定义：SWE 持久化；加载/保存路径

M2 反照率与可视化
- [ ] physics/energy：在陆地分支融合 C_snow 到 α_surface_eff
- [ ] TrueColor：按 SWE 渲染开关（QD_TRUECOLOR_SNOW_BY_SWE）
- [ ] 雪线/覆盖叠加到状态图

M3 诊断与验证
- [ ] SnowDiag 日志：⟨SWE⟩、⟨P_snow⟩、⟨M_snow⟩、⟨C_snow⟩
- [ ] 与能量/水量闭合联测（docs/06/09）
- [ ] 河网/湖泊叠加验证“春融基流”（docs/14）

M4 文档与运行指引
- [ ] docs/04 增补变量目录
- [ ] README：在“运行 GCM”添加本功能示例
- [ ] 变更记录与默认参数建议

---

## 6. 测试计划

类型与范围
- 单元：
  - 相态 Sigmoid：ΔT 扫描，硬阈值对照
  - SWE/M_snow 守恒：典型温度轨迹（日/年）与限幅
- 集成：
  - 山地带实验：开启/关闭 lapse 与 orographic，比较 P_snow、SWE、R/flow_accum
  - 季节场景（年际）：⟨SWE⟩ 周期；春季 M_snow 峰；基流高原
- 回归：
  - 能量/水量闭合阈值
  - 性能：启用后步时开销可控（<3–5% 额外）

可视化检查
- 雪线随季节迁移（纬带与山地）
- TrueColor 与 α_total 中的雪一致
- 河网主干在雪融期增强，路由 mass_closure 误差小

---

## 7. 风险与对策

- 反馈放大（高反照率→更冷→滞融）：使用 degree_day（温和）、α_snow 衰减、SWE_ref 调节
- 网格锯齿（硬阈值）：已用 Sigmoid 平滑，ΔT 可调
- 双计能量：仅在陆地分支融合雪反照率，海冰仍走 P007
- 性能：新增计算为 O(N) 线性，避免高频昂贵卷积；诊断频率与绘图可下采样
- 重启一致性：SWE 与能量/水量诊断时钟同步；保存/加载在 run_simulation 的现有时机执行

---

## 8. 运行示例（精简）

基础（按 docs/18）
```bash
export QD_LAPSE_ENABLE=1
export QD_LAPSE_K_KPM=6.5
export QD_SWE_ENABLE=1
export QD_SNOW_T_BAND=1.5
export QD_SNOW_MELT_MODE=degree_day
export QD_SNOW_DDF_MM_PER_K_DAY=3.0
export QD_SWE_REF_MM=15
export QD_TRUECOLOR_SNOW_BY_SWE=1
python3 -m scripts.run_simulation
```

叠加地形降水增强
```bash
export QD_OROG=1
export QD_OROG_K=7e-4
python3 -m scripts.run_simulation
```

---

## 9. 变更记录（Changelog）

- 2025‑09‑27：v1 提案（设计与知识文档完成；待代码落地与参数目录更新）
