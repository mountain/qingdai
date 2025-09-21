# 9. 行星水循环闭合（Hydrology Closure: E–P–R–Reservoirs）

本章系统化迁移项目 P009 的方案，建立最小但自洽的“行星水循环闭合”框架：蒸发（E）→ 大气水汽 → 凝结/降水（P）→ 地表水库（海/陆/冰/雪）→ 径流回海（R），并与能量（P006）、海洋/海冰（P007/011）、湿度（P008）保持能量—水量一致。

状态（2025-09-21）
- 已实现：水库定义（大气水汽/海冰/陆地贮水/雪）、陆地“桶”模型与径流 R、相态分配（雨/雪）、全局闭合诊断
- 待完善：长期标定（τ_runoff、雪阈值/融化率）、图像输出套件

关联文档
- docs/06-energy-framework.md：潜热能量一致性（LH/LH_release）
- docs/07-ocean-and-sea-ice.md：海冰相变与 α_ice、SST 耦合
- docs/08-humidity-and-clouds.md：E 与 P_cond 的计算与耦合
- docs/04-runtime-config.md：运行时环境变量总表（含水文条目）
- projects/009-planetary-hydrology.md：原始实现与路线图


## 9.1 目标与原则

- 定义最小水库集合并演化：大气水汽（CWV）、海冰质量（M_ice）、陆地贮水（W_land）、积雪（S_snow）
- 建立陆地“桶”模型和径流 R 以闭合 E–P–R 收支，并将 R 记账回海洋
- 在多年平均/稳定期实现质量守恒：⟨E⟩ ≈ ⟨P⟩ + ⟨R⟩；能量潜热一致：⟨LH⟩ ≈ ⟨LH_release⟩


## 9.2 水库与连续方程（面积权重 w=cosφ）

大气水汽（CWV）
- 近地层单层近似：CWV ~ ρ_a h_mbl q
- 全星积分趋势：d⟨CWV⟩/dt ≈ ⟨E⟩ − ⟨P⟩（平流散度的全球积分趋近 0）

海冰（M_ice）
- 相变由能量（Q_net）主导，质量通量 d⟨M_ice⟩/dt = ρ_i ⟨dh_ice/dt⟩（冻结正、融化负）
- 对水量而言：冻结将海水转入冰库；融化将冰回归海水

陆地贮水（W_land）
- dW_land/dt = P_land − E_land − R
- 线性径流 R = W_land / τ_runoff；可选容量上限 W_cap 超限溢流

积雪（S_snow）
- P_snow 在 T_s < T_thresh 时累积；M_snow 以固定/温度调制速率融化
- 融雪进入 W_land 或直接 R（简化）

全局守恒（期望长期平均）
- d/dt [⟨CWV⟩ + ⟨M_ice/ρ_w⟩ + ⟨W_land⟩ + ⟨S_snow⟩] ≈ ⟨E⟩ − ⟨P⟩ − ⟨R⟩ → 0  
- 稳态：⟨E⟩ ≈ ⟨P⟩ + ⟨R⟩


## 9.3 相态分配：雨/雪

- 阈值温度 T_thresh（默认 273.15 K）
- 若 T_s < T_thresh：P → P_snow；否则 P → P_rain
- P_snow 累入 S_snow，P_rain 入 W_land（或直接径流一部分，见桶模型）


## 9.4 陆地“桶”模型与径流 R

- 状态：W_land（单位 mm 或 kg m⁻²，内部有单位换算）
- 源汇：+P_land、−E_land、−R
- 排水：R = W_land / τ_runoff（线性），可扩展“双时标排水”
- 容量（可选）：W_cap，超出部分立即作为快流排出
- 记账：R 全部计入“回海”闭合，不创建独立海洋淡水库（最小化）


## 9.5 与能量/湿度/海冰的一致性

- 能量潜热一致（见 docs/06 与 docs/08）：  
  LH = L_v·E（地表负项），LH_release = L_v·P_cond（大气正项），长期 ⟨LH⟩ ≈ ⟨LH_release⟩
- 海冰相变（见 docs/07）：冻结/融化的质量与能量一致；冰面蒸发显著减弱；α_ice 提高反照率
- 水量闭合：E、P（含 P_cond 与动力诊断融合）与 R 在 WaterDiag 中统一核算


## 9.6 接口与实现位置（pygcm/hydrology.py，scripts/run_simulation.py）

关键函数（建议/现状）
- partition_precip_phase(P, T_s, T_thresh) → (P_rain, P_snow)
- snow_step(S_snow, P_snow, melt_rate, dt) → (S_next, M_snow)
- update_land_bucket(W_land, P_land, E_land, dt, τ_runoff, W_cap=None) → (W_next, R)
- diagnose_water_closure(q, E, P, R, h_ice, W_land, S_snow, grid) → dict  
  返回 ⟨E⟩、⟨P⟩、⟨R⟩、d/dt[各水库] 与闭合误差

主循环接入顺序（建议）
1) 云/降水（动力+湿度一致性）得到 P  
2) 相态分配（雨/雪）与雪库更新  
3) 陆地桶/径流更新，R 记账回海  
4) 诊断输出（WaterDiag）  
5) 与能量/湿度诊断共同检查长期守恒


## 9.7 环境变量（节选；详见 docs/04-runtime-config.md）

- QD_RUNOFF_TAU_DAYS（默认 10）：陆地径流时标（天）
- QD_WLAND_CAP（可选，单位 mm）：桶容量上限
- QD_SNOW_THRESH（默认 273.15 K）：雨/雪阈值
- QD_SNOW_MELT_RATE（默认 5 mm/day）：融雪速率
- QD_RHO_W（默认 1000）：必要时的质量↔厚度换算
- QD_WATER_DIAG（默认 1）：打印水量闭合诊断


## 9.8 推荐默认与运行示例

```bash
export QD_ENERGY_W=1
export QD_ENERGY_DIAG=1
export QD_HUMIDITY_DIAG=1
export QD_WATER_DIAG=1

# 水文参数（温和默认）
export QD_RUNOFF_TAU_DAYS=10
export QD_SNOW_THRESH=273.15
export QD_SNOW_MELT_RATE=5
# export QD_WLAND_CAP=50  # 可选

python3 -m scripts.run_simulation
```


## 9.9 诊断与验收标准（建议）

- 守恒：长期平均 ⟨E⟩ ≈ ⟨P⟩ + ⟨R⟩；d/dt[CWV + M_ice/ρ_w + W_land + S_snow] ≈ 0  
- 能量一致：长期平均 ⟨LH⟩ ≈ ⟨LH_release⟩  
- 物理合理：冷区/高地雪季合理；融雪季 R 增强；海陆差异清晰（海上 E 大、冰上 E ≪）  
- 输出：WaterDiag 打印收支与水库变化；图像/时间序列可选保存


## 9.10 变更记录（Changelog）

- 2025‑09‑20：水库/相态/桶模型/闭合诊断建立，联通能量与湿度路径  
- 2025‑09‑21：与海冰/海洋耦合细化；诊断打印完善  
- 2025‑09‑21：文档迁移与整合至 docs/09‑hydrology‑closure.md；与 04/06/07/08 交叉引用对齐
