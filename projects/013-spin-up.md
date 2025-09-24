# Project 013: 快速自旋与平衡初始化（Spin-up & Equilibrium Initialization）

## 0. 状态（2025-09-21）
- [x] 文档定稿（本文件）
- [x] 策略与指标体系设计（两阶段 SOP、能量/水量/气候态稳定性指标）
- [ ] 与运行脚本对齐的参数与重启接口固化（QD_RESTART_IN/OUT、阶段切换）
- [ ] 首次基准 Spin-up 套件运行与验收（相同分辨率/参数组）
- [ ] 平衡态重启文件发布与元数据标准化

关联模块：P006（能量收支）、P007（平板海洋/海冰）、P008（湿度）、P009（水循环）、P010（反噪）、P011（动态海洋）、P012（极区处理）

---

## 1. 背景与问题

- 冷启动（静止、等温）耦合气候模型会因海洋巨大热容导致极长的调整时间（数千年量级），使调试、标定与科学试验成本高昂。
- 当前 GCM 已具备能量收支（P006）、湿度闭环（P008）、水量闭合（P009）、海冰相变（P007）、更好的数值稳定性（P010）以及风驱动浅水海洋（P011）。需要一套统一、可复现的 Spin-up 协议，在不破坏物理一致性的前提下显著缩短至准平衡的时间。

目标：
1) 将达到气候准平衡的时间从“数千年”缩短到“数百年”；  
2) 生成通过诊断验收的“平衡态重启文件（Equilibrium Restart）”；  
3) 为所有后续试验提供一致的“热启动（warm start）”初值。

---

## 2. 加速核心策略

本协议采用“物理合理初始场 + 暂时性低热惯性”的两阶段策略。

### 2.1 构造物理合理初始场（避免等温/静止）
- 为关键场设置分带结构（以纬向为主）而非常数：
  - 大气温度 Ta(φ) 与海表温度 Ts(φ)：建议余弦平方分布，例如
    T(φ) = T_pole + (T_equator − T_pole) · cos²(φ)
  - 初始风：可用弱幅度 Hadley/副热带急流雏形或设为零风但允许快速自旋。
- 目的：一开始即提供合理的南北温度梯度，加速环流形成，降低“从无到有”的过渡耗时。

建议默认：T_equator ≈ 295 K，T_pole ≈ 265 K（可按青黛能量平衡调参）。

### 2.2 暂时性浅混合层（降低海洋热容量）
- 临时将海洋混合层厚度 H_mld 减小（如 50 m → 5 m），等效热容量 C_s,ocean ∝ H_mld 随之下降一个数量级。
- 效果：SST 对辐射与湍流通量的响应加快，热力平衡时间尺度近似按比例缩短。
- 注意：此为“扭曲物理”的加速步骤，仅用于阶段一；阶段二恢复真实 H_mld。

---

## 3. 标准执行协议（SOP）

Spin-up 分两阶段连续执行；阶段二以阶段一的重启文件作为初值。

### 3.1 阶段一：快速热力平衡（Rapid Equilibration）
- 配置：
  - 使用分带初值（Ta/Ts）；
  - 设置浅混合层：QD_MLD_M=5；
  - 开启能量/湿度/水文/海冰/反噪/动态海洋的常规耦合路径（推荐与 P006–P012 的默认一致）；
  - 维持数值稳定：P010 超扩散（QD_FILTER_TYPE=combo，QD_SIGMA4≈0.02）。
- 建议时长：100–200 “青黛年”（视漂移而定）。
- 验收标准（阶段终止条件，计算十年滑动平均趋势）：
  - |d⟨Ts⟩/dt| < 0.05 K/十年；
  - |d⟨TOA_net⟩/dt| < 0.1 W m⁻²/十年；
  - 能量/水量诊断无系统性漂移。
- 产出：restart_phase1.nc（包含所有必要的预报/诊断态，用作阶段二的初值）。

### 3.2 阶段二：恢复真实物理（Full Physics Adjustment）
- 配置：
  - 从阶段一重启（QD_RESTART_IN=restart_phase1.nc）；
  - 恢复 H_mld 至物理值：QD_MLD_M=50（或实验设定真实值）；
  - 其它参数保持一致。
- 建议时长：≥100 “青黛年”。
- 验收标准（更严格的稳态要求）：
  - |⟨TOA_net⟩| ≤ 0.2 W m⁻²（多年平均）；
  - |d⟨Ts⟩/dt| < 0.02 K/十年；
  - 冰面积/体积（若启用）、E–P–R 水量闭合的漂移在阈值内。
- 产出：restart_equilibrium.nc（标准“热启动”文件）。

---

## 4. 诊断与最终验收

- 能量守恒（P006）：TOA/SFC/ATM 净通量多年平均接近 0；LH 与 LH_release 长期平均一致。
- 水量闭合（P009）：长期平均 ⟨E⟩ ≈ ⟨P⟩ + ⟨R⟩；各水库（CWV/海冰/陆面/雪）总量无系统漂移。
- 气候态稳定：全球平均 Ts/海冰面积/洋流动能无长期趋势；环流结构物理合理。
- 数值稳定：未出现“羽状云/条纹”（P010 已启用选择性耗散与必要滤波）；极区数值一致（P012 已启用极点平均回填）。

建议统一打印的诊断条目（每 N 步）：
- [EnergyDiag] TOA_net, SFC_net, ATM_net, ⟨Ts⟩, ⟨Ta⟩；
- [HumidityDiag] ⟨E⟩, ⟨P_cond⟩, ⟨LH⟩, ⟨LH_release⟩；
- [WaterDiag] ⟨E⟩, ⟨P⟩, ⟨R⟩, d/dt[CWV + M_ice/ρ_w + W_land + S_snow]；
- [OceanDiag] KEo, |Uo|max, ⟨SST⟩, CFL；
- Land fraction、Albedo/Friction 统计与地形来源（如使用外部 NetCDF）。

---

## 5. 参数与环境变量（建议默认）

核心控制：
- QD_MLD_M：混合层厚度（阶段一 5；阶段二 50）
- QD_ENERGY_W=1：启用能量收支耦合
- QD_USE_SEAICE=1：启用海冰最小方案（P007）
- QD_USE_OCEAN=1：启用动态海洋（P011，项目近期默认开启）
- QD_OCEAN_USE_QNET=1：SST 吸收净热通量
- QD_FILTER_TYPE=combo、QD_SIGMA4=0.02、QD_SHAPIRO_EVERY=6、QD_SHAPIRO_N=2：反噪稳定（P010）
- QD_ENERGY_DIAG=1、QD_HUMIDITY_DIAG=1、QD_WATER_DIAG=1、QD_OCEAN_DIAG=1：诊断打印
- QD_TS_MIN/QD_TS_MAX（开发期保护，示例 260–305 K）
- 外部地形（可选）：QD_TOPO_NC、QD_USE_TOPO_ALBEDO=1

重启/运行：
- QD_RESTART_IN：阶段二的输入重启文件路径
- QD_RESTART_OUT：当前运行输出重启文件路径
- QD_SIM_DAYS / QD_TOTAL_YEARS：控制总时长（任选其一，视脚本支持；QD_SIM_DAYS 单位：行星日，QD_TOTAL_YEARS 单位：行星年）
- QD_PLOT_EVERY_DAYS：出图间隔

（若脚本尚未实现 QD_RESTART_IN/OUT 或 QD_TOTAL_YEARS，请按本项目落地时补齐，保持与 README 的运行风格一致。）

---

## 6. 任务拆解与里程碑

- T1 初值与脚本接口
  - [ ] 在 scripts/run_simulation.py 中支持：分带初值生成（Ta/Ts 可选）与参数化；
  - [ ] 增加 QD_RESTART_IN/OUT、QD_TOTAL_YEARS 接口；固化保存/加载变量集合。
- T2 阶段一运行（浅混合层）
  - [ ] H_mld=5 m，长期运行；监测 ⟨Ts⟩、⟨TOA_net⟩ 十年滑动趋势；
  - [ ] 产出 restart_phase1.nc。
- T3 阶段二运行（真实混合层）
  - [ ] 从 restart_phase1.nc 重启，H_mld=50 m；
  - [ ] 多年运行与严格标准验收；产出 restart_equilibrium.nc。
- T4 诊断与报告
  - [ ] 统一打印与 NetCDF 输出诊断；给出时间序列图与收敛曲线；
  - [ ] 平衡态文件的元数据完善（参数、网格、版本、时间戳、海陆比分等）。
- T5 默认参数组固化与发布
  - [ ] 标准 spin-up 参数与流程写入 README；
  - [ ] 平衡态重启文件随版本产出，便于复现实验。

---

## 7. 与其他项目的交叉引用

- P006 能量收支：Spin-up 验收的关键（TOA/SFC/ATM 守恒与夜侧温度下限 T_floor 等保护）。
- P007 平板海洋与海冰：H_mld 与海冰反馈直接决定平衡时间尺度与热惯性。
- P008 湿度：LH/LH_release 能量一致性；E 的时空结构与 SST 反馈。
- P009 水循环闭合：E–P–R 与水库闭合的长期稳定，是平衡态的重要验证维度。
- P010 反噪：数值抑噪消除条纹伪迹，避免降水/云的高频噪音污染长期趋势判断。
- P011 动态海洋：风驱动热输送收敛至稳定格局；SST 平流强化热量再分配，加速靠近稳态。
- P012 极区处理：极圈 SST/流速一致化，避免极点数值伪差影响能量与水量诊断。

---

## 8. 运行示例

- 阶段一（Rapid Equilibration；浅混合层 5 m）：
```bash
export QD_TOPO_NC=$(ls -t data/*.nc | head -n1)
export QD_USE_TOPO_ALBEDO=1

export QD_MLD_M=5
export QD_ENERGY_W=1
export QD_USE_SEAICE=1
export QD_USE_OCEAN=1
export QD_OCEAN_USE_QNET=1

export QD_FILTER_TYPE=combo
export QD_SIGMA4=0.02
export QD_SHAPIRO_EVERY=6
export QD_SHAPIRO_N=2

export QD_ENERGY_DIAG=1
export QD_HUMIDITY_DIAG=1
export QD_WATER_DIAG=1
export QD_OCEAN_DIAG=1

export QD_RESTART_OUT="restart_phase1.nc"
# 可选：总时长（若脚本支持年单位）
# export QD_TOTAL_YEARS=150

python3 -m scripts.run_simulation
```

- 阶段二（Full Physics；恢复 50 m）：
```bash
export QD_TOPO_NC=$(ls -t data/*.nc | head -n1)
export QD_USE_TOPO_ALBEDO=1

export QD_MLD_M=50
export QD_ENERGY_W=1
export QD_USE_SEAICE=1
export QD_USE_OCEAN=1
export QD_OCEAN_USE_QNET=1

export QD_FILTER_TYPE=combo
export QD_SIGMA4=0.02
export QD_SHAPIRO_EVERY=6
export QD_SHAPIRO_N=2

export QD_ENERGY_DIAG=1
export QD_HUMIDITY_DIAG=1
export QD_WATER_DIAG=1
export QD_OCEAN_DIAG=1

export QD_RESTART_IN="restart_phase1.nc"
export QD_RESTART_OUT="restart_equilibrium.nc"
# export QD_TOTAL_YEARS=120

python3 -m scripts.run_simulation
```

---

## 9. 输出与元数据建议

为保证可复现性与后续试验一致性，平衡态重启文件建议包含以下元数据：
- 代码版本（git hash）、生成日期 UTC；
- 网格尺寸与行星常数摘要；
- 关键环境变量快照（H_mld、能量/湿度/水文/海冰/反噪/海洋等开关与参数）；
- 外部地形来源（文件名、海陆比、反照率/摩擦统计）；
- 关键诊断（最后 10 年的 ⟨TOA_net⟩、⟨Ts⟩、冰面积等时间平均与标准差）。

---

## 10. 风险与注意事项

- 阶段一扭曲了海洋热容量，仅用于加速逼近平衡态；阶段二必须恢复真实参数并再次验收。
- 若能量或水量诊断长期偏离 0，应优先检查 P006/P008/P009 的耦合路径与参数（如 ε_lw、Bowen、E→LH、P_cond→LH_release、一致的单位换算）。
- 反噪参数过弱可能保留条纹伪迹，过强会过度平滑大尺度；建议使用 P010 的默认“combo”方案并调整 σ₄ 小幅扫描（0.02–0.04）。
- 极区处理（P012）应保持启用，以避免极点数值异常影响诊断与可视化判读。

---

本项目完成后，仓库将提供一个通过验收的标准“平衡态重启文件”，作为全部科学试验的统一初始条件；同时 SOP 与参数接口将统一记录在 README 的“运行 GCM”一节，确保用户以最小的步骤复现热启动结果并进行参数试验或敏感性分析。
