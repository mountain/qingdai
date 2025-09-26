# 11. 快速自旋与重启（Spin-up & Restarts）

本章整合项目 P013 的策略，提供可复现的两阶段 Spin-up 协议与重启（restart）文件的接口规范，显著缩短耦合气候系统达到准平衡所需的时间，并生成“热启动”初值，服务后续科学试验。

状态（2025-09-21）
- 已完成：两阶段 SOP、验收指标体系与脚本示例（spin-up.sh）
- 待固化：QD_RESTART_IN/OUT、QD_TOTAL_YEARS 在脚本/代码层的统一接口与默认产物发布

关联文档
- docs/06-energy-framework.md（能量收支、诊断指标）
- docs/07-ocean-and-sea-ice.md（混合层深度 H_mld、海冰相变与 α_ice）
- docs/08-humidity-and-clouds.md（E/LH 与 P_cond/LH_release 闭环）
- docs/09-hydrology-closure.md（E–P–R 与水库闭合）
- docs/10-numerics-and-stability.md（数值抑噪与稳定策略）
- projects/013-spin-up.md（原始项目文档）

## 11.1 背景与目标

- 冷启动（静止、等温）下，海洋巨大热容量导致系统需要数千年才能接近平衡。
- 目标：在不破坏物理一致性的前提下，将逼近平衡的时间尺度降至数百年，并产出可复现的“平衡态重启文件（equilibrium restart）”。

## 11.2 加速核心：两阶段策略

阶段 1：快速热力平衡（Rapid Equilibration）
- 物理合理初值：采用分带 Ta(φ)/Ts(φ)，例如 T(φ)=T_pole+(T_eq−T_pole)·cos²φ（典型 T_eq≈295 K、T_pole≈265 K）
- 降低海洋热容量：临时设置浅混合层 H_mld=5 m（海洋响应显著加快）
- 全链路耦合：启用能量/湿度/水文/海冰/动态海洋与数值抑噪的“温和默认”

阶段 2：恢复真实物理（Full Physics）
- 从阶段 1 的重启文件继续；
- 恢复 H_mld 至物理值（默认 50 m）；
- 继续积分到严格稳态标准。

注意：阶段 1“扭曲物理”的仅是热惯性（H_mld），作为加速手段；最终以阶段 2 的稳态为准。

## 11.3 标准执行协议（SOP）

阶段 1（建议 100–200 “青黛年”，或按趋稳指标提前结束）
- 配置要点：H_mld=5；能量诊断与数值抑噪开启；极区/夜侧保护；
- 验收（十年滑动平均趋势）：
  - |d⟨Ts⟩/dt| < 0.05 K/十年
  - |d⟨TOA_net⟩/dt| < 0.1 W m⁻²/十年
  - 能量/水量诊断无系统漂移（见 11.5）

- 产出：restart_phase1.nc（包含必要的预报/诊断态）

阶段 2（建议 ≥100 “青黛年”）
- 配置要点：从 restart_phase1.nc 重启；H_mld=50（或实验既定值）
- 验收（多年平均）：
  - |⟨TOA_net⟩| ≤ 0.2 W m⁻²
  - |d⟨Ts⟩/dt| < 0.02 K/十年
  - 冰面积/体积、水量闭合等无系统漂移

- 产出：restart_equilibrium.nc（标准“热启动”初值）

## 11.4 脚本与运行示例

建议使用仓库提供的 spin-up.sh，也可用环境变量直接运行 scripts/run_simulation.py。

阶段 1（Rapid Equilibration；浅混合层 5 m）
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
# 可选：若脚本支持，以年为单位
# export QD_TOTAL_YEARS=150

python3 -m scripts.run_simulation
```

阶段 2（Full Physics；恢复 50 m）
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

也可直接使用脚本：
```bash
chmod +x spin-up.sh
./spin-up.sh
```

## 11.5 统一诊断与验收指标

能量守恒（docs/06）
- [EnergyDiag] 打印 TOA_net、SFC_net、ATM_net、⟨Ts⟩、⟨Ta⟩
- 稳态多年平均：|⟨TOA_net⟩|、|⟨SFC_net⟩|、|⟨ATM_net⟩| → 0（阈值 ~2 W m⁻²）

水量闭合（docs/09）
- [WaterDiag] 打印 ⟨E⟩、⟨P⟩、⟨R⟩、d/dt[CWV + M_ice/ρ_w + W_land + S_snow]
- 稳态多年平均：⟨E⟩ ≈ ⟨P⟩ + ⟨R⟩，总水库变化 ~ 0

潜热一致（docs/08）
- ⟨LH⟩（SFC） ≈ ⟨LH_release⟩（ATM）

海洋与极区（docs/07）
- OceanDiag：KE、|U_o|max、⟨SST⟩、CFL
- 极圈：启用极点一致化与极区 sponge，SST/流速无经向“条带”

可视化
- Ts/Ta/SST 统一色表；降水图为“行星日累计”；TrueColor 支持云/冰渲染

## 11.6 环境变量与接口

重启与时长
- QD_RESTART_IN：重启输入路径（NetCDF）
- QD_RESTART_OUT：重启输出路径（NetCDF）
- QD_TOTAL_YEARS：按“青黛年”控制时长（优先于 QD_SIM_DAYS）
- QD_SIM_DAYS：按“青黛日”控制时长
- QD_PLOT_EVERY_DAYS：出图间隔
- QD_ORBIT_EPOCH_SECONDS（可选）：仿真起始“天文纪元”t₀（单位：秒）
- QD_ORBIT_EPOCH_DAYS（可选）：仿真起始“天文纪元”t₀（单位：行星日；与 QD_ORBIT_EPOCH_SECONDS 二选一）
  - 纪元优先级：若从 QD_RESTART_IN 或 data/restart_autosave.nc 载入且文件内含 t_seconds，则以文件内的 t_seconds 为准（覆盖上述 EPOCH 变量）。当 QD_AUTOSAVE_LOAD=1 且 autosave 存在时脚本会自动续跑。

初值（分带）
- QD_INIT_BANDED（脚本默认 1；代码默认 0）
- QD_INIT_T_EQ（默认 295 K）、QD_INIT_T_POLE（默认 265 K）

主模块开关与常用参数
- 能量/湿度/水文：见 docs/06、docs/08、docs/09
- 海洋/海冰：见 docs/07
- 数值抑噪与稳定：见 docs/10
- 夜侧温度下限：QD_T_FLOOR（默认 150 K）

## 11.7 输出与元数据建议

为复现与共享，restart 文件应包含：
- 代码版本（git hash）、生成日期（UTC）
- 网格尺寸与行星常数摘要
- 关键环境变量快照（H_mld、能量/湿度/水文/海冰/反噪/海洋等）
- 外部地形来源（文件名、海陆比、反照率/摩擦统计）
- 诊断摘要（最后 N 年的 ⟨TOA_net⟩、⟨Ts⟩、冰面积等均值/标准差）

## 11.8 风险与注意事项

- 阶段 1 的浅混合层仅用于加速；必须进入阶段 2 并通过严标准验收。
- 若 TOA/SFC/ATM 长期偏差，优先检查 docs/06 的参数与云–辐射一致性（α_total、ε_eff 等）。
- 反噪参数过弱会残留条纹伪迹；过强会过度平滑。建议从 docs/10 默认“combo”起步微调。
- 极区一致化（docs/07）建议保持启用，避免极点数值伪差影响守恒诊断与图像判读。

## 11.9 最小操作清单（Checklist）

- [ ] 阶段 1：H_mld=5，自旋到趋稳，输出 restart_phase1.nc
- [ ] 阶段 2：从阶段 1 重启，H_mld=50，自旋到稳态，输出 restart_equilibrium.nc
- [ ] 验收：能量/水量/潜热一致性与极区/可视化整体通过
- [ ] 发布：随 restart_equilibrium.nc 提供元数据与参数快照

## 11.10 变更记录（Changelog）

- 2025‑09‑21：迁移与整合到 docs/11‑spin‑up‑and‑restarts.md；与 06/07/08/09/10 章节交叉引用对齐
