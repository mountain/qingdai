# 青黛世界气候模拟 (PyGCM for Qingdai)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

本项目旨在为原创科幻世界“青黛”开发一个基于物理的简化全球气候模型（GCM）。我们致力于通过科学计算，构建一个具有内在逻辑自洽性的虚构行星气候系统。

## 📖 项目简介

“青黛”是一个环绕“和光”双星系统运行的虚构岩石行星。其独特的天文环境导致了复杂的气候节律。本项目通过 `PyGCM for Qingdai` 这个 Python 软件包，模拟该行星表面的**水、光、热**三大核心生态要素的分布与变化。

更多信息请参阅我们的项目启动文档：
- **[项目启动：青黛行星气候模拟](./projects/001-genesis.md)**

## 📂 目录结构

```
.
├── AGENTS.md               # 项目参与者（人类与AI）的角色定义
├── docs/                   # 核心知识库与技术文档
│   ├── 01-astronomical-setting.md
│   ├── 02-orbital-dynamics.md
│   ├── 03-climate-model.md
│   ├── 04-runtime-config.md   # 运行配置与环境变量目录
│   ├── 05-surface-topography-and-albedo.md
│   ├── 06-energy-framework.md
│   ├── 07-ocean-and-sea-ice.md
│   ├── 08-humidity-and-clouds.md
│   ├── 09-hydrology-closure.md
│   ├── 10-numerics-and-stability.md
│   ├── 11-spin-up-and-restarts.md
│   └── 12-code-architecture-and-apis.md
├── projects/               # 项目高级规划与里程碑
│   └── 001-genesis.md
├── pyproject.toml          # Python 项目配置文件 (待定)
└── README.md               # 本文档
```

## 🚀 快速开始

当前项目已进入原型实现阶段。你可以直接生成地形并导出标准化 NetCDF，再进行可视化检查。

- 安装依赖：
  - `python3 -m ensurepip --upgrade`
  - `python3 -m pip install -r requirements.txt`
- 生成地形与多字段 NetCDF（默认 181x360, seed=42, 目标陆比 0.29）：
  - `python3 -m scripts.generate_topography`
  - 输出目录：`data/`，示例：`topography_qingdai_181x360_seed42_YYYYMMDDTHHMMSSZ.nc`
- 基本可视化（自动选择 data 下最新 nc）：
  - `python3 -m scripts.plot_topography`
  - 将在 `data/` 生成对应的 `*_overview.png`

- 运行 GCM（使用外部地形 NetCDF 与可选地形降水）：
  - 使用 data 下最新 topography：
    - `export QD_TOPO_NC=$(ls -t data/*.nc | head -n1)`
    - `export QD_USE_TOPO_ALBEDO=1`
    - （可选）开启地形降水增强：
      - `export QD_OROG=1`
      - `export QD_OROG_K=7e-4`
    - 运行：
      - `python3 -m scripts.run_simulation`
  - 不使用外部 NetCDF（回退到内置生成）：
    - 不设置 `QD_TOPO_NC`，直接运行：
      - `python3 -m scripts.run_simulation`
  - 其它运行控制（环境变量）：
    - `QD_SIM_DAYS`：模拟时长（单位：行星日，默认 ≈5 个公转周期）
    - `QD_PLOT_EVERY_DAYS`：出图间隔（单位：行星日，默认 10）
    - `QD_DT_SECONDS`：积分步长（秒）
    - 云与降水参数：`QD_CMAX`、`QD_PREF`、`QD_W_MEM`、`QD_W_P`、`QD_W_SRC`
    - 能量框架（P006）：`QD_ENERGY_W`（0..1，能量收支权重）、`QD_ENERGY_DIAG`（能量诊断）、`QD_T_FLOOR`（夜侧温度下限）
    - 湿度–云一致性（P008 M4）：`QD_CLOUD_COUPLE`（启用耦合）、`QD_RH0`、`QD_K_Q`、`QD_K_P`、`QD_PCOND_REF`
    - 水文闭合与径流（P009）：`QD_WATER_DIAG`（水量诊断）、`QD_RUNOFF_TAU_DAYS`（径流时标/天）、`QD_WLAND_CAP`（陆地水库容量/毫米，可选）、`QD_SNOW_THRESH`（雨雪阈值/K）、`QD_SNOW_MELT_RATE`（融雪速率/毫米·天⁻¹）
    - 动力学反噪（P010）：`QD_FILTER_TYPE`（`hyper4|shapiro|spectral|combo`，默认 `combo`）、`QD_SIGMA4`（∇⁴ 自适应强度，默认 0.02）、`QD_DIFF_EVERY`（施加频率，默认 1）、`QD_K4_NSUB`（超扩散子步，默认 1）、`QD_SHAPIRO_N`（默认 2）、`QD_SHAPIRO_EVERY`（默认 6）、`QD_SPEC_EVERY`（谱带阻频率，默认 0=关闭）、`QD_SPEC_CUTOFF`（默认 0.75）、`QD_SPEC_DAMP`（默认 0.5）、`QD_DIFF_FACTOR`（温和全局扩散，默认 0.998）
    - True Color 可视化：`QD_TRUECOLOR_ICE_FRAC`（冰显示阈值，默认 0.15）、`QD_TRUECOLOR_CLOUD_ALPHA`（云不透明度，默认 0.60）、`QD_TRUECOLOR_CLOUD_WHITE`（云白度，默认 0.95）、`QD_TRUECOLOR_SNOW_BY_TS`（是否按温度渲染陆地积雪，默认 0）
    - 说明：脚本启动时会打印地形来源、海陆比例、反照率/摩擦统计等日志，便于检查。


参考阅读：
1.  了解世界观与时间节律：阅读 [docs/01-astronomical-setting.md](./docs/01-astronomical-setting.md)
2.  轨道与气候模型框架：浏览 [docs/02-orbital-dynamics.md](./docs/02-orbital-dynamics.md) 与 [docs/03-climate-model.md](./docs/03-climate-model.md)
3.  运行配置与环境变量目录： [docs/04-runtime-config.md](./docs/04-runtime-config.md)
4.  地形与接入（P004/P005）：[docs/05-surface-topography-and-albedo.md](./docs/05-surface-topography-and-albedo.md)（设计细节参见 [projects/004](./projects/004-topography-generation.md)、[projects/005](./projects/005-topography-integration-into-gcm.md)）
5.  能量收支（P006）：[docs/06-energy-framework.md](./docs/06-energy-framework.md)（方案详见 [projects/006](./projects/006-energy-budget.md)）
6.  海洋与海冰/动态洋流/极点处理（P007/P011/P012）：[docs/07-ocean-and-sea-ice.md](./docs/07-ocean-and-sea-ice.md)（详见 [projects/007](./projects/007-slab-ocean.md)、[projects/011](./projects/011-ocean-model.md)、[projects/012](./projects/012-polar-treatment.md)）
7.  湿度与云–辐射耦合（P003/P008）：[docs/08-humidity-and-clouds.md](./docs/08-humidity-and-clouds.md)（方案详见 [projects/003](./projects/003-cloud-precipitation-albedo.md)、[projects/008](./projects/008-humidity.md)）
8.  水循环闭合（P009）：[docs/09-hydrology-closure.md](./docs/09-hydrology-closure.md)（详见 [projects/009](./projects/009-planetary-hydrology.md)）
9.  数值稳定与反噪（P010）：[docs/10-numerics-and-stability.md](./docs/10-numerics-and-stability.md)（详见 [projects/010](./projects/010-better-dynamics.md)）
10. 快速自旋与重启（P013）：[docs/11-spin-up-and-restarts.md](./docs/11-spin-up-and-restarts.md)（详见 [projects/013](./projects/013-spin-up.md)）
11. 开发者指南/代码架构与 API（P002 + 实现）：[docs/12-code-architecture-and-apis.md](./docs/12-code-architecture-and-apis.md)（参见 [projects/002](./projects/002-physics-core.md)）
12. 地表水文与径流路由（P014）：[projects/014-surface-hydrology.md](./projects/014-surface-hydrology.md)（运行参数见 [docs/04-runtime-config.md](./docs/04-runtime-config.md) 第 10 节）

## 🤝 贡献

本项目采用人机协作的开发模式。关于协作流程的详细信息，请参阅 [AGENTS.md](./AGENTS.md)。

## 📜 许可证

本项目采用 [MIT 许可证](./LICENSE)。
