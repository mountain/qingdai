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
│   └── 03-climate-model.md
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

参考阅读：
1.  **了解世界观**: 阅读 [docs/01-astronomical-setting.md](./docs/01-astronomical-setting.md)
2.  **理解物理模型**: 浏览 [docs/02-orbital-dynamics.md](./docs/02-orbital-dynamics.md) 和 [docs/03-climate-model.md](./docs/03-climate-model.md)
3.  **查看项目规划**: 阅读 [projects/001-genesis.md](./projects/001-genesis.md)
4.  **当前项目状态与设计细节（P004）**: [projects/004-topography-generation.md](./projects/004-topography-generation.md)

## 🤝 贡献

本项目采用人机协作的开发模式。关于协作流程的详细信息，请参阅 [AGENTS.md](./AGENTS.md)。

## 📜 许可证

本项目采用 [MIT 许可证](./LICENSE)。
