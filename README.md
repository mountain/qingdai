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

- 生成水文路由网络（P014，一次性）：
  - 使用 data 下最新地形（包含 elevation/land_mask）：
    - `export QD_TOPO_NC=$(ls -t data/topography_*.nc | head -n1)`
    - `python3 -m scripts.generate_hydrology_maps --topo "$QD_TOPO_NC" --out data/hydrology_network.nc`
  - 未提供外部地形时脚本会回退生成 `land_mask` 并以平坦高程计算 D8（可运行但河网不真实）

- 启用在线径流路由与河网/湖泊叠加（P014 M2–M4）：
  - 最小环境变量：
    - `export QD_HYDRO_ENABLE=1`
    - `export QD_HYDRO_NETCDF=data/hydrology_network.nc`
    - `export QD_HYDRO_DT_HOURS=6`    # 路由步长（小时）
    - `export QD_PLOT_RIVERS=1`       # 在状态图与 TrueColor 叠加河网/湖泊
  - 可选可视化参数：
    - `export QD_RIVER_MIN_KGPS=1e6`  # 河网阈值（kg/s），仅显示主干
    - `export QD_RIVER_ALPHA=0.35`    # 状态图河网透明度（TrueColor 内部为 0.45）
    - `export QD_LAKE_ALPHA=0.40`     # 湖泊透明度
  - 运行：
    - `python3 -m scripts.run_simulation`

- 启用生态模块（P015 M1，小时级回耦）：
  - 最小环境变量（项目当前 M1 约定：NB=16，每物理步子步、即时回耦）
    - `export QD_ECO_ENABLE=1`
    - `export QD_ECO_SUBDAILY_ENABLE=1`
    - `export QD_ECO_SUBSTEP_EVERY_NPHYS=1`      # 每 N 个物理步调用 1 次子步，这里为每步
    - `export QD_ECO_FEEDBACK_MODE=instant`      # 子步立即回写带反照率用于下一物理步
    - `export QD_ECO_ALBEDO_COUPLE=1`            # 开启生态反照率回写
    - `export QD_ECO_SPECTRAL_BANDS=16`          # 光谱带数（建议 16）
    - （可选）TOA→Surface 光谱调制：`export QD_ECO_TOA_TO_SURF_MODE=rayleigh`
  - 运行（与上文相同）：
    - `python3 -m scripts.run_simulation`

参考阅读（含生态进展与参数）：
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
13. 项目状态与进展对齐（2025‑09‑26）：[docs/STATUS-2025-09-26.md](./docs/STATUS-2025-09-26.md)
14. 地形递减率与雪线/雪被（P019）：[docs/18-orography-lapse-and-snowpack.md](./docs/18-orography-lapse-and-snowpack.md)（设计与任务见 [projects/019](./projects/019-orography-lapse-and-snowpack.md)）

## 🤝 贡献

本项目采用人机协作的开发模式。关于协作流程的详细信息，请参阅 [AGENTS.md](./AGENTS.md)。

## 📜 许可证

本项目采用 [MIT 许可证](./LICENSE)。

---

## 👩‍💻 开发者指引（Phase 0）

本节面向开发者，介绍本地开发环境、质量检查与测试流程，以及 P020 Phase 0 的 OO 骨架开关。

### 1) 环境与安装

- 推荐 Python 3.11+（CI 同时在 3.11/3.12/3.13 上验证）
- 克隆代码后，进入仓库根目录，安装开发依赖：
  - 使用 pip：
    ```bash
    python3 -m venv .venv && source .venv/bin/activate
    python3 -m pip install --upgrade pip
    pip install -e '.[dev]'
    ```
  - 或使用 uv：
    ```bash
    uv venv && source .venv/bin/activate
    uv pip install -e '.[dev]'
    ```

说明：
- 项目包本体并非必须安装即可运行工具链（black/ruff/mypy/pytest），但安装 `.[dev]` 可一次性获得这些工具。
- 若遇到 Python 版本约束，可直接使用工具链命令（无需安装包本体）。

### 2) 质量检查与测试（本地）

当前 Phase 0 仅对 `pygcm/world` 与 `tests` 目录执行格式/静态/类型检查（后续阶段将逐步扩大覆盖面）：

```bash
# 代码风格（不修改代码）
ruff check pygcm/world tests

# 格式化（仅检查；如需自动修复可去掉 --check）
black --check pygcm/world tests

# 类型检查
mypy pygcm/world tests

# 单元测试
pytest -q
```

提示：
- `tests/conftest.py` 会将默认网格缩小（10×20）、禁用 autosave 与重型子系统（海洋/生态/路由/海色），并设置 `MPLBACKEND=Agg`，确保测试快速稳定。
- 若需要在单测中启用相关子系统，可在测试用例中使用 `monkeypatch` 单独开启。

### 3) CI

仓库已配置 GitHub Actions（`.github/workflows/ci.yml`）：
- 平台矩阵：Ubuntu/macOS × Python 3.11/3.12/3.13
- 步骤顺序：`ruff check` → `black --check` → `mypy` → `pytest`
- 同样仅针对 `pygcm/world` 与 `tests` 范围做静态与类型检查

### 4) P020 Phase 0：OO 骨架与开关

Phase 0 提供最小 OO 骨架（`pygcm/world`），默认保持 legacy 路径不变：
- 激活 façade（仍运行 legacy 引擎）：
  ```bash
  export QD_USE_OO=1
  python3 -m scripts.run_simulation
  ```
- 仅运行 façade 桩函数（跳过 legacy）：
  ```bash
  export QD_USE_OO=1
  export QD_USE_OO_STRICT=1
  python3 -m scripts.run_simulation
  ```

后续 Phase 1–5 将逐步将模块迁移进 OO 世界对象（参数对象化、纯函数化 forcing/physics、JAX 优先、生态/水文等子系统契约化）。

### 5) 测试策略与范围

- Phase 0：已提供基础用例与 `scripts/test_orbital_module.py` 的 pytest 迁移版
  - `tests/test_phase0_basics.py`
  - `tests/test_orbital_module.py`
- 生态 autosave 的测试迁移将在 Phase 1 引入“契约级测试”（不锁死文件格式细节，以便重构）
- 若需添加更多测试，建议优先选择“物理/数学不变量、契约接口”类型用例，避免对实现细节形成过早约束

### 6) 贡献与协作

- 代码提交前建议先运行本地工具链（ruff/black/mypy/pytest）
- 若提交涉及 CI/工具配置（pyproject/CI/workflows），请说明变更范围与动机
- 参考协作流程：参见 [AGENTS.md](./AGENTS.md)
