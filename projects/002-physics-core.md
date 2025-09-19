# 项目任务：P002 - 核心物理引擎实现

本文档承接 `projects/001-genesis.md` 的规划，旨在将“青黛气候模拟”项目**第一阶段**的设想具体化，形成一份可执行的软件开发计划。

## 1. 目标

本阶段的核心目标是构建 `PyGCM for Qingdai` 软件包的物理引擎基础。这包括三个关键部分：
1.  **轨道力学模块**: 精确计算行星“青黛”随时间变化的能量接收。
2.  **大气动力学核心**: 在球坐标下模拟大尺度的大气环流。
3.  **热力驱动模块**: 将轨道能量输入转化为驱动大气运动的物理力。

## 2. 软件包结构 (`pygcm/`)

为保证代码的模块化与可扩展性，我们将采用以下包结构，所有核心代码位于 `pygcm/` 目录下：

```
qingdai/
├── pygcm/
│   ├── __init__.py
│   ├── constants.py      # 物理与天文常量
│   ├── grid.py           # 球面网格与坐标系统
│   ├── orbital.py        # 轨道力学与能量通量计算
│   ├── forcing.py        # 热力强迫（辐射平衡温度）计算
│   └── dynamics.py       # 大气动力学核心（浅水波模型）
├── scripts/
│   └── run_simulation.py   # 主仿真运行脚本
└── ...
```

## 3. 任务分解与模块API定义

### 任务 3.1: 常量与网格模块 (`constants.py`, `grid.py`)

这两个基础模块需要最先建立。

-   **`constants.py`**:
    -   **目标**: 集中管理所有物理和天文参数。
    -   **内容**:
        -   物理常量 (G, $\sigma$, etc.)。
        -   天文单位 (AU, M_SUN, L_SUN)。
        -   “和光”与“青黛”的核心参数（质量、光度、轨道半长轴、行星半径、自转周期、轴倾角等），数据源为 `docs/01-astronomical-setting.md`。

-   **`grid.py`**:
    -   **目标**: 定义全球经纬度网格。
    -   **API**:
        -   `class SphericalGrid`:
            -   `__init__(self, n_lat, n_lon)`: 初始化一个具有指定分辨率的网格。
            -   `self.lat`, `self.lon`: 存储经纬度坐标的一维数组。
            -   `self.coriolis_param`: 存储科里奥利参数 `f = 2 * Omega * sin(lat)` 的二维数组。

### 任务 3.2: 轨道力学模块 (`orbital.py`)

-   **目标**: 将 `docs/02-orbital-dynamics.md` 中的计算脚本封装成一个可复用的模块。
-   **API**:
    -   `class OrbitalSystem`:
        -   `__init__(self, params)`: 接收一个包含所有天文参数的字典或对象。
        -   `calculate_stellar_positions(self, t)`: 计算并返回 `t` 时刻两颗恒星相对于质心的坐标。
        -   `calculate_total_flux(self, t)`: 计算并返回 `t` 时刻行星接收到的总能量通量 $S_{total}$。这是该模块的核心输出。

### 任务 3.3: 热力驱动模块 (`forcing.py`)

-   **目标**: 计算驱动大气运动的辐射平衡温度场 $T_{eq}$。
-   **依赖**: `orbital.py`, `grid.py`, `constants.py`。
-   **API**:
    -   `class ThermalForcing`:
        -   `__init__(self, grid, orbital_system, planet_params)`: 初始化。
        -   `calculate_insolation(self, t)`: 计算 `t` 时刻每个网格点 `(lat, lon)` 接收的瞬时太阳辐射 $I(\text{lat, lon, t})$。
        -   `calculate_equilibrium_temp(self, t)`: 根据 `docs/03-climate-model.md` 中的公式，计算 `t` 时刻的全球辐射平衡温度场 $T_{eq}(\text{lat, lon, t})$。

### 任务 3.4: 大气动力学核心 (`dynamics.py`)

-   **目标**: 实现基于浅水波方程的GCM核心。
-   **方法**: 采用有限差分法在 `SphericalGrid` 上求解。
-   **API**:
    -   `class ShallowWaterModel`:
        -   `__init__(self, grid, initial_state)`: 初始化模型状态（风场 `u`, `v` 和高度场 `h`）。
        -   `time_step(self, Teq_field, dt)`: 根据当前状态和输入的 $T_{eq}$ 场，使用牛顿冷却方案计算热力强迫，并向前积分一个时间步 `dt`，更新 `u`, `v`, `h`。

## 4. 集成与主仿真循环 (`scripts/run_simulation.py`)

主脚本将协调所有模块，执行完整的仿真流程。

1.  **初始化**:
    -   加载 `constants`。
    -   创建 `SphericalGrid` 实例。
    -   创建 `OrbitalSystem` 实例。
    -   创建 `ThermalForcing` 实例。
    -   创建 `ShallowWaterModel` 实例并设置初始条件（如静止大气）。

2.  **时间积分循环**:
    -   `for t in time_steps:`
        1.  `Teq = forcing_module.calculate_equilibrium_temp(t)`
        2.  `gcm.time_step(Teq, dt)`
        3.  （可选）存储或可视化当前 GCM 状态。

## 5. 下一步行动

1.  **创建包结构**: 在项目根目录下创建 `pygcm/` 目录及上述所有 `.py` 空文件。
2.  **实现 `constants.py`**: 将所有已知参数填入。
3.  **实现 `orbital.py`**: 按照API定义，将现有脚本代码迁移并封装。

完成以上步骤后，我们将拥有一个能够输出动态能量通量的核心模块，为后续GCM的开发奠定基础。
