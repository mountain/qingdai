# scripts/verify_calculation.py

import numpy as np
import matplotlib.pyplot as plt

# --- 1. 定义常量与参数 ---
# 物理常量 (SI units)
G = 6.67430e-11      # 引力常量 (m^3 kg^-1 s^-2)
M_SUN = 1.989e30     # 太阳质量 (kg)
L_SUN = 3.828e26     # 太阳光度 (W)
AU = 1.496e11        # 天文单位 (m)

# 参考值
SOLAR_CONSTANT = 1361 # 地球接收的太阳常数 (W/m^2)

# 系统参数 (源自 01-astronomical-setting.md)
M_A = 1.0 * M_SUN    # 主星质量
M_B = 0.8 * M_SUN    # 伴星质量
L_A = 1.0 * L_SUN    # 主星光度
L_B = 0.4 * L_SUN    # 伴星光度
a_bin = 0.5 * AU     # 双星轨道半长轴
a_p = 1.5 * AU       # 青黛行星轨道半长轴

# --- 2. 计算轨道周期 (基于开普勒第三定律) ---
M_total = M_A + M_B

# 双星互绕周期
T_bin = 2 * np.pi * np.sqrt(a_bin**3 / (G * M_total))

# 青黛行星公转周期
T_p = 2 * np.pi * np.sqrt(a_p**3 / (G * M_total))

# 将周期从秒转换为地球日
T_bin_days = T_bin / (3600 * 24)
T_p_days = T_p / (3600 * 24)

print(f"双星互绕周期 (脉冲季): {T_bin_days:.2f} 地球日")
print(f"青黛行星公转年 (年): {T_p_days:.2f} 地球日")

# --- 3. 模拟天体运动轨迹 ---
# 模拟时长为两个行星年，以观察多个周期
sim_time = 2 * T_p
# 时间步长: 1个地球日
dt = 3600 * 24
t = np.arange(0, sim_time, dt)

# 计算角速度
omega_bin = 2 * np.pi / T_bin
omega_p = 2 * np.pi / T_p

# 计算双星各自的轨道半径
r_A = a_bin * (M_B / M_total)
r_B = a_bin * (M_A / M_total)

# 计算每个时间点上各天体的笛卡尔坐标 (质心为原点)
x_A = r_A * np.cos(omega_bin * t)
y_A = r_A * np.sin(omega_bin * t)
x_B = -r_B * np.cos(omega_bin * t)
y_B = -r_B * np.sin(omega_bin * t)
x_p = a_p * np.cos(omega_p * t)
y_p = a_p * np.sin(omega_p * t)

# --- 4. 计算实时距离与能量通量 ---
# 利用勾股定理计算行星到每颗恒星的距离
d_A = np.sqrt((x_p - x_A)**2 + (y_p - y_A)**2)
d_B = np.sqrt((x_p - x_B)**2 + (y_p - y_B)**2)

# 计算来自每颗恒星的能量通量
S_A = L_A / (4 * np.pi * d_A**2)
S_B = L_B / (4 * np.pi * d_B**2)

# 计算总能量通量
S_total = S_A + S_B

# --- 5. 结果分析 (数值) ---
mean_flux = np.mean(S_total)
min_flux = np.min(S_total)
max_flux = np.max(S_total)
flux_variation_percentage = ((max_flux - min_flux) / mean_flux) * 100

print(f"\n能量通量分析:")
print(f"  - 平均值: {mean_flux:.2f} W/m^2")
print(f"  - 最小值: {min_flux:.2f} W/m^2")
print(f"  - 最大值: {max_flux:.2f} W/m^2")
print(f"  - 波动幅度: {flux_variation_percentage:.2f}%")
