import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import matplotlib.pyplot as plt

import pygcm.constants as c

# --- 1. 定义常量与参数 ---
# (Copied from docs/02-orbital-dynamics.md)
G = c.G      # 引力常量 (m^3 kg^-1 s^-2)
M_SUN = c.M_SUN     # 太阳质量 (kg)
L_SUN = c.L_SUN     # 太阳光度 (W)
AU = c.AU        # 天文单位 (m)

# 系统参数 (源自 01-astronomical-setting.md)
M_A = c.M_A    # 主星质量
M_B = c.M_B    # 伴星质量
L_A = c.L_A    # 主星光度
L_B = c.L_B    # 伴星光度
a_bin = c.A_BINARY   # 双星轨道半长轴
a_p = c.A_PLANET     # 青黛行星轨道半长轴

# --- 2. 计算轨道周期 (基于开普勒第三定律) ---
M_total = M_A + M_B
T_bin = 2 * np.pi * np.sqrt(a_bin**3 / (G * M_total))
T_p = 2 * np.pi * np.sqrt(a_p**3 / (G * M_total))

# --- 3. 模拟天体运动轨迹 ---
sim_time = 2 * T_p
dt = 3600 * 24
t_seconds = np.arange(0, sim_time, dt)
t_days = t_seconds / (3600 * 24) # Time in Earth days for plotting

omega_bin = 2 * np.pi / T_bin
omega_p = 2 * np.pi / T_p
r_A = a_bin * (M_B / M_total)
r_B = a_bin * (M_A / M_total)

x_A = r_A * np.cos(omega_bin * t_seconds)
y_A = r_A * np.sin(omega_bin * t_seconds)
x_B = -r_B * np.cos(omega_bin * t_seconds)
y_B = -r_B * np.sin(omega_bin * t_seconds)
x_p = a_p * np.cos(omega_p * t_seconds)
y_p = a_p * np.sin(omega_p * t_seconds)

# --- 4. 计算实时距离与能量通量 ---
d_A = np.sqrt((x_p - x_A)**2 + (y_p - y_A)**2)
d_B = np.sqrt((x_p - x_B)**2 + (y_p - y_B)**2)
S_A = L_A / (4 * np.pi * d_A**2)
S_B = L_B / (4 * np.pi * d_B**2)
S_total = S_A + S_B

# --- 5. 可视化 ---
# Create directories if they don't exist
output_dir = "docs/images"
os.makedirs(output_dir, exist_ok=True)

# Plot 1: Energy Flux vs. Time
plt.style.use('seaborn-v0_8-darkgrid')
fig1, ax1 = plt.subplots(figsize=(12, 6))
ax1.plot(t_days, S_total, label="Total Flux (S_total)", color='crimson')
ax1.set_title('Fig 2.1: Stellar Energy Flux on Qingdai Over Time', fontsize=16)
ax1.set_xlabel('Time (Earth Days)', fontsize=12)
ax1.set_ylabel('Energy Flux (W/m²)', fontsize=12)
ax1.grid(True)

# Add horizontal lines for mean, min, max
mean_flux = np.mean(S_total)
min_flux = np.min(S_total)
max_flux = np.max(S_total)
ax1.axhline(mean_flux, color='gray', linestyle='--', linewidth=1, label=f'Mean: {mean_flux:.2f} W/m²')
ax1.axhline(min_flux, color='blue', linestyle=':', linewidth=1, label=f'Min: {min_flux:.2f} W/m²')
ax1.axhline(max_flux, color='red', linestyle=':', linewidth=1, label=f'Max: {max_flux:.2f} W/m²')
ax1.legend()
fig1.tight_layout()
plot1_path = os.path.join(output_dir, "02-flux-over-time.png")
fig1.savefig(plot1_path)
print(f"Saved plot to {plot1_path}")


# Plot 2: Orbital Paths
fig2, ax2 = plt.subplots(figsize=(10, 10))
# Plot orbits
ax2.plot(x_A / AU, y_A / AU, label='Star A Orbit', color='gold')
ax2.plot(x_B / AU, y_B / AU, label='Star B Orbit', color='orangered')
ax2.plot(x_p / AU, y_p / AU, label='Qingdai Orbit', color='royalblue', linestyle='--')

# Plot initial positions
ax2.plot(x_A[0] / AU, y_A[0] / AU, 'o', color='gold', markersize=10, label='Star A Initial Position')
ax2.plot(x_B[0] / AU, y_B[0] / AU, 'o', color='orangered', markersize=8, label='Star B Initial Position')
ax2.plot(x_p[0] / AU, y_p[0] / AU, 'o', color='royalblue', markersize=6, label='Qingdai Initial Position')

# Plot barycenter
ax2.plot(0, 0, 'k+', markersize=10, label='Barycenter')

ax2.set_title('Fig 2.2: Orbits of Harmony System and Planet Qingdai', fontsize=16)
ax2.set_xlabel('X-coordinate (AU)', fontsize=12)
ax2.set_ylabel('Y-coordinate (AU)', fontsize=12)
ax2.set_aspect('equal', adjustable='box')
ax2.grid(True)
ax2.legend()
fig2.tight_layout()
plot2_path = os.path.join(output_dir, "02-orbital-paths.png")
fig2.savefig(plot2_path)
print(f"Saved plot to {plot2_path}")

plt.close('all')
print("Successfully generated and saved both plots.")
