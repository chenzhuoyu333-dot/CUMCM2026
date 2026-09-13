# -*- coding: utf-8 -*-
"""问题一完整独立程序：定物性预热阶段传热传质有限体积求解。

运行：python code/problem1_complete.py
输出：figure/Q1/*.png，figure/Q1/result1_temperature.csv，result1_moisture.csv
"""
from pathlib import Path
import csv
import math
import numpy as np
import matplotlib.pyplot as plt
from xlsx_output import write_table

# ------------------------------ 模型参数与输出路径
R = 0.02                       # 圆柱半径/m
T_INIT = 28.0 + 273.15         # 初始温度/K
C_INIT = 2.55                  # 初始干基含水率/(kg/kg)
H = 25.0                       # 对流换热系数/(W/(m2 K))
HM = 8.0e-7                    # 对流传质系数/(m/s)
K, RHO, CP = 0.36, 820.0, 2600.0
N, DT, T_END = 200, 1.0, 1800.0
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figure" / "Q1"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
                     "axes.unicode_minus": False, "figure.dpi": 150,
                     "savefig.dpi": 220})

# 根据 essay.tex 表1的升温过程及4 h末值构造分段线性环境边界。
AIR_T = np.array([0, 300, 600, 900, 1200, 1500, 1800], float)
AIR_TEMP = np.array([28.0, 29.4, 31.8, 34.8, 38.1, 40.5, 41.7])
AIR_C = np.full(AIR_T.size, 0.04986)


def ambient(t):
    """附件1边界数据的分段线性插值。"""
    return (float(np.interp(t, AIR_T, AIR_TEMP)),
            float(np.interp(t, AIR_T, AIR_C)))


def diffusion(c):
    """附录2水分扩散系数。"""
    return 7.0e-9 * math.exp(-0.89 / max(float(c), 1.0e-8))


def thomas(lower, diag, upper, rhs):
    """O(N) Thomas追赶法求解三对角方程组。"""
    b, d = diag.copy(), rhs.copy()
    for i in range(1, len(b)):
        factor = lower[i] / b[i - 1]
        b[i] -= factor * upper[i - 1]
        d[i] -= factor * d[i - 1]
    x = np.empty_like(d)
    x[-1] = d[-1] / b[-1]
    for i in range(len(b) - 2, -1, -1):
        x[i] = (d[i] - upper[i] * x[i + 1]) / b[i]
    return x


def implicit_fvm_step(old, capacity, conductivity, boundary_h, boundary_value, dt, dr):
    """柱坐标守恒型有限体积全隐式推进一步。"""
    n = len(old)
    face = np.arange(1, n + 1, dtype=float) * dr
    volume = math.pi * (face**2 - np.r_[0.0, face[:-1]**2])
    coef = np.asarray(conductivity, float)
    if coef.ndim == 0:
        coef = np.full(n, float(coef))
    face_coef = 0.5 * (coef[:-1] + coef[1:])
    internal = 2.0 * math.pi * face[:-1] * face_coef / dr
    west = np.r_[0.0, internal]
    east = np.r_[internal, 2.0 * math.pi * R * boundary_h]
    cap = np.asarray(capacity, float)
    if cap.ndim == 0:
        cap = np.full(n, float(cap))
    lower, upper = -west, -east
    diag = volume * cap / dt + west + east
    rhs = volume * cap / dt * old
    rhs[-1] += east[-1] * boundary_value
    return thomas(lower, diag, upper, rhs)


def solve():
    dr = R / N
    temp = np.full(N, T_INIT)
    moisture = np.full(N, C_INIT)
    times = [0.0]
    temp_hist, moist_hist = [temp.copy()], [moisture.copy()]
    for step in range(1, int(T_END / DT) + 1):
        t = step * DT
        air_temp, air_c = ambient(t)
        temp = implicit_fvm_step(temp, RHO * CP, K, H, air_temp + 273.15, DT, dr)
        diffusivity = np.array([diffusion(c) for c in moisture])
        moisture = implicit_fvm_step(moisture, 1.0, diffusivity, HM, air_c, DT, dr)
        if True:
            times.append(t)
            temp_hist.append(temp.copy())
            moist_hist.append(moisture.copy())
    centers_cm = (np.arange(N) + 0.5) * dr * 100.0
    return np.array(times), centers_cm, np.array(temp_hist) - 273.15, np.array(moist_hist)


def sample(profile, centers_cm, positions):
    return np.interp(positions, centers_cm, profile, left=profile[0], right=profile[-1])


def write_csv(path, times, centers, values):
    positions = np.arange(0.0, 2.01, 0.1)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(["time_s"] + [f"r={x:.1f}cm" for x in positions])
        for t, row in zip(times, values):
            writer.writerow([f"{t:.0f}"] + [f"{x:.6f}" for x in sample(row, centers, positions)])


def save(fig, name):
    fig.tight_layout()
    fig.savefig(OUT / name, bbox_inches="tight")
    plt.close(fig)


def draw_all(times, r, temp, moisture):
    selected = [600, 1200, 1800]
    colors = ["#1d6996", "#38a6a5", "#ed553b"]
    for data, ylabel, filename in [(temp, "温度 (°C)", "q1_T_profile.png"),
                                   (moisture, "含水率 (kg/kg)", "q1_C_profile.png")]:
        fig, ax = plt.subplots(figsize=(5.7, 3.6))
        for t, color in zip(selected, colors):
            k = int(np.argmin(abs(times - t)))
            ax.plot(r, data[k], lw=2.2, color=color, label=f"{t} s")
        ax.set(xlabel="到中心距离 r (cm)", ylabel=ylabel)
        ax.grid(alpha=0.22); ax.legend(frameon=False)
        save(fig, filename)

    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    image = ax.imshow(temp, origin="lower", aspect="auto", cmap="turbo",
                      extent=[0, 2, 0, T_END / 3600])
    ax.set(xlabel="到中心距离 r (cm)", ylabel="时间 (h)")
    fig.colorbar(image, ax=ax, label="温度 (°C)")
    save(fig, "q1_field.png")

    fig, ax = plt.subplots(figsize=(5.7, 3.6))
    ax.plot(times / 60, moisture[:, 0], lw=2.2, label="中心")
    ax.plot(times / 60, moisture[:, -1], lw=2.2, label="表面")
    ax.set(xlabel="时间 (min)", ylabel="含水率 (kg/kg)")
    ax.grid(alpha=0.22); ax.legend(frameon=False)
    save(fig, "q1_center_surface.png")

    fig, ax = plt.subplots(figsize=(6.0, 3.3))
    theta = np.linspace(0, 2 * np.pi, 400)
    ax.fill(np.cos(theta), np.sin(theta), color="#f6d55c", alpha=0.78,
            ec="#8a5a00", lw=2)
    ax.annotate("", xy=(1, 0), xytext=(0, 0),
                arrowprops=dict(arrowstyle="<->", lw=1.6, color="#12355b"))
    ax.text(0.5, 0.07, "R=2 cm", ha="center", color="#12355b")
    ax.text(0, -1.22, "轴心零通量；表面采用第三类对流边界", ha="center")
    ax.set_aspect("equal"); ax.set_xlim(-1.35, 1.35); ax.set_ylim(-1.35, 1.15); ax.axis("off")
    save(fig, "q1_schematic.png")


def main():
    times, r, temp, moisture = solve()
    write_csv(OUT / "result1_temperature.csv", times, r, temp)
    write_csv(OUT / "result1_moisture.csv", times, r, moisture)
    pos = np.arange(0.0, 2.01, 0.1)
    header = ["time_s"] + [f"r={x:.1f}cm" for x in pos]
    rows_t = [[float(tt)] + list(sample(row, r, pos)) for tt, row in zip(times, temp)]
    rows_c = [[float(tt)] + list(sample(row, r, pos)) for tt, row in zip(times, moisture)]
    write_table(ROOT / "result1.xlsx", [("温度", header, rows_t), ("水分浓度", header, rows_c)])
    draw_all(times, r, temp, moisture)
    positions = [0, 0.5, 1.0, 1.5, 2.0]
    print("问题一计算完成；1800 s抽样结果")
    print("r/cm:", positions)
    print("T/°C:", np.round(sample(temp[-1], r, positions), 4))
    print("C:", np.round(sample(moisture[-1], r, positions), 4))
    print("输出目录：", OUT)


if __name__ == "__main__":
    main()
