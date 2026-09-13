# -*- coding: utf-8 -*-
"""
收缩（移动边界）模型的不同写法与敏感性计算

问题四的半径 R(t) 给定，但"收缩如何进入控制方程"存在不同约定，
本模块用同一套有限体积格式分别求解，用于
  1) 论文 8 章的鲁棒性/敏感性分析；
  2) 解释本文算法与 COMSOL 结果差异的来源。

约定：设 s = R(t)/R0（R0 = 2 cm 为初始半径）
  scale_D : 扩散项放大系数（材料坐标下）
  scale_B : 表面传热/传质项放大系数
  code  : scale_D = 1/s^2, scale_B = 1/s   —— 本文采用的写法
  phys  : scale_D = 1,     scale_B = s     —— 按"当前表面积"推导的写法
  mixed : scale_D = 1/s^2, scale_B = s     —— 只把表面项改成物理面积
  none  : scale_D = 1,     scale_B = 1     —— 完全忽略收缩（对照）
"""

import math

import dataio as d

R0 = 0.02          # 初始半径 (m)
H_CONV = 25.0      # 对流换热系数 W/(m^2·K)
H_MASS = 8e-7      # 对流传质系数 m/s
T0 = 301.15        # 初始温度 K
C0 = 2.55          # 初始干基含水率

_pr = d.props(4)
_, _Ta, _Ca = d.air_data()
_, _, _RF = d.radius_curve()


def _interp(dic):
    ts = sorted(dic)

    def f(t):
        if t <= ts[0]:
            return dic[ts[0]]
        if t >= ts[-1]:
            return dic[ts[-1]]
        for k in range(len(ts) - 1):
            if ts[k] <= t <= ts[k + 1]:
                w = (t - ts[k]) / (ts[k + 1] - ts[k])
                return dic[ts[k]] + (dic[ts[k + 1]] - dic[ts[k]]) * w
        return dic[ts[-1]]
    return f


env_T = _interp(_Ta)
env_C = _interp(_Ca)

SCHEMES = {
    "code": (lambda s: 1 / s ** 2, lambda s: 1 / s),
    "phys": (lambda s: 1.0, lambda s: s),
    "mixed": (lambda s: 1 / s ** 2, lambda s: s),
    "none": (lambda s: 1.0, lambda s: 1.0),
}


def _thomas(a, b, c, e, N):
    cp = [0.0] * N
    dp = [0.0] * N
    cp[0] = c[0] / b[0]
    dp[0] = e[0] / b[0]
    for i in range(1, N):
        m = b[i] - a[i] * cp[i - 1]
        cp[i] = c[i] / m if i < N - 1 else 0.0
        dp[i] = (e[i] - a[i] * dp[i - 1]) / m
    x = [0.0] * N
    x[N - 1] = dp[N - 1]
    for i in range(N - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x


def solve(scheme="code", N=200, dt=60.0, t_end=3.0e5, dT=0.0, record_every=600.0):
    """
    求解问题四的收缩模型。
    scheme: code / phys / mixed / none
    dT    : 人为给扩散系数中的温度加一个偏置（敏感性分析用，单位 K）
    返回 dict(t=[s...], Cen=[...], Surf=[...], Tcen=[K...], t_star_h=... 或 None)
    """
    fD, fB = SCHEMES[scheme]
    ds = R0 / N
    rf = [(i + 1) * ds for i in range(N)]
    V = [math.pi * (rf[i] ** 2 - (rf[i - 1] ** 2 if i > 0 else 0.0))
         for i in range(N)]
    T = [T0] * N
    C = [C0] * N
    out = dict(t=[0.0], Cen=[C0], Surf=[C0], Tcen=[T0], t_star_h=None)
    t_star = None
    n = 0
    while n * dt < t_end - 1e-9:
        n += 1
        tt = min(n * dt, t_end)
        s = (_RF(tt) / 100.0) / R0            # 附件 2 的半径单位为 cm
        aD, aB = fD(s), fB(s)
        Ta, Ca = env_T(tt) + 273.15, env_C(tt)
        # ---- 温度 ----
        kk = [_pr["k"](c) * aD for c in C]
        rc = [_pr["rho"](c) * _pr["cp"](c) for c in C]
        kf = [0.5 * (kk[i] + kk[i + 1]) for i in range(N - 1)] + [kk[-1]]
        kb = [0.0] + [0.5 * (kk[i - 1] + kk[i]) for i in range(1, N)]
        A = [2 * math.pi * rf[i - 1] * kb[i] / ds if i > 0 else 0.0
             for i in range(N)]
        B = [2 * math.pi * rf[i] * kf[i] / ds for i in range(N - 1)] + \
            [2 * math.pi * R0 * H_CONV * aB]
        a = [0.0] * N
        b = [0.0] * N
        c = [0.0] * N
        e = [0.0] * N
        for i in range(N):
            b[i] = V[i] * rc[i] / dt + A[i] + B[i]
            e[i] = V[i] * rc[i] / dt * T[i]
            if i > 0:
                a[i] = -A[i]
            if i < N - 1:
                c[i] = -B[i]
            else:
                e[i] += B[i] * Ta
        T = _thomas(a, b, c, e, N)
        # ---- 水分 ----
        aD, aB = fD(s), fB(s)
        DD = [_pr["D"](C[i], T[i] + dT) * aD for i in range(N)]
        Df = [0.5 * (DD[i] + DD[i + 1]) for i in range(N - 1)] + [DD[-1]]
        Db = [0.0] + [0.5 * (DD[i - 1] + DD[i]) for i in range(1, N)]
        A = [2 * math.pi * rf[i - 1] * Db[i] / ds if i > 0 else 0.0
             for i in range(N)]
        B = [2 * math.pi * rf[i] * Df[i] / ds for i in range(N - 1)] + \
            [2 * math.pi * R0 * H_MASS * aB]
        a = [0.0] * N
        b = [0.0] * N
        c = [0.0] * N
        e = [0.0] * N
        for i in range(N):
            b[i] = V[i] / dt + A[i] + B[i]
            e[i] = V[i] / dt * C[i]
            if i > 0:
                a[i] = -A[i]
            if i < N - 1:
                c[i] = -B[i]
            else:
                e[i] += B[i] * Ca
        C = _thomas(a, b, c, e, N)
        if t_star is None and C[0] < 0.15:
            t_star = tt
        if abs(tt % record_every) < 1e-9:
            out["t"].append(tt)
            out["Cen"].append(C[0])
            out["Surf"].append(C[-1])
            out["Tcen"].append(T[0])
        if t_star is not None:
            break
    out["t_star_h"] = t_star / 3600.0 if t_star else None
    return out
