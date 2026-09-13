# -*- coding: utf-8 -*-
# 给论文验证章用的全部数字：两法对比 + 守恒残差
import os
import sys

import numpy as np

sys.path.insert(0, r"C:\Users\Wil\Documents\ChatGPT\CUMCM\code2")
import dataio as d

R0 = 0.02
H_M = 8e-7


def cmp_line(name, a, b):
    print("  %-14s 自主 %10.4f   COMSOL %10.4f   相对偏差 %6.3f%%"
          % (name, a, b, 100 * abs(a - b) / abs(b)))


def comsol_p1(field):
    return d.read_comsol_matrix(os.path.join(d.DATA_ROOT, "Problem1",
                                             "comsol_p1_%s.csv" % field))


print("=== 1. 问题一 1800 s")
hdr, T1 = d.result(1)["温度"]
_, C1 = d.result(1)["水分浓度"]
r1T, t1 = comsol_p1("T")
r1C, _ = comsol_p1("C")
j = int(np.argmin(np.abs(np.array(t1) - 1800.0)))
cmp_line("中心温度/℃", T1[-1, 1], r1T[0, 1 + j])
cmp_line("表面温度/℃", T1[-1, -1], r1T[400, 1 + j])
cmp_line("中心含水率", C1[-1, 1], r1C[0, 1 + j])
cmp_line("表面含水率", C1[-1, -1], r1C[400, 1 + j])


print("")
print("=== 2. 问题二 10800 s")
hdr2, T2 = d.result(2)["温度"]
_, C2 = d.result(2)["水分浓度"]
r2T, t2 = d.read_comsol_matrix(os.path.join(d.DATA_ROOT, "Problem2",
                                            "comsol_p2_T.csv"))
r2C, _ = d.read_comsol_matrix(os.path.join(d.DATA_ROOT, "Problem2",
                                           "comsol_p2_C.csv"))
j2 = int(np.argmin(np.abs(np.array(t2) - 10800.0)))
cmp_line("中心温度/℃", T2[-1, 1], r2T[0, 1 + j2])
cmp_line("中心含水率", C2[-1, 1], r2C[0, 1 + j2])
cmp_line("表面含水率", C2[-1, -1], r2C[400, 1 + j2])


print("")
print("=== 3. 烘干时长")
hdr3, D3 = d.result(3)["Sheet1"]
t3, c3 = D3[:, 0], D3[:, 1]
v3, ts3 = d.read_comsol_matrix(os.path.join(d.DATA_ROOT, "Problem3",
                                            "comsol_p3_C.csv"))
k3 = int(np.argmax(c3 < 0.15))
k3r = int(np.argmax(np.round(c3, 4) < 0.15))
k3c = int(np.argmax(v3[0, 1:] < 0.15))
print("  问题三 自主 严格判据      : %8.0f s = %.3f h" % (t3[k3], t3[k3] / 3600))
print("  问题三 自主 四位小数判据  : %8.0f s = %.3f h" % (t3[k3r], t3[k3r] / 3600))
print("  问题三 COMSOL 严格判据    : %8.0f s = %.3f h" % (ts3[k3c], ts3[k3c] / 3600))
print("  相对偏差 %.2f%%" % (100 * abs(t3[k3] - ts3[k3c]) / ts3[k3c]))
hdr4, D4 = d.result(4)["Sheet1"]
t4, c4 = D4[:, 0], D4[:, 1]
tc4, cur4 = d.comsol_table(4, "current")
k4 = int(np.argmax(c4 < 0.15))
k4c = int(np.argmax(cur4[:, 0] < 0.15))
print("  问题四 自主   : %8.0f s = %.3f h" % (t4[k4], t4[k4] / 3600))
print("  问题四 COMSOL : %8.0f s = %.3f h" % (tc4[k4c], tc4[k4c] / 3600))
print("  相对偏差 %.2f%%" % (100 * abs(t4[k4] - tc4[k4c]) / tc4[k4c]))


print("")
print("=== 4. 误差指标（全部时刻 x 21 个位置）")
_, V1T, _ = d.comsol_matrix(1, "T")
_, V1C, _ = d.comsol_matrix(1, "C")
mT = d.metrics(T1[:, 1:], V1T[::20, 1:1801].T)
mC = d.metrics(C1[:, 1:], V1C[::20, 1:1801].T)
print("  问题一 温度   MAE=%.4f ℃  RMSE=%.4f ℃  最大相对=%.3f%%"
      % (mT[0], mT[1], mT[3]))
print("  问题一 含水率 MAE=%.5f   RMSE=%.5f   最大相对=%.3f%%"
      % (mC[0], mC[1], mC[3]))
_, V2T, _ = d.comsol_matrix(2, "T")
_, V2C, _ = d.comsol_matrix(2, "C")
mT2 = d.metrics(T2[:, 1:], V2T[::20, 1:10801].T)
mC2 = d.metrics(C2[:, 1:], V2C[::20, 1:10801].T)
print("  问题二 温度   MAE=%.4f ℃  RMSE=%.4f ℃  最大相对=%.3f%%"
      % (mT2[0], mT2[1], mT2[3]))
print("  问题二 含水率 MAE=%.5f   RMSE=%.5f   最大相对=%.3f%%"
      % (mC2[0], mC2[1], mC2[3]))
_, V3C, _ = d.comsol_matrix(3, "C")
m3 = d.metrics(D3[:3429, 1:], V3C[::20, :3429].T)
print("  问题三 含水率 MAE=%.5f   RMSE=%.5f" % (m3[0], m3[1]))
th4 = t4[:3146] / 3600.0
m4c = d.metrics(c4[:3146], cur4[1:3147, 0])
_, sur4 = d.comsol_table(4, "surface")
m4s = d.metrics(D4[:3146, -1], sur4[1:3147, 0])
print("  问题四 中心 MAE=%.5f  表面 MAE=%.5f" % (m4c[0], m4s[0]))


def balance(path, label):
    _, _, Ca = d.air_data()
    a, times = d.read_comsol_matrix(path)
    r = a[:, 0]
    V = a[:, 1:]
    tarr = np.array(times)
    mid = 0.5 * (r[1:] + r[:-1])
    w = np.zeros(len(r))
    w[0] = 0.5 * mid[0] ** 2
    w[-1] = 0.5 * (r[-1] ** 2 - mid[-1] ** 2)
    w[1:-1] = 0.5 * (mid[1:] ** 2 - mid[:-1] ** 2)
    I = V.T @ w
    tC2 = np.array(sorted(set(Ca)))
    print("  [%s]" % label)
    dt = tarr[1] - tarr[0]
    for tt in (600.0, 1800.0, 3600.0, 10800.0, 36000.0, 72000.0):
        k = int(np.argmin(np.abs(tarr - tt)))
        if k == 0 or k == len(tarr) - 1 or tarr[k] != tt:
            continue
        lhs = (I[k + 1] - I[k - 1]) / (2.0 * dt)
        rhs = -R0 * H_M * (V[-1, k] - np.interp(tt, tC2, [Ca[t] for t in tC2]))
        print("    t=%7.0f s  dI/dt=%12.6e  通量项=%12.6e  残差=%+6.2f%%"
              % (tt, lhs, rhs, 100 * (lhs - rhs) / abs(rhs)))


print("")
print("=== 5. 守恒恒等式检验")
balance(os.path.join(d.DATA_ROOT, "Problem1", "comsol_p1_C.csv"), "问题一")
balance(os.path.join(d.DATA_ROOT, "Problem3", "comsol_p3_C.csv"), "问题三")
