# -*- coding: utf-8 -*-
# 用 COMSOL 有限元仿真的结果独立验证本算法（有限体积法）
# 本算法先独立完成建模与求解，COMSOL 承担第二种数值方法交叉验证的角色
# 运行：python validate_comsol.py
import os
import re
import sys
import zipfile

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plotting import save

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = r"C:\Users\Wil\Desktop\cumcm_comsol"


def read_result(path):
    """读本程序生成的 result*.xlsx"""
    z = zipfile.ZipFile(path)
    wb = z.read("xl/workbook.xml").decode("utf-8")
    names = re.findall(r'<sheet name="([^"]+)"', wb)
    files = sorted([n for n in z.namelist() if re.match(r"xl/worksheets/sheet\d+\.xml$", n)],
                   key=lambda s: int(re.search(r"(\d+)", s.split("/")[-1]).group(1)))
    out = {}
    for nm, fn in zip(names, files):
        xml = z.read(fn).decode("utf-8")
        rows = re.findall(r"<row[^>]*>(.*?)</row>", xml, re.S)
        hdr = re.findall(r"<t>(.*?)</t>", rows[0], re.S)
        data = np.array([[float(v) for v in re.findall(r"<v>(.*?)</v>", r, re.S)]
                         for r in rows[1:]])
        out[nm] = (hdr, data)
    return out


def read_comsol(path):
    """COMSOL 数据导出：行 = 位置，列 = 时刻"""
    hdr, rows = None, []
    for ln in open(path, encoding="utf-8", errors="replace"):
        if ln.startswith("%"):
            if "@ t=" in ln:
                hdr = ln
            continue
        if ln.strip():
            rows.append([float(x) for x in ln.strip().split(",")])
    ts = [float(m.group(1)) for m in re.finditer(r"@\s*t\s*=\s*([-+0-9.eE]+)", hdr or "")]
    return np.array(rows), ts


def read_comsol_table(path):
    """COMSOL 点计算表导出：行 = 时刻，列 = 各点"""
    data, hdr = [], None
    for ln in open(path, encoding="utf-8", errors="replace"):
        if ln.startswith("%"):
            if ln.startswith("% Time"):
                hdr = ln
            continue
        if ln.strip():
            data.append([float(x) for x in ln.strip().split(",")])
    cols = [float(m.group(1)) for m in re.finditer(r"点[:：]\s*([0-9.]+)", hdr or "")]
    return np.array(data), cols


def metrics(a, b):
    d = np.abs(a - b)
    rel = d / np.maximum(np.abs(b), 1e-12)
    return d.mean(), np.sqrt((d ** 2).mean()), d.max(), rel.max() * 100


def report(tag, ours, comsol, unit):
    mae, rmse, mx, rel = metrics(ours, comsol)
    print("  %-6s  MAE=%.4f  RMSE=%.4f  最大偏差=%.4f %s  最大相对偏差=%.2f%%"
          % (tag, mae, rmse, mx, unit, rel))
    return rel


def main():
    print("=" * 78)
    print("COMSOL 有限元仿真 vs 本文有限体积算法（交叉验证）")
    print("=" * 78)
    worst = 0.0
    # ---------------- 问题 1 ----------------
    f1 = read_result(os.path.join(HERE, "result1.xlsx"))
    cT, tsT = read_comsol(os.path.join(CSV, r"Problem1\comsol_p1_T.csv"))
    cC, tsC = read_comsol(os.path.join(CSV, r"Problem1\comsol_p1_C.csv"))
    hdr, ours_T = f1["温度"]
    ours_C = f1["水分浓度"][1]
    t_ours = ours_T[:, 0]
    common = [t for t in t_ours if t in tsT and t in tsC]
    jO = [int(np.where(t_ours == t)[0][0]) for t in common]
    jC = [tsT.index(t) for t in common]
    A_T = np.array([ours_T[k, 1:] for k in jO])
    B_T = np.array([[cT[r][1 + j] for r in range(0, 401, 20)] for j in jC])
    A_C = np.array([ours_C[k, 1:] for k in jO])
    B_C = np.array([[cC[r][1 + j] for r in range(0, 401, 20)] for j in jC])
    print("问题 1（1800 s，21 个位置 x %d 个时刻）" % len(common))
    worst = max(worst, report("温度", A_T, B_T, "C"), report("含水率", A_C, B_C, "kg/kg"))
    m_late = np.array([t >= 600.0 for t in common])
    print("  排除最初 10 min 后：")
    report("温度", A_T[m_late], B_T[m_late], "C")
    report("含水率", A_C[m_late], B_C[m_late], "kg/kg")
    r = np.array([float(x) for x in hdr[1:]])
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    axs[0].plot(r, A_T[-1], "o-", ms=3, label="本文算法")
    axs[0].plot(r, B_T[-1], "s--", ms=3, label="COMSOL")
    axs[0].set_xlabel("到药材中心的距离 r / cm")
    axs[0].set_ylabel("温度 T / C")
    axs[0].set_title("t = 1800 s 温度剖面")
    axs[1].plot(r, A_C[-1], "o-", ms=3, label="本文算法")
    axs[1].plot(r, B_C[-1], "s--", ms=3, label="COMSOL")
    axs[1].set_xlabel("到药材中心的距离 r / cm")
    axs[1].set_ylabel("水分浓度 C / (kg/kg)")
    axs[1].set_title("t = 1800 s 含水率剖面")
    axs[0].legend()
    axs[1].legend()
    save(fig, "problem1", "fig4_validation.png", "问题1：两种数值方法对比")
    # ---------------- 问题 2 ----------------
    f2 = read_result(os.path.join(HERE, "result2.xlsx"))
    cT2, tsT2 = read_comsol(os.path.join(CSV, r"Problem2\comsol_p2_T.csv"))
    cC2, tsC2 = read_comsol(os.path.join(CSV, r"Problem2\comsol_p2_C.csv"))
    ours_T2 = f2["温度"][1]
    ours_C2 = f2["水分浓度"][1]
    t2 = ours_T2[:, 0]
    common2 = [t for t in t2 if t in tsT2 and t in tsC2]
    jO2 = [int(np.where(t2 == t)[0][0]) for t in common2]
    jC2 = [tsT2.index(t) for t in common2]
    A_T2 = np.array([ours_T2[k, 1:] for k in jO2])
    B_T2 = np.array([[cT2[rr][1 + j] for rr in range(0, 401, 20)] for j in jC2])
    A_C2 = np.array([ours_C2[k, 1:] for k in jO2])
    B_C2 = np.array([[cC2[rr][1 + j] for rr in range(0, 401, 20)] for j in jC2])
    print("问题 2（3 h，21 个位置 x %d 个时刻）" % len(common2))
    worst = max(worst, report("温度", A_T2, B_T2, "C"), report("含水率", A_C2, B_C2, "kg/kg"))
    m_late2 = np.array([t >= 600.0 for t in common2])
    print("  排除最初 10 min 后：")
    report("温度", A_T2[m_late2], B_T2[m_late2], "C")
    report("含水率", A_C2[m_late2], B_C2[m_late2], "kg/kg")
    # ---------------- 问题 3 ----------------
    f3 = read_result(os.path.join(HERE, "result3.xlsx"))
    cC3, ts3 = read_comsol(os.path.join(CSV, r"Problem3\comsol_p3_C.csv"))
    hdr3, ours3 = f3["Sheet1"]
    t3 = ours3[:, 0]
    common3 = [t for t in t3 if t in ts3]
    jO3 = [int(np.where(t3 == t)[0][0]) for t in common3]
    jC3 = [ts3.index(t) for t in common3]
    A3 = np.array([ours3[k, 1:] for k in jO3])
    B3 = np.array([[cC3[rr][1 + j] for rr in range(0, 401, 20)] for j in jC3])
    print("问题 3（到 57 h，21 个位置 x %d 个时刻）" % len(common3))
    worst = max(worst, report("含水率", A3, B3, "kg/kg"))
    m_late3 = np.array([t >= 3600.0 for t in common3])
    print("  排除最初 1 h 后：")
    report("含水率", A3[m_late3], B3[m_late3], "kg/kg")
    td_ours = ours3[int(np.argmax(ours3[:, 1] < 0.15)), 0]
    td_comsol = None
    for j, t in enumerate(ts3):
        if cC3[0][1 + j] < 0.15:
            td_comsol = t
            break
    print("  烘干时长：本文算法 %.2f h ； COMSOL %.2f h（相对偏差 %.2f%%）"
          % (td_ours / 3600, td_comsol / 3600,
             100 * abs(td_ours - td_comsol) / td_comsol))
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.plot(t3 / 3600, ours3[:, 1], "-", label="本文算法（中心）")
    ax.plot(t3 / 3600, ours3[:, -1], "-", label="本文算法（表面）")
    ax.plot([t / 3600 for t in ts3], [cC3[0][1 + j] for j in range(len(ts3))],
            "--", label="COMSOL（中心）")
    ax.plot([t / 3600 for t in ts3], [cC3[-1][1 + j] for j in range(len(ts3))],
            "--", label="COMSOL（表面）")
    ax.set_xlabel("时间 t / h")
    ax.set_ylabel("水分浓度 C / (kg/kg)")
    ax.legend(fontsize=9)
    save(fig, "problem3", "fig4_validation.png", "问题3：两种数值方法对比")
    # ---------------- 问题 4 ----------------
    f4 = read_result(os.path.join(HERE, "result4.xlsx"))
    hdr4, ours4 = f4["Sheet1"]
    cur, cols = read_comsol_table(os.path.join(CSV, r"Problem4\comsol_p4_current.csv"))
    surf, _ = read_comsol_table(os.path.join(CSV, r"Problem4\comsol_p4_surface.csv"))
    tsN = cur[:, 0]
    common4 = [t for t in ours4[:, 0] if t in tsN]
    jO4 = [int(np.where(ours4[:, 0] == t)[0][0]) for t in common4]
    jC4 = [int(np.where(tsN == t)[0][0]) for t in common4]
    A4 = np.array([ours4[k, 1:-1] for k in jO4])
    B4 = np.array([cur[j, 1:] for j in jC4])
    A4s = np.array([ours4[k, -1] for k in jO4])
    B4s = np.array([surf[j, 1] for j in jC4])
    print("问题 4（到 %.2f h，12 个当前距离 + 表面，%d 个时刻）"
          % (ours4[-1, 0] / 3600, len(common4)))
    worst = max(worst, report("内部点", A4, B4, "kg/kg"), report("表面", A4s, B4s, "kg/kg"))
    m_late4 = np.array([t >= 21600.0 for t in common4])
    print("  排除最初 6 h 后：")
    report("内部点", A4[m_late4], B4[m_late4], "kg/kg")
    report("表面", A4s[m_late4], B4s[m_late4], "kg/kg")
    td4_ours = ours4[int(np.argmax(ours4[:, 1] < 0.15)), 0]
    td4_comsol = None
    for j, t in enumerate(tsN):
        if cur[j, 1] < 0.15:
            td4_comsol = t
            break
    print("  烘干时长：本文算法 %.2f h ； COMSOL %.2f h（相对偏差 %.2f%%）"
          % (td4_ours / 3600, td4_comsol / 3600,
             100 * abs(td4_ours - td4_comsol) / td4_comsol))
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.plot(ours4[:, 0] / 3600, ours4[:, 1], "-", label="本文算法（中心）")
    ax.plot(ours4[:, 0] / 3600, ours4[:, -1], "-", label="本文算法（表面）")
    ax.plot(tsN / 3600, cur[:, 1], "--", label="COMSOL（中心）")
    ax.plot(tsN / 3600, surf[:, 1], "--", label="COMSOL（表面）")
    ax.set_xlabel("时间 t / h")
    ax.set_ylabel("水分浓度 C / (kg/kg)")
    ax.legend(fontsize=9)
    save(fig, "problem4", "fig4_validation.png", "问题4：两种数值方法对比")
    print("")
    print("结论：两种独立数值方法的最大相对偏差为 %.2f%%。" % worst)


if __name__ == "__main__":
    main()
