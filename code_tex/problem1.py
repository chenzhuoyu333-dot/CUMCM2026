# -*- coding: utf-8 -*-
"""
问题 1：预热平衡阶段（附录 2 物性，常数物性）

本程序由模型方程出发，自主用有限体积法求解 0~1800 s 的温度场与水分浓度场，
输出论文表 1、表 2 与结果文件 result1.xlsx。
（COMSOL 仿真仅作为后续独立验证手段，见 validate_comsol.py。）

运行：python problem1.py
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (R0, read_air, make_interp, props_appendix2, solve,
                    sample_fixed, write_xlsx)

ATT = r"C:\Users\Wil\Desktop\CUMCM2026Problems\A题\附件"
T_END = 1800.0
DT = 1.0
N = 400
DISTS = [i / 10.0 for i in range(21)]          # 0, 0.1, ..., 2.0 cm
TAB_TIMES = [100, 300, 600, 900, 1200, 1500, 1800]
TAB_DISTS = [0.0, 0.5, 1.0, 1.5, 2.0]


def main():
    Ta, Ca = read_air(os.path.join(ATT, "附件1.xlsx"))
    env_T = make_interp(Ta)
    env_C = make_interp(Ca)
    print("附件 1：%d 个时间点，t = %d ~ %d s" % (len(Ta), min(Ta), max(Ta)))

    res = solve(dict(props=props_appendix2(), env_T=env_T, env_C=env_C, t_end=T_END),
                N=N, dt=DT, save_times=[float(t) for t in range(1, int(T_END) + 1)])
    print("自主求解完成：%d 个时刻" % len(res["t"]))

    # ---------------- 论文表 1、表 2 ----------------
    idx = {t: k for k, t in enumerate(res["t"])}
    print("")
    print("表 1  30 分钟内药材的温度（单位：℃）")
    print("  时间/s     0.0      0.5      1.0      1.5      2.0")
    for t in TAB_TIMES:
        T = [v - 273.15 for v in sample_fixed(res["T"][idx[float(t)]], TAB_DISTS, N)]
        print("  %-6d  %s" % (t, "  ".join("%8.4f" % v for v in T)))
    print("")
    print("表 2  30 分钟内药材的水分浓度（单位：kg/kg）")
    print("  时间/s     0.0      0.5      1.0      1.5      2.0")
    for t in TAB_TIMES:
        C = sample_fixed(res["C"][idx[float(t)]], TAB_DISTS, N)
        print("  %-6d  %s" % (t, "  ".join("%8.4f" % v for v in C)))

    # ---------------- result1.xlsx ----------------
    hdr = ["时间\\到药材中心的距离"] + ["%.1f" % d for d in DISTS]
    rows_T, rows_C = [], []
    for k, t in enumerate(res["t"]):
        rows_T.append([t] + [v - 273.15 for v in sample_fixed(res["T"][k], DISTS, N)])
        rows_C.append([t] + sample_fixed(res["C"][k], DISTS, N))
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "result1.xlsx")
    write_xlsx(out, [("温度", hdr, rows_T), ("水分浓度", hdr, rows_C)])
    print("")
    print("已写出：%s（温度 / 水分浓度 两个工作表，%d 行 × %d 列）"
          % (out, len(rows_T), len(hdr)))


if __name__ == "__main__":
    main()
