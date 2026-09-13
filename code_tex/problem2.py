# -*- coding: utf-8 -*-
"""
问题 2：整个烘干过程的前 3 小时（附录 3 变物性，温度—水分双向耦合）

自主用有限体积法求解 0~10800 s，输出论文表 3、表 4 与 result2.xlsx。
（COMSOL 仿真仅作为独立验证手段，见 validate_comsol.py。）

运行：python problem2.py
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import make_interp, props_appendix3, read_air, sample_fixed, solve, write_xlsx

ATT = r"C:\Users\Wil\Desktop\CUMCM2026Problems\A题\附件"
T_END = 10800.0
DT = 1.0
N = 400
DISTS = [i / 10.0 for i in range(21)]
TAB_TIMES = [1800, 3600, 5400, 7200, 9000, 10800]
TAB_DISTS = [0.0, 0.5, 1.0, 1.5, 2.0]


def main():
    Ta, Ca = read_air(os.path.join(ATT, "附件1.xlsx"))
    env_T, env_C = make_interp(Ta), make_interp(Ca)
    res = solve(dict(props=props_appendix3(), env_T=env_T, env_C=env_C, t_end=T_END),
                N=N, dt=DT, save_times=[float(t) for t in range(1, int(T_END) + 1)])
    idx = {t: k for k, t in enumerate(res["t"])}
    print("自主求解完成：%d 个时刻（1 s 间隔）" % len(res["t"]))

    print("")
    print("表 3  3 小时内药材的温度（单位：℃）")
    print("  时间/h     0.0      0.5      1.0      1.5      2.0")
    for t in TAB_TIMES:
        T = [v - 273.15 for v in sample_fixed(res["T"][idx[float(t)]], TAB_DISTS, N)]
        print("  %-8.1f  %s" % (t / 3600.0, "  ".join("%8.4f" % v for v in T)))
    print("")
    print("表 4  3 小时内药材的水分浓度（单位：kg/kg）")
    print("  时间/h     0.0      0.5      1.0      1.5      2.0")
    for t in TAB_TIMES:
        C = sample_fixed(res["C"][idx[float(t)]], TAB_DISTS, N)
        print("  %-8.1f  %s" % (t / 3600.0, "  ".join("%8.4f" % v for v in C)))

    hdr = ["时间\\到药材中心的距离"] + ["%.1f" % d for d in DISTS]
    rows_T, rows_C = [], []
    for k, t in enumerate(res["t"]):
        rows_T.append([t] + [v - 273.15 for v in sample_fixed(res["T"][k], DISTS, N)])
        rows_C.append([t] + sample_fixed(res["C"][k], DISTS, N))
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "result2.xlsx")
    write_xlsx(out, [("温度", hdr, rows_T), ("水分浓度", hdr, rows_C)])
    print("")
    print("已写出：%s（%d 行 × %d 列，两个工作表）" % (out, len(rows_T), len(hdr)))


if __name__ == "__main__":
    main()
