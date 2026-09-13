# -*- coding: utf-8 -*-
"""
问题 3：确定烘干所需时间（附录 3 物性，几何不变）

判据：题目要求"药材各处的水分浓度均低于 0.15 kg/kg"。
      表面先干、中心最后干（水分由中心向表面单向扩散，含水率沿半径单调递减），
      故只需判断中心含水率；烘干时长 t_dry = 中心含水率首次低于 0.15 的时刻。

输出：论文表 5、烘干时长，以及 result3.xlsx（每 60 s × 每 0.1 cm）。
（COMSOL 仿真仅作为独立验证手段，见 validate_comsol.py。）

运行：python problem3.py
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import C0, make_interp, props_appendix3, read_air, sample_fixed, solve, write_xlsx

ATT = r"C:\Users\Wil\Desktop\CUMCM2026Problems\A题\附件"
T_END = 206400.0          # 覆盖到烘干结束（57 h 左右）
DT = 10.0
N = 400
DISTS = [i / 10.0 for i in range(21)]
TAB_DISTS = [0.0, 0.5, 1.0, 1.5, 2.0]


def main():
    Ta, Ca = read_air(os.path.join(ATT, "附件1.xlsx"))
    env_T, env_C = make_interp(Ta), make_interp(Ca)
    save = [float(t) for t in range(60, int(T_END) + 1, 60)]
    res = solve(dict(props=props_appendix3(), env_T=env_T, env_C=env_C, t_end=T_END),
                N=N, dt=DT, save_times=save)
    print("自主求解完成：%d 个输出时刻（60 s 间隔，到 %.2f h）"
          % (len(res["t"]), res["t"][-1] / 3600.0))

    # 烘干时长：中心含水率首次低于 0.15
    t_dry = None
    for k, t in enumerate(res["t"]):
        if res["C"][k][0] < 0.15:
            t_dry = t
            break
    print("")
    print("烘干时长（中心含水率首次 < 0.15）：%.0f s = %.2f h" % (t_dry, t_dry / 3600.0))
    print("  参考：0.15 前一个输出时刻中心 = %.5f"
          % res["C"][res["t"].index(t_dry - 60)][0])

    # 论文表 5
    idx = {t: k for k, t in enumerate(res["t"])}
    print("")
    print("表 5  药材烘干过程的水分浓度（单位：kg/kg）")
    print("  时间/h     0.0      0.5      1.0      1.5      2.0")
    hours = [6 * i for i in range(1, 10)]
    # 最后一行取"四位小数显示上已明确达标"的最早时刻，避免 0.1500 的显示歧义
    t_last = None
    for k, t in enumerate(res["t"]):
        if round(res["C"][k][0], 4) < 0.15:
            t_last = t
            break
    for h in hours:
        t = h * 3600.0
        if t in idx:
            C = sample_fixed(res["C"][idx[t]], TAB_DISTS, N)
            print("  %-8d  %s" % (h, "  ".join("%8.4f" % v for v in C)))
    C = sample_fixed(res["C"][idx[t_last]], TAB_DISTS, N)
    print("  %-8s  %s   ← 烘干结束时间 %.2f h"
          % ("结束", "  ".join("%8.4f" % v for v in C), t_last / 3600.0))

    hdr = ["时间\\到药材中心的距离"] + ["%.1f" % d for d in DISTS]
    rows = [[t] + sample_fixed(res["C"][k], DISTS, N) for k, t in enumerate(res["t"])]
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "result3.xlsx")
    write_xlsx(out, [("Sheet1", hdr, rows)])
    print("")
    print("已写出：%s（单工作表 Sheet1，%d 行 × %d 列）" % (out, len(rows), len(hdr)))


if __name__ == "__main__":
    main()
