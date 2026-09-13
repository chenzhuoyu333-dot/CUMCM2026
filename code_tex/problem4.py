# -*- coding: utf-8 -*-
# 问题 4：药材尺寸随水分流失收缩（附录 4 物性）
#
# 模型要点：半径按附件 2 的 R(t) 从 2.000 cm 收缩到 1.198 cm。
# 取材料坐标 X in [0, R0]（质点的初始径向位置），当前物理位置 r = X*R(t)/R0，
# 方程化为 dC/dt = (R0/R)^2 * (1/X) d/dX ( X D dC/dX )，
# 边界 -D dC/dX|_{X=R0} = (R/R0) hm (C - C_air)，
# 即扩散项放大 (R0/R)^2 倍、表面通量项系数变为 (R0^2/R)*hm。
#
# 输出：论文表 6（距离列按“当前离中心距离”取 0、0.5、1.0 cm，末列药材表面）、
#       烘干时长，以及 result4.xlsx（0~1.1 cm + 药材表面，每 60 s）。
# COMSOL 仿真仅用于独立验证，见 validate_comsol.py。
# 运行：python problem4.py
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (make_interp, make_radius_func, props_appendix4, read_air,
                    read_radius_cm, sample_shrink, solve, write_xlsx)

ATT = r"C:\Users\Wil\Desktop\CUMCM2026Problems\A题\附件"
T_END = 188760.0
DT = 10.0
N = 400
DISTS_FILE = [i / 10.0 for i in range(12)]
TAB_DISTS = [0.0, 0.5, 1.0]


def mark_dry_time(res):
    """最后一行取"四位小数显示上已明确达标"的最早时刻，避免 0.1500 的显示歧义"""
    for k, t in enumerate(res["t"]):
        if round(res["C"][k][0], 4) < 0.15:
            return t
    return res["t"][-1]


def main():
    Ta, Ca = read_air(os.path.join(ATT, "附件1.xlsx"))
    rcm = read_radius_cm(os.path.join(ATT, "附件2.xlsx"))
    env_T, env_C = make_interp(Ta), make_interp(Ca)
    R = make_radius_func(rcm)
    print("附件 2：半径 %.3f cm -> %.3f cm（t = 0 -> %d s）"
          % (rcm[min(rcm)], rcm[max(rcm)], max(rcm)))
    save = [float(t) for t in range(60, int(T_END) + 1, 60)]
    res = solve(dict(props=props_appendix4(), env_T=env_T, env_C=env_C,
                     R=R, t_end=T_END),
                N=N, dt=DT, save_times=save, shrink=True)
    print("自主求解完成：%d 个输出时刻，R(末) = %.5f m"
          % (len(res["t"]), res["R"][-1]))
    idx = {t: k for k, t in enumerate(res["t"])}
    t_dry = None
    for k, t in enumerate(res["t"]):
        if res["C"][k][0] < 0.15:
            t_dry = t
            break
    print("烘干时长（中心含水率首次 < 0.15）：%.0f s = %.2f h"
          % (t_dry, t_dry / 3600.0))
    print("")
    print("表 6  药材烘干过程的水分浓度（单位：kg/kg，距离为当前离中心距离）")
    print("  时间/h     0.0      0.5      1.0    药材表面")
    for h in [6 * i for i in range(1, 9)]:
        t = float(h * 3600)
        if t in idx:
            k = idx[t]
            C = sample_shrink(res["C"][k], TAB_DISTS, res["R"][k], N)
            print("  %-8d  %s  %8.4f"
                  % (h, "  ".join("%8.4f" % v for v in C), res["C"][k][-1]))
    MARK_DRY = mark_dry_time(res)
    k = idx[MARK_DRY]
    C = sample_shrink(res["C"][k], TAB_DISTS, res["R"][k], N)
    print("  %-8s  %s  %8.4f   <- 烘干结束时间 %.2f h"
          % ("结束", "  ".join("%8.4f" % v for v in C),
             res["C"][k][-1], MARK_DRY / 3600.0))
    hdr = ["时间\\到药材中心的距离"] + ["%g" % d for d in DISTS_FILE] + ["药材表面"]
    rows = []
    for k, t in enumerate(res["t"]):
        rows.append([t] + sample_shrink(res["C"][k], DISTS_FILE, res["R"][k], N)
                    + [res["C"][k][-1]])
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "result4.xlsx")
    write_xlsx(out, [("Sheet1", hdr, rows)])
    print("")
    print("已写出：%s（单工作表 Sheet1，%d 行 × %d 列）"
          % (out, len(rows), len(hdr)))


if __name__ == "__main__":
    main()
