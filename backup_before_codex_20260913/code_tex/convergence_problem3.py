# -*- coding: utf-8 -*-
"""问题 3 的网格/时间步收敛性实验（用 code/common.py 的同一套算法）"""
import os
import sys
import time

sys.path.insert(0, r"C:\Users\Wil\Desktop\cumcm_comsol\code")
from common import make_interp, props_appendix3, read_air, solve   # noqa: E402

ATT = r"C:\Users\Wil\Desktop\CUMCM2026Problems\A题\附件"
T_END = 206400.0
Ta, Ca = read_air(os.path.join(ATT, "附件1.xlsx"))
env_T, env_C = make_interp(Ta), make_interp(Ca)
save = [float(t) for t in range(60, int(T_END) + 1, 60)]


def run(N, dt):
    t0 = time.time()
    res = solve(dict(props=props_appendix3(), env_T=env_T, env_C=env_C, t_end=T_END),
                N=N, dt=dt, save_times=save)
    t_dry = None
    for k, t in enumerate(res["t"]):
        if res["C"][k][0] < 0.15:
            t_dry = t
            break
    surf = res["C"][res["t"].index(1800.0)][-1]
    print("N=%3d dt=%5.1f s -> t_dry = %8.0f s = %7.3f h ; 1800 s 表面 C = %.4f  (%.0f s)"
          % (N, dt, t_dry, t_dry / 3600.0, surf, time.time() - t0), flush=True)
    return t_dry / 3600.0


print("问题 3 收敛性（同一算法，只改网格与时间步）", flush=True)
for N, dt in ((400, 20.0), (400, 10.0), (400, 5.0), (200, 10.0), (100, 10.0)):
    run(N, dt)
