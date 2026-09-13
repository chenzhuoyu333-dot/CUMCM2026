# -*- coding: utf-8 -*-
"""问题四稳定运行入口，同时输出 result4.xlsx。"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
import problem4_complete as core
from xlsx_output import write_table

_PLOTS = core.plots

def safe_plots(t, X, C, RH, dry, fixed, coarse):
    tf, Xf, Tf, Cf, Rf, dryf = fixed
    if dryf is None:
        fixed = (tf, Xf, Tf, Cf, Rf, tf[-1])
    _PLOTS(t, X, C, RH, dry, fixed, coarse)

def main():
    result = core.solve(True, 300)
    fixed = core.solve(False, 220)
    coarse = core.solve(True, 150)
    t, X, T, C, RH, dry = result
    core.write_result(t, X, C, RH)
    pos = __import__('numpy').arange(0.0, 1.11, 0.1)
    header = ['time_s', 'radius_cm'] + [f'r={x:.1f}cm' for x in pos] + ['surface']
    rows = [[float(tt), float(Rt * 100)] + list(core.physical_sample(row, X, pos, Rt)) + [float(row[-1])]
            for tt, row, Rt in zip(t, C, RH)]
    write_table(core.ROOT / 'result4.xlsx', [('水分浓度', header, rows)])
    safe_plots(t, X, C, RH, dry, fixed, coarse)
    print(f'问题四计算完成；考虑收缩={dry/3600:.4f} h，已写出 {core.ROOT / "result4.xlsx"}')

if __name__ == '__main__':
    main()
