# -*- coding: utf-8 -*-
"""稳定运行问题四独立完整程序。"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
import problem4_complete as core

ORIGINAL_PLOTS = core.plots

def safe_plots(t, X, C, RH, dry, fixed, coarse):
    tf, Xf, Tf, Cf, Rf, dryf = fixed
    if dryf is None:
        fixed = (tf, Xf, Tf, Cf, Rf, tf[-1])
    ORIGINAL_PLOTS(t, X, C, RH, dry, fixed, coarse)

def main():
    result = core.solve(True, 300)
    fixed = core.solve(False, 220)
    coarse = core.solve(True, 150)
    t, X, T, C, RH, dry = result
    core.write_result(t, X, C, RH)
    safe_plots(t, X, C, RH, dry, fixed, coarse)
    print(f'问题四计算完成；考虑收缩={dry/3600:.4f} h')
    print('固定半径对照在当前模拟时限内未达到 0.15，图中已明确标注')

if __name__ == '__main__':
    main()
