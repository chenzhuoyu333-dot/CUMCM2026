# -*- coding: utf-8 -*-
"""灵敏度柱状图（使用论文统一样式）"""
import os
import sys

sys.path.insert(0, r"C:\Users\Wil\Documents\ChatGPT\CUMCM\code2")
import paper_style as ps
import matplotlib.pyplot as plt

names = ["传热系数 h\n(±20%)", "传质系数 hm\n(±20%)", "初始含水率 C0\n(±5%)",
         "扩散系数前置系数\n(±20%)", "环境温度\n(+1 K)", "收缩项口径\n(问题四)"]
vals = [0.0015, 0.1216, 0.0850, 0.9188, 0.0300 * 0 + 0.032, 0.088]
cols = [ps.C_GRAY, ps.C_GOLD, ps.C_PURPLE, ps.C_RED, ps.C_GREEN, ps.C_BLUE]

fig, ax = ps.new_fig(6.4, 3.6)
x = range(len(names))
ax.bar(x, vals, width=0.62, color=cols, alpha=0.9, edgecolor="#333333", lw=0.7)
for i, v in enumerate(vals):
    ax.text(i, v + 0.025, "%.3f" % v, ha="center", va="bottom", fontsize=8.6)
ax.set_xticks(list(x))
ax.set_xticklabels(names, fontsize=8.4)
ax.set_ylabel("无量纲灵敏度（绝对值）")
ax.set_ylim(0, 1.06)
ax.axhline(0.1, color=ps.C_GRAY, ls=":", lw=0.9)
ax.text(5.45, 0.115, "0.1（分界线）", fontsize=8, color=ps.C_GRAY, ha="right")
ps.style(ax)
ps.save(fig, "paper", "sens_bar")
