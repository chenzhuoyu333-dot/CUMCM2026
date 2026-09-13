# -*- coding: utf-8 -*-
"""
药材烘干问题 —— 自主求解算法核心模块（问题 1~4 共用）

物理模型（一维轴对称，长 25 cm、半径 R=2 cm 的圆柱）：
    水分： dC/dt = (1/r) * d/dr ( r * D(C,T) * dC/dr )
    温度： rho(C)*cp(C)*dT/dt = (1/r) * d/dr ( r * k(C) * dT/dr )
    边界： r=R:  -k dT/dr = h*(T - T_air(t)),  -D dC/dr = h_m*(C - C_air(t))
          r=0:  对称（零通量）
    初始： T=28 ℃ = 301.15 K, C=2.55 kg/kg

数值方法：
    * 空间：守恒型有限体积法（控制体柱面权重，轴心 r=0 处面积为零，奇点自然消解）
    * 时间：全隐式（无条件稳定）
    * 求解：三对角方程组的 Thomas（追赶）法
    * 非线性：物性用上一时间步的值滞后一步（误差 O(dt)）

问题 4 的尺寸收缩：
    半径按附件 2 的 R(t) 收缩。取材料坐标 X∈[0,R0]（质点的初始径向位置），
    当前物理位置 r = X*R(t)/R0。方程变为
        dC/dt = (R0/R)^2 * (1/X) * d/dX ( X * D * dC/dX )
        -D dC/dX |_{X=R0} = (R/R0) * h_m * (C - C_air)
    即：扩散项被放大 (R0/R)^2 倍、表面通量项系数变为 (R0^2/R)*h_m。
"""

import math
import re
import zipfile
import xml.etree.ElementTree as ET

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}

R0 = 0.02          # 药材初始半径 m
L = 0.25           # 药材长度 m
H_CONV = 25.0      # 对流换热系数 W/(m^2*K)
H_MASS = 8e-7      # 对流传质系数 m/s
T0 = 301.15        # 初始温度 K（28 ℃）
C0 = 2.55          # 初始干基含水率 kg/kg


# ----------------------------------------------------------------- 读附件
def read_xlsx_columns(path):
    """读取 xlsx 第一张工作表，返回二维列表（字符串/浮点混合）"""
    z = zipfile.ZipFile(path)
    strs = []
    if "xl/sharedStrings.xml" in z.namelist():
        root = ET.fromstring(z.read("xl/sharedStrings.xml"))
        for si in root.findall("m:si", NS):
            strs.append("".join(t.text or "" for t in si.iter(
                "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}t")))
    sh = ET.fromstring(z.read("xl/worksheets/sheet1.xml"))
    rows = []
    for row in sh.findall(".//m:sheetData/m:row", NS):
        cells = {}
        for c in row.findall("m:c", NS):
            v = c.find("m:v", NS)
            if v is None:
                continue
            col = c.get("r")[0]
            cells[col] = strs[int(v.text)] if c.get("t") == "s" else float(v.text)
        rows.append(cells)
    return rows


def read_air(path):
    """附件 1：时间 s、烘房温度 ℃、烘房水分浓度 kg/kg"""
    rows = read_xlsx_columns(path)
    Ta, Ca = {}, {}
    for cells in rows:
        try:
            t = int(float(cells["A"]))
            Ta[t] = float(cells["B"])
            Ca[t] = float(cells["C"])
        except (KeyError, ValueError, TypeError):
            pass
    return Ta, Ca


def read_radius_cm(path):
    """附件 2：时间 s、半径 cm"""
    rows = read_xlsx_columns(path)
    out = {}
    for cells in rows:
        try:
            out[int(float(cells["A"]))] = float(cells["B"])
        except (KeyError, ValueError, TypeError):
            pass
    return out


def make_interp(dic, tmax=None):
    """线性插值 + 常数外推（烘房温度/水分浓度、半径都用它）"""
    ts = sorted(dic)
    hi = tmax if tmax is not None else ts[-1]

    def f(t):
        if t <= ts[0]:
            return dic[ts[0]]
        if t >= hi:
            return dic[hi] if hi in dic else dic[ts[-1]]
        for k in range(len(ts) - 1):
            if ts[k] <= t <= ts[k + 1]:
                if ts[k + 1] == ts[k]:
                    return dic[ts[k]]
                w = (t - ts[k]) / (ts[k + 1] - ts[k])
                return dic[ts[k]] + (dic[ts[k + 1]] - dic[ts[k]]) * w
        return dic[ts[-1]]

    return f


def make_radius_func(rcm):
    """R(t)（单位 m）；t 超过附件范围后取末值（常数外推）"""
    ts = sorted(rcm)
    tmax = ts[-1]

    def R(t):
        if t >= tmax:
            return rcm[tmax] / 100.0
        return make_interp(rcm)(t) / 100.0

    return R


# ------------------------------------------------------------- 物性（附录）
def props_appendix2():
    """问题 1：附录 2"""
    return dict(
        name="附录 2",
        k=lambda c: 0.36,
        rho=lambda c: 820.0,
        cp=lambda c: 2600.0,
        D=lambda c, T: 7e-9 * math.exp(-0.89 / max(c, 1e-6)),
    )


def props_appendix3():
    """问题 2、3：附录 3"""
    return dict(
        name="附录 3",
        k=lambda c: 0.21 + 0.38 * c / (c + 1),
        rho=lambda c: 650.0 + 128.0 * c,
        cp=lambda c: 1450.0 + 2736.0 * c / (c + 1),
        D=lambda c, T: 2.4e-3 * math.exp(-0.45 / max(c, 1e-6)) * math.exp(-3850.0 / T),
    )


def props_appendix4():
    """问题 4：附录 4"""
    return dict(
        name="附录 4",
        k=lambda c: 0.12 + 0.20 * c / (c + 1),
        rho=lambda c: 760.0 + 90.0 * c,
        cp=lambda c: 1850.0 + 2150.0 * c / (c + 1),
        D=lambda c, T: 4.2e-4 * math.exp(-0.30 / max(c, 1e-6)) * math.exp(-3850.0 / T),
    )


# ------------------------------------------------------------------ 追赶法
def thomas(a, b, c, d, N):
    cp = [0.0] * N
    dp = [0.0] * N
    cp[0] = c[0] / b[0]
    dp[0] = d[0] / b[0]
    for i in range(1, N):
        m = b[i] - a[i] * cp[i - 1]
        cp[i] = c[i] / m if i < N - 1 else 0.0
        dp[i] = (d[i] - a[i] * dp[i - 1]) / m
    x = [0.0] * N
    x[N - 1] = dp[N - 1]
    for i in range(N - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x


def _step_fixed(N, dr, T, C, dt, props, Ta, Ca, h=H_CONV, hm=H_MASS, R=None):
    """固定半径的一个时间步（隐式），返回新的 T、C"""
    R = R0 if R is None else R
    rf = [(i + 1) * dr for i in range(N)]
    V = [math.pi * (rf[i] ** 2 - (rf[i - 1] ** 2 if i > 0 else 0.0)) for i in range(N)]
    # 温度
    kk = [props["k"](c) for c in C]
    rc = [props["rho"](c) * props["cp"](c) for c in C]
    kf = [0.5 * (kk[i] + kk[i + 1]) for i in range(N - 1)] + [kk[-1]]
    kb = [0.0] + [0.5 * (kk[i - 1] + kk[i]) for i in range(1, N)]
    A = [2 * math.pi * rf[i - 1] * kb[i] / dr if i > 0 else 0.0 for i in range(N)]
    B = [2 * math.pi * rf[i] * kf[i] / dr for i in range(N - 1)] + [2 * math.pi * R * h]
    a = [0.0] * N; b = [0.0] * N; c = [0.0] * N; e = [0.0] * N
    for i in range(N):
        b[i] = V[i] * rc[i] / dt + A[i] + B[i]
        e[i] = V[i] * rc[i] / dt * T[i]
        if i > 0:
            a[i] = -A[i]
        if i < N - 1:
            c[i] = -B[i]
        else:
            e[i] += B[i] * Ta
    T = thomas(a, b, c, e, N)
    # 水分
    DD = [props["D"](C[i], T[i]) for i in range(N)]
    Df = [0.5 * (DD[i] + DD[i + 1]) for i in range(N - 1)] + [DD[-1]]
    Db = [0.0] + [0.5 * (DD[i - 1] + DD[i]) for i in range(1, N)]
    A = [2 * math.pi * rf[i - 1] * Db[i] / dr if i > 0 else 0.0 for i in range(N)]
    B = [2 * math.pi * rf[i] * Df[i] / dr for i in range(N - 1)] + [2 * math.pi * R * hm]
    a = [0.0] * N; b = [0.0] * N; c = [0.0] * N; e = [0.0] * N
    for i in range(N):
        b[i] = V[i] / dt + A[i] + B[i]
        e[i] = V[i] / dt * C[i]
        if i > 0:
            a[i] = -A[i]
        if i < N - 1:
            c[i] = -B[i]
        else:
            e[i] += B[i] * Ca
    C = thomas(a, b, c, e, N)
    return T, C


def solve(case, N=400, dt=10.0, t_end=None, save_every=None, save_times=None,
          shrink=False, verbose=True):
    """
    求解主程序。
    case: dict(env_T=..., env_C=..., props=..., R=None)
    shrink=True 时按材料坐标做收缩模型（问题 4）。
    返回 dict(t=[...], T=[[...]], C=[[...]], R=[...])，只保存需要输出的时刻。
    """
    props = case["props"]
    env_T = case["env_T"]
    env_C = case["env_C"]
    R = case.get("R")
    t_end = case.get("t_end", t_end)
    out = {"t": [], "T": [], "C": [], "R": []}
    want = sorted(save_times) if save_times else None

    if not shrink:
        dr = R0 / N
        T = [T0] * N
        C = [C0] * N
        n = 0
        if want and want[0] == 0:
            out["t"].append(0.0); out["T"].append(list(T)); out["C"].append(list(C))
            out["R"].append(R0)
        while n * dt < t_end - 1e-9:
            n += 1
            tt = n * dt
            if tt > t_end:
                tt = t_end
            T, C = _step_fixed(N, dr, T, C, dt, props,
                               env_T(tt) + 273.15, env_C(tt))
            if want is not None:
                if tt in want:
                    out["t"].append(tt); out["T"].append(list(T)); out["C"].append(list(C))
                    out["R"].append(R0)
            elif save_every and abs(tt % save_every) < 1e-9:
                out["t"].append(tt); out["T"].append(list(T)); out["C"].append(list(C))
                out["R"].append(R0)
        return out

    # ---------------- 收缩模型（材料坐标） ----------------
    ds = R0 / N
    sf = [(i + 1) * ds for i in range(N)]
    V = [math.pi * (sf[i] ** 2 - (sf[i - 1] ** 2 if i > 0 else 0.0)) for i in range(N)]
    T = [T0] * N
    C = [C0] * N
    n = 0
    if want and want[0] == 0:
        out["t"].append(0.0); out["T"].append(list(T)); out["C"].append(list(C))
        out["R"].append(R(0.0))
    while n * dt < t_end - 1e-9:
        n += 1
        tt = min(n * dt, t_end)
        Rt = R(tt)
        sc = (R0 / Rt) ** 2                     # 扩散项放大系数
        bc = 2 * math.pi * R0 * R0 / Rt         # 表面通量项系数
        Ta, Ca = env_T(tt) + 273.15, env_C(tt)
        # 温度
        kk = [props["k"](c) * sc for c in C]
        rc = [props["rho"](c) * props["cp"](c) for c in C]
        kf = [0.5 * (kk[i] + kk[i + 1]) for i in range(N - 1)] + [kk[-1]]
        kb = [0.0] + [0.5 * (kk[i - 1] + kk[i]) for i in range(1, N)]
        A = [2 * math.pi * sf[i - 1] * kb[i] / ds if i > 0 else 0.0 for i in range(N)]
        B = [2 * math.pi * sf[i] * kf[i] / ds for i in range(N - 1)] + [bc * H_CONV]
        a = [0.0] * N; b = [0.0] * N; c = [0.0] * N; e = [0.0] * N
        for i in range(N):
            b[i] = V[i] * rc[i] / dt + A[i] + B[i]
            e[i] = V[i] * rc[i] / dt * T[i]
            if i > 0:
                a[i] = -A[i]
            if i < N - 1:
                c[i] = -B[i]
            else:
                e[i] += B[i] * Ta
        T = thomas(a, b, c, e, N)
        # 水分
        DD = [props["D"](C[i], T[i]) * sc for i in range(N)]
        Df = [0.5 * (DD[i] + DD[i + 1]) for i in range(N - 1)] + [DD[-1]]
        Db = [0.0] + [0.5 * (DD[i - 1] + DD[i]) for i in range(1, N)]
        A = [2 * math.pi * sf[i - 1] * Db[i] / ds if i > 0 else 0.0 for i in range(N)]
        B = [2 * math.pi * sf[i] * Df[i] / ds for i in range(N - 1)] + [bc * H_MASS]
        a = [0.0] * N; b = [0.0] * N; c = [0.0] * N; e = [0.0] * N
        for i in range(N):
            b[i] = V[i] / dt + A[i] + B[i]
            e[i] = V[i] / dt * C[i]
            if i > 0:
                a[i] = -A[i]
            if i < N - 1:
                c[i] = -B[i]
            else:
                e[i] += B[i] * Ca
        C = thomas(a, b, c, e, N)
        if want is not None:
            if tt in want:
                out["t"].append(tt); out["T"].append(list(T)); out["C"].append(list(C))
                out["R"].append(Rt)
        elif save_every and abs(tt % save_every) < 1e-9:
            out["t"].append(tt); out["T"].append(list(T)); out["C"].append(list(C))
            out["R"].append(Rt)
    return out


# ------------------------------------------------------- 取样与输出工具
def sample_fixed(C, dists_cm, N=400):
    """在"到中心距离"的固定网格上取样（问题 1~3：坐标即为材料坐标）"""
    dr = R0 / N
    xs = [(i + 0.5) * dr for i in range(N)]        # 控制体中心
    out = []
    for d in dists_cm:
        x = d / 100.0
        if x <= 0:
            out.append(C[0])
        elif x >= R0:
            out.append(C[-1])
        else:
            k = min(range(N - 1), key=lambda i: 1 if not (xs[i] <= x <= xs[i + 1]) else 0)
            for i in range(N - 1):
                if xs[i] <= x <= xs[i + 1]:
                    w = (x - xs[i]) / (xs[i + 1] - xs[i])
                    out.append(C[i] + (C[i + 1] - C[i]) * w)
                    break
            else:
                out.append(C[-1])
    return out


def sample_shrink(C, dists_cm, R_t, N=400):
    """问题 4：按"当前离中心距离"取样（材料坐标 X = d*R0/R(t)）"""
    ds = R0 / N
    xs = [(i + 0.5) * ds for i in range(N)]
    out = []
    for d in dists_cm:
        if d >= R_t * 100 - 1e-9:
            out.append(C[-1])
            continue
        X = (d / 100.0) * R0 / R_t
        if X <= 0:
            out.append(C[0])
            continue
        for i in range(N - 1):
            if xs[i] <= X <= xs[i + 1]:
                w = (X - xs[i]) / (xs[i + 1] - xs[i])
                out.append(C[i] + (C[i + 1] - C[i]) * w)
                break
        else:
            out.append(C[-1])
    return out


def col_letter(i):
    s = ""
    i += 1
    while i:
        i, r = divmod(i - 1, 26)
        s = chr(65 + r) + s
    return s


def write_xlsx(path, sheets):
    """把 [(表名, 表头, 数据行), ...] 写成 xlsx（不依赖第三方库）"""
    ct = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
          '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">',
          '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>',
          '<Default Extension="xml" ContentType="application/xml"/>',
          '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>',
          '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>']
    for i in range(len(sheets)):
        ct.append('<Override PartName="/xl/worksheets/sheet%d.xml" '
                  'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>' % (i + 1))
    ct.append("</Types>")
    rels = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
            "</Relationships>")
    wb = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
          '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
          'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>']
    for i, (name, _, _) in enumerate(sheets):
        wb.append('<sheet name="%s" sheetId="%d" r:id="rId%d"/>' % (name, i + 1, i + 1))
    wb.append("</sheets></workbook>")
    wbrels = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
              '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">']
    for i in range(len(sheets)):
        wbrels.append('<Relationship Id="rId%d" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet%d.xml"/>' % (i + 1, i + 1))
    wbrels.append('<Relationship Id="rId%d" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>' % (len(sheets) + 1))
    wbrels.append("</Relationships>")
    styles = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
              '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
              '<fonts count="1"><font><sz val="11"/><name val="Calibri"/></font></fonts>'
              '<fills count="1"><fill><patternFill patternType="none"/></fill></fills>'
              '<borders count="1"><border/></borders><cellStyleXfs count="1"><xf/></cellStyleXfs>'
              '<cellXfs count="1"><xf/></cellXfs></styleSheet>')
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", "".join(ct))
        z.writestr("_rels/.rels", rels)
        z.writestr("xl/workbook.xml", "".join(wb))
        z.writestr("xl/_rels/workbook.xml.rels", "".join(wbrels))
        z.writestr("xl/styles.xml", styles)
        for i, (_, header, rows) in enumerate(sheets):
            out = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
                   '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData>']
            out.append('<row r="1">')
            for j, h in enumerate(header):
                out.append('<c r="%s1" t="inlineStr"><is><t>%s</t></is></c>' % (col_letter(j), str(h)))
            out.append("</row>")
            for k, row in enumerate(rows):
                rr = k + 2
                out.append('<row r="%d">' % rr)
                for j, v in enumerate(row):
                    if j == 0:
                        out.append('<c r="%s%d"><v>%d</v></c>' % (col_letter(j), rr, int(v)))
                    else:
                        out.append('<c r="%s%d"><v>%.4f</v></c>' % (col_letter(j), rr, v))
                out.append("</row>")
            out.append("</sheetData></worksheet>")
            z.writestr("xl/worksheets/sheet%d.xml" % (i + 1), "".join(out))
