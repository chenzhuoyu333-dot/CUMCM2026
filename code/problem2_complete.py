# -*- coding: utf-8 -*-
"""问题二完整独立程序：变物性温度-水分双向耦合有限体积求解。"""
from pathlib import Path
import csv, math
import numpy as np
import matplotlib.pyplot as plt

R, T0, C0, H, HM = 0.02, 301.15, 2.55, 25.0, 8e-7
N, DT, T_END = 200, 2.0, 10800.0
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figure" / "Q2"; OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Microsoft YaHei","SimHei","DejaVu Sans"],
                     "axes.unicode_minus":False,"figure.dpi":150,"savefig.dpi":220})

AIR_T=np.array([0,300,600,900,1200,1500,1800,2400,3600,5400,7200,9000,10800,14400.])
AIR_TEMP=np.array([28,29.4,31.8,34.8,38.1,40.5,41.7,44.5,47,49.3,50,50.18,50.165,50.165])
AIR_C=np.full(AIR_T.size,0.04986)

def ambient(t):
    return float(np.interp(t,AIR_T,AIR_TEMP)),float(np.interp(t,AIR_T,AIR_C))

# 附录3状态相关物性
def k_fun(c): return 0.21+0.38*c/(c+1.0)
def rho_fun(c): return 650.0+128.0*c
def cp_fun(c): return 1450.0+2736.0*c/(c+1.0)
def d_fun(c,T): return 2.4e-3*math.exp(-0.45/max(float(c),1e-8))*math.exp(-3850.0/max(float(T),250.0))

def thomas(a,b,c,d):
    b,d=b.copy(),d.copy()
    for i in range(1,len(b)):
        m=a[i]/b[i-1]; b[i]-=m*c[i-1]; d[i]-=m*d[i-1]
    x=np.empty_like(d); x[-1]=d[-1]/b[-1]
    for i in range(len(b)-2,-1,-1): x[i]=(d[i]-c[i]*x[i+1])/b[i]
    return x

def assemble_and_solve(old,capacity,conductivity,h,value,dt,dr):
    """组装带柱面权重的全隐式三对角系统。"""
    n=len(old); faces=np.arange(1,n+1)*dr
    volume=math.pi*(faces**2-np.r_[0.0,faces[:-1]**2])
    cap=np.broadcast_to(np.asarray(capacity,float),old.shape)
    coef=np.broadcast_to(np.asarray(conductivity,float),old.shape)
    internal=2*math.pi*faces[:-1]*0.5*(coef[:-1]+coef[1:])/dr
    west=np.r_[0.0,internal]; east=np.r_[internal,2*math.pi*R*h]
    diagonal=volume*cap/dt+west+east
    rhs=volume*cap/dt*old; rhs[-1]+=east[-1]*value
    return thomas(-west,diagonal,-east,rhs)

def solve():
    dr=R/N; T=np.full(N,T0); C=np.full(N,C0)
    times=[0.0]; TH=[T.copy()]; CH=[C.copy()]
    for step in range(1,int(T_END/DT)+1):
        t=step*DT; Ta,Ca=ambient(t)
        # 顺序Picard：先用旧含水率更新热物性，再以新温度更新扩散系数。
        k=np.array([k_fun(x) for x in C]); capacity=np.array([rho_fun(x)*cp_fun(x) for x in C])
        T=assemble_and_solve(T,capacity,k,H,Ta+273.15,DT,dr)
        D=np.array([d_fun(c,temp) for c,temp in zip(C,T)])
        C=assemble_and_solve(C,1.0,D,HM,Ca,DT,dr)
        if step%30==0:
            times.append(t); TH.append(T.copy()); CH.append(C.copy())
    r=(np.arange(N)+0.5)*dr*100
    return np.array(times),r,np.array(TH)-273.15,np.array(CH)

def sample(v,r,x): return np.interp(x,r,v,left=v[0],right=v[-1])

def write_result(name,t,r,data):
    x=np.arange(0,2.01,.1)
    with (OUT/name).open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.writer(f);w.writerow(['time_s']+[f'r={p:.1f}cm' for p in x])
        for tt,row in zip(t,data): w.writerow([f'{tt:.0f}']+[f'{z:.6f}' for z in sample(row,r,x)])

def save(fig,name): fig.tight_layout();fig.savefig(OUT/name,bbox_inches='tight');plt.close(fig)

def plots(t,r,T,C):
    selected=[1800,5400,10800]; colors=['#1d6996','#38a6a5','#ed553b']
    for data,ylabel,name in [(T,'温度 (°C)','q2_T_profile.png'),(C,'含水率 (kg/kg)','q2_C_profile.png')]:
        fig,ax=plt.subplots(figsize=(5.7,3.6))
        for tt,color in zip(selected,colors):
            j=np.argmin(abs(t-tt));ax.plot(r,data[j],lw=2.2,color=color,label=f'{tt/3600:g} h')
        ax.set(xlabel='到中心距离 r (cm)',ylabel=ylabel);ax.grid(alpha=.22);ax.legend(frameon=False);save(fig,name)
    fig,ax=plt.subplots(figsize=(6.2,3.6));im=ax.imshow(T,origin='lower',aspect='auto',extent=[0,2,0,3],cmap='turbo');ax.set(xlabel='到中心距离 r (cm)',ylabel='时间 (h)');fig.colorbar(im,ax=ax,label='温度 (°C)');save(fig,'q2_field.png')
    fig,ax=plt.subplots(figsize=(5.7,3.6));ax.plot(t/3600,C[:,0],lw=2.2,label='中心');ax.plot(t/3600,C[:,-1],lw=2.2,label='表面');ax.set(xlabel='时间 (h)',ylabel='含水率 (kg/kg)');ax.grid(alpha=.22);ax.legend(frameon=False);save(fig,'q2_center_surface.png')
    c=np.linspace(.05,2.55,250); kval=k_fun(c);rho=rho_fun(c);cp=cp_fun(c)
    fig,ax=plt.subplots(figsize=(5.8,3.6));ax.plot(c,kval,lw=2.2,color='#38a6a5',label='导热系数 k');ax.set(xlabel='含水率 C (kg/kg)',ylabel='k (W/(m·K))');ax.grid(alpha=.22)
    ax2=ax.twinx();ax2.plot(c,rho,lw=2,color='#ed553b',label='密度 ρ');ax2.plot(c,cp,lw=2,ls='--',color='#1d6996',label='比热容 cp');ax2.set_ylabel('ρ (kg/m³)，cp (J/(kg·K))');
    lines=ax.lines+ax2.lines;ax.legend(lines,[x.get_label() for x in lines],frameon=False,loc='center right');save(fig,'q2_properties.png')

def main():
    t,r,T,C=solve();write_result('result2_temperature.csv',t,r,T);write_result('result2_moisture.csv',t,r,C);plots(t,r,T,C)
    x=[0,.5,1,1.5,2]
    print('问题二计算完成；3 h抽样结果');print('r/cm:',x);print('T/°C:',np.round(sample(T[-1],r,x),4));print('C:',np.round(sample(C[-1],r,x),4));print('输出目录：',OUT)

if __name__=='__main__': main()
