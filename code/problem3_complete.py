# -*- coding: utf-8 -*-
"""问题三完整独立程序：固定半径下以全域最大含水率确定烘干终点。"""
from pathlib import Path
import csv, math
import numpy as np
import matplotlib.pyplot as plt
from xlsx_output import write_table

R,T0,C0,H,HM=0.02,301.15,2.55,25.0,8e-7
N,DT,T_LIMIT,THRESHOLD=300,30.0,230000.0,0.15
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'figure'/'Q3';OUT.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],
                     'axes.unicode_minus':False,'figure.dpi':150,'savefig.dpi':220})
AT=np.array([0,300,600,900,1200,1500,1800,2400,3600,5400,7200,9000,10800,14400.])
TEMP=np.array([28,29.4,31.8,34.8,38.1,40.5,41.7,44.5,47,49.3,50,50.18,50.165,50.165]);AC=np.full(AT.size,.04986)

def ambient(t): return float(np.interp(t,AT,TEMP)),float(np.interp(t,AT,AC))
def kfun(c): return .21+.38*c/(c+1)
def rhofun(c): return 650+128*c
def cpfun(c): return 1450+2736*c/(c+1)
def dfun(c,T): return 2.4e-3*math.exp(-.45/max(float(c),1e-8))*math.exp(-3850/max(float(T),250))

def thomas(a,b,c,d):
    b,d=b.copy(),d.copy()
    for i in range(1,len(b)):
        f=a[i]/b[i-1];b[i]-=f*c[i-1];d[i]-=f*d[i-1]
    x=np.empty_like(d);x[-1]=d[-1]/b[-1]
    for i in range(len(b)-2,-1,-1):x[i]=(d[i]-c[i]*x[i+1])/b[i]
    return x

def fvm(old,cap,coef,h,external,dt,dr):
    n=len(old);face=np.arange(1,n+1)*dr;vol=math.pi*(face**2-np.r_[0.,face[:-1]**2])
    cap=np.broadcast_to(np.asarray(cap,float),old.shape);coef=np.broadcast_to(np.asarray(coef,float),old.shape)
    internal=2*math.pi*face[:-1]*.5*(coef[:-1]+coef[1:])/dr
    west=np.r_[0.,internal];east=np.r_[internal,2*math.pi*R*h]
    rhs=vol*cap/dt*old;rhs[-1]+=east[-1]*external
    return thomas(-west,vol*cap/dt+west+east,-east,rhs)

def solve():
    dr=R/N;T=np.full(N,T0);C=np.full(N,C0);times=[0.];TH=[T.copy()];CH=[C.copy()]
    previous_center=C[0];dry_time=None
    for step in range(1,int(T_LIMIT/DT)+1):
        t=step*DT;Ta,Ca=ambient(t)
        cap=np.array([rhofun(c)*cpfun(c) for c in C]);k=np.array([kfun(c) for c in C])
        T=fvm(T,cap,k,H,Ta+273.15,DT,dr)
        D=np.array([dfun(c,temp) for c,temp in zip(C,T)]);C=fvm(C,1.,D,HM,Ca,DT,dr)
        if previous_center>THRESHOLD>=C[0]:
            fraction=(previous_center-THRESHOLD)/(previous_center-C[0]);dry_time=(t-DT)+fraction*DT
        previous_center=C[0]
        if step%2==0:times.append(t);TH.append(T.copy());CH.append(C.copy())
        if dry_time is not None and step%2==0:break
    r=(np.arange(N)+.5)*dr*100
    return np.array(times),r,np.array(TH)-273.15,np.array(CH),dry_time

def sample(v,r,x):return np.interp(x,r,v,left=v[0],right=v[-1])
def save(fig,name):fig.tight_layout();fig.savefig(OUT/name,bbox_inches='tight');plt.close(fig)

def write_result(t,r,C):
    x=np.arange(0,2.01,.1)
    with (OUT/'result3_moisture.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.writer(f);w.writerow(['time_s']+[f'r={z:.1f}cm' for z in x])
        for tt,row in zip(t,C):w.writerow([f'{tt:.0f}']+[f'{v:.6f}' for v in sample(row,r,x)])

def plots(t,r,C,dry_time):
    fig,ax=plt.subplots(figsize=(5.8,3.6));
    for hour,color in zip([6,18,36,54],['#1d6996','#38a6a5','#f6a21a','#ed553b']):
        j=np.argmin(abs(t-hour*3600));ax.plot(r,C[j],lw=2,label=f'{hour} h',color=color)
    ax.set(xlabel='到中心距离 r (cm)',ylabel='含水率 (kg/kg)');ax.grid(alpha=.22);ax.legend(frameon=False);save(fig,'q3_profiles.png')
    fig,ax=plt.subplots(figsize=(6.2,3.6));im=ax.imshow(C,origin='lower',aspect='auto',extent=[0,2,0,t[-1]/3600],cmap='turbo_r');ax.set(xlabel='到中心距离 r (cm)',ylabel='时间 (h)');fig.colorbar(im,ax=ax,label='含水率 (kg/kg)');save(fig,'q3_field.png')
    fig,ax=plt.subplots(figsize=(5.9,3.6));hours=t/3600;ax.plot(hours,C[:,0],lw=2.3,label='中心');ax.plot(hours,C[:,-1],lw=2.3,label='表面');ax.axhline(.15,color='#ed553b',ls='--',label='达标阈值 0.15');ax.axvline(dry_time/3600,color='#555',ls=':');ax.annotate(f'{dry_time/3600:.2f} h',xy=(dry_time/3600,.15),xytext=(-55,35),textcoords='offset points',arrowprops=dict(arrowstyle='->'));ax.set(xlabel='时间 (h)',ylabel='含水率 (kg/kg)');ax.grid(alpha=.22);ax.legend(frameon=False);save(fig,'q3_center_surface.png')
    rate=-np.gradient(C[:,0],t)*3600;fig,ax=plt.subplots(figsize=(5.9,3.6));ax.plot(hours,rate,lw=2.2,color='#1d6996');ax.fill_between(hours,0,rate,color='#38a6a5',alpha=.3);ax.set(xlabel='时间 (h)',ylabel='中心失水速率 (kg/(kg·h))');ax.grid(alpha=.22);save(fig,'q3_dryingrate.png')

def main():
    t,r,T,C,dry=solve();write_result(t,r,C)
    pos=np.arange(0.0,2.01,0.1); header=["time_s"]+[f"r={x:.1f}cm" for x in pos]
    rows=[[float(tt)]+list(sample(row,r,pos)) for tt,row in zip(t,C)]
    write_table(ROOT/"result3.xlsx",[("水分浓度",header,rows)]);plots(t,r,C,dry)
    x=[0,.5,1,1.5,2];print(f'问题三计算完成；烘干时间={dry/3600:.4f} h');print('结束时C:',np.round(sample(C[-1],r,x),4));print('输出目录：',OUT)

if __name__=='__main__':main()
