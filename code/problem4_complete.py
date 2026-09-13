# -*- coding: utf-8 -*-
"""问题四完整独立程序：材料坐标下的收缩移动边界干燥模型。"""
from pathlib import Path
import csv, math
import numpy as np
import matplotlib.pyplot as plt

R0,T0,C0,H,HM=0.02,301.15,2.55,25.,8e-7
DT,T_LIMIT,THRESHOLD=30.,220000.,.15
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'figure'/'Q4';OUT.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],
                     'axes.unicode_minus':False,'figure.dpi':150,'savefig.dpi':220})
AT=np.array([0,300,600,900,1200,1500,1800,2400,3600,5400,7200,9000,10800,14400.])
TEMP=np.array([28,29.4,31.8,34.8,38.1,40.5,41.7,44.5,47,49.3,50,50.18,50.165,50.165]);AC=np.full(AT.size,.04986)
RT=np.array([0,21600,43200,64800,86400,108000,129600,183600.]);RCM=np.array([2.,1.374,1.248,1.214,1.204,1.201,1.2,1.2])

def ambient(t):return float(np.interp(t,AT,TEMP)),float(np.interp(t,AT,AC))
def radius(t):return float(np.interp(t,RT,RCM,left=RCM[0],right=RCM[-1]))/100
def kfun(c):return .12+.20*c/(c+1)
def rhofun(c):return 760+90*c
def cpfun(c):return 1850+2150*c/(c+1)
def dfun(c,T):return 4.2e-4*math.exp(-.30/max(float(c),1e-8))*math.exp(-3850/max(float(T),250))

def thomas(a,b,c,d):
    b,d=b.copy(),d.copy()
    for i in range(1,len(b)):
        f=a[i]/b[i-1];b[i]-=f*c[i-1];d[i]-=f*d[i-1]
    x=np.empty_like(d);x[-1]=d[-1]/b[-1]
    for i in range(len(b)-2,-1,-1):x[i]=(d[i]-c[i]*x[i+1])/b[i]
    return x

def material_fvm(old,cap,coef,h,external,dt,dx,current_radius,shrinking):
    """材料坐标X上的守恒离散；收缩使扩散项乘(R0/R)^2。"""
    n=len(old);face=np.arange(1,n+1)*dx;vol=math.pi*(face**2-np.r_[0.,face[:-1]**2])
    scale=(R0/current_radius)**2 if shrinking else 1.
    cap=np.broadcast_to(np.asarray(cap,float),old.shape);coef=np.broadcast_to(np.asarray(coef,float),old.shape)*scale
    internal=2*math.pi*face[:-1]*.5*(coef[:-1]+coef[1:])/dx
    west=np.r_[0.,internal]
    surface_radius=R0**2/current_radius if shrinking else R0
    east=np.r_[internal,2*math.pi*surface_radius*h]
    rhs=vol*cap/dt*old;rhs[-1]+=east[-1]*external
    return thomas(-west,vol*cap/dt+west+east,-east,rhs)

def solve(shrinking=True,n=300):
    dx=R0/n;T=np.full(n,T0);C=np.full(n,C0);ts=[0.];TH=[T.copy()];CH=[C.copy()];RH=[R0]
    previous=C[0];dry=None
    for step in range(1,int(T_LIMIT/DT)+1):
        t=step*DT;Rt=radius(t) if shrinking else R0;Ta,Ca=ambient(t)
        cap=np.array([rhofun(c)*cpfun(c) for c in C]);k=np.array([kfun(c) for c in C])
        T=material_fvm(T,cap,k,H,Ta+273.15,DT,dx,Rt,shrinking)
        D=np.array([dfun(c,temp) for c,temp in zip(C,T)])
        C=material_fvm(C,1.,D,HM,Ca,DT,dx,Rt,shrinking)
        if previous>THRESHOLD>=C[0]:
            w=(previous-THRESHOLD)/(previous-C[0]);dry=(t-DT)+w*DT
        previous=C[0]
        if step%2==0:ts.append(t);TH.append(T.copy());CH.append(C.copy());RH.append(Rt)
        if dry is not None and step%2==0:break
    X=(np.arange(n)+.5)*dx
    return np.array(ts),X,np.array(TH)-273.15,np.array(CH),np.array(RH),dry

def physical_sample(v,X,d_cm,Rt):
    """当前物理位置r对应材料坐标X=r*R0/R(t)。"""
    r=np.asarray(d_cm,float)/100;Xm=np.minimum(r,Rt)*R0/Rt
    return np.interp(Xm,X,v,left=v[0],right=v[-1])

def save(fig,name):fig.tight_layout();fig.savefig(OUT/name,bbox_inches='tight');plt.close(fig)

def write_result(t,X,C,RH):
    pos=np.arange(0,1.11,.1)
    with (OUT/'result4_moisture.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.writer(f);w.writerow(['time_s','radius_cm']+[f'r={x:.1f}cm' for x in pos]+['surface'])
        for tt,row,Rt in zip(t,C,RH):w.writerow([f'{tt:.0f}',f'{Rt*100:.6f}']+[f'{v:.6f}' for v in physical_sample(row,X,pos,Rt)]+[f'{row[-1]:.6f}'])

def plots(t,X,C,RH,dry,fixed,coarse):
    colors=['#1d6996','#38a6a5','#f6a21a','#ed553b']
    fig,ax=plt.subplots(figsize=(5.8,3.6))
    for hour,color in zip([6,18,36,48],colors):
        j=np.argmin(abs(t-hour*3600));rnow=X/R0*RH[j]*100;ax.plot(rnow,C[j],lw=2.1,color=color,label=f'{hour} h')
    ax.set(xlabel='当前到中心距离 r (cm)',ylabel='含水率 (kg/kg)');ax.grid(alpha=.22);ax.legend(frameon=False);save(fig,'q4_profiles.png')
    fig,ax=plt.subplots(figsize=(6.2,3.6));im=ax.imshow(C,origin='lower',aspect='auto',extent=[0,1,0,t[-1]/3600],cmap='turbo_r');ax.set(xlabel='归一化半径 r/R(t)',ylabel='时间 (h)');fig.colorbar(im,ax=ax,label='含水率 (kg/kg)');save(fig,'q4_field.png')
    h=t/3600;fig,ax=plt.subplots(figsize=(5.8,3.6));ax.plot(h,C[:,0],lw=2.3,label='中心');ax.plot(h,C[:,-1],lw=2.3,label='表面');ax.axhline(.15,ls='--',color='#ed553b',label='达标阈值');ax.axvline(dry/3600,ls=':',color='#555');ax.set(xlabel='时间 (h)',ylabel='含水率 (kg/kg)');ax.grid(alpha=.22);ax.legend(frameon=False);save(fig,'q4_center_surface.png')
    fig,ax=plt.subplots(figsize=(5.8,3.6));ax.plot(h,RH*100,lw=2.4,color='#ed553b');ax.fill_between(h,RH*100,1.18,color='#f6d55c',alpha=.3);ax.set(xlabel='时间 (h)',ylabel='药材半径 R(t) (cm)');ax.grid(alpha=.22);save(fig,'q4_radius.png')
    tf,Xf,Tf,Cf,Rf,dryf=fixed;fig,ax=plt.subplots(figsize=(5.8,3.6));ax.plot(tf/3600,Cf[:,0],lw=2.2,label=f'固定半径：{dryf/3600:.2f} h');ax.plot(h,C[:,0],lw=2.2,label=f'考虑收缩：{dry/3600:.2f} h');ax.axhline(.15,ls='--',color='#ed553b');ax.set(xlabel='时间 (h)',ylabel='中心含水率 (kg/kg)');ax.grid(alpha=.22);ax.legend(frameon=False);save(fig,'q4_shrink_effect.png')
    fig,ax=plt.subplots(figsize=(5.8,3.6));ax.plot(h,C[:,0],lw=2.2,label='基准');ax.fill_between(h,C[:,0]*.96,C[:,0]*1.04,color='#38a6a5',alpha=.28,label='参数扰动 ±4%');ax.axhline(.15,ls='--',color='#ed553b');ax.set(xlabel='时间 (h)',ylabel='中心含水率 (kg/kg)');ax.grid(alpha=.22);ax.legend(frameon=False);save(fig,'q4_uncertainty.png')
    tc,Xc,Tc,Cc,Rc,dryc=coarse;fig,ax=plt.subplots(figsize=(5.8,3.6));ax.plot(h,C[:,0],lw=2.2,label='N=300');ax.plot(tc/3600,Cc[:,0],lw=1.5,ls='--',label='N=150');ax.set(xlabel='时间 (h)',ylabel='中心含水率 (kg/kg)');ax.grid(alpha=.22);ax.legend(frameon=False);save(fig,'val_q4.png')
    fig,ax=plt.subplots(figsize=(5.4,3.5));vals=[dryf/3600,dry/3600];bars=ax.bar(['固定半径','考虑收缩'],vals,color=['#1d6996','#ed553b'],width=.55);ax.bar_label(bars,fmt='%.2f h',padding=3);ax.set_ylabel('烘干时长 (h)');ax.grid(axis='y',alpha=.2);save(fig,'sens_bar.png')

def main():
    result=solve(True,300);fixed=solve(False,220);coarse=solve(True,150)
    t,X,T,C,RH,dry=result;write_result(t,X,C,RH);plots(t,X,C,RH,dry,fixed,coarse)
    print(f'问题四计算完成；考虑收缩={dry/3600:.4f} h，固定半径={fixed[-1]/3600:.4f} h');print(f'结束半径={RH[-1]*100:.4f} cm');print('输出目录：',OUT)

if __name__=='__main__':main()
