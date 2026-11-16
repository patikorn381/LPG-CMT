"""Figure 4 of the paper, drawn from results/fom_map.csv (written by fom_map.m).
Run from this folder:  python make_fig8.py   (needs numpy, pandas, matplotlib)"""
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
d=pd.read_csv('results/fom_map.csv')
cols={5:'#2a6fdb',7:'#e8603c',9:'#1aa37a',11:'#f2a20c',13:'#e86fa6',15:'#1b7d1b',17:'#6a4fd8'}
mk={5:'o',7:'s',9:'^',11:'D',13:'v',15:'P',17:'X'}
plt.rcParams.update({'font.size':10,'font.family':'DejaVu Sans'})
fig,ax=plt.subplots(2,2,figsize=(13.2,10))
s=d[d.n3==1.33]
for nu in sorted(cols):
    m=s[s.nu==nu].sort_values('lambda_r_nm')
    for br in (m[m.dn_g>0], m[m.dn_g<0]):
        br=br.dropna(subset=['FOM_sim'])
        if br.empty: continue
        kw=dict(color=cols[nu],lw=1.8)
        lab=f'ν = {nu}' if br is not None and (br.dn_g.iloc[0]>0) else None
        ax[0,0].plot(br.lambda_r_nm,br.S_nm_RIU.abs(),marker=mk[nu],markevery=8,ms=5,label=lab,**kw)
        ax[0,1].plot(br.lambda_r_nm,br.FWHM_nm,marker=mk[nu],markevery=8,ms=5,**kw)
        ok=br.dn_g.abs()>=1.5e-3
        y=br.FOM_sim.where(ok)
        ax[1,0].plot(br.lambda_r_nm,y,marker=mk[nu],markevery=8,ms=5,**kw)
        ax[1,0].plot(br.lambda_r_nm,br.FOM_closedform,ls=':',color=cols[nu],lw=1.5)
for a,t,yl in [(ax[0,0],'(a) Sensitivity','|S| (nm/RIU)'),(ax[0,1],'(b) Bandwidth','FWHM (nm)'),(ax[1,0],'(c) FOM: simulated (solid), Eq. (7) (dotted)','FOM = |S|/FWHM (RIU$^{-1}$)')]:
    a.set_yscale('log'); a.set_xlim(700,1800); a.set_xlabel('Resonant wavelength $\\lambda_r$ (nm)'); a.set_ylabel(yl)
    a.set_title(t,loc='left'); a.grid(alpha=.3)
ax[0,0].legend(ncol=2,frameon=False,loc='upper left')
v=d.dropna(subset=['FOM_sim','FOM_closedform'])
near=v[v.dn_g.abs()<1.5e-3]; far=v[v.dn_g.abs()>=1.5e-3]
a=ax[1,1]
a.plot([1e-2,20],[1e-2,20],color='gray',lw=1)
a.scatter(near.FOM_closedform,near.FOM_sim,s=18,facecolors='none',edgecolors='#e8603c',label='|Δn$_g$| < 1.5×10$^{-3}$ (near TAP)')
a.scatter(far.FOM_closedform,far.FOM_sim,s=14,color='#2a6fdb',label='|Δn$_g$| ≥ 1.5×10$^{-3}$')
a.set_xscale('log'); a.set_yscale('log'); a.set_xlim(1e-2,20); a.set_ylim(1e-2,20)
a.set_xlabel('FOM, Eq. (7) (RIU$^{-1}$)'); a.set_ylabel('FOM, simulated (RIU$^{-1}$)')
a.set_title('(d) Eq. (7) vs simulation, SRI 1.00/1.33/1.40',loc='left'); a.grid(alpha=.3); a.legend(frameon=False,loc='upper left')
from matplotlib.ticker import FuncFormatter
import math
f=FuncFormatter(lambda v,_: f'{v:g}' if round(v/10**math.floor(math.log10(v)+1e-9),3) in (1,2,3,5) else '')
ax[1,0].yaxis.set_major_formatter(f); ax[1,0].yaxis.set_minor_formatter(f)
plt.tight_layout(); plt.savefig('results/fig8_fom_map.png',dpi=100)
print('designs:',len(v),' |dn_g|>=1.5e-3:',len(far),' near TAP:',len(near))
