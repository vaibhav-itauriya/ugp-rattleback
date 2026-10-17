"""Figures 1-7 of the paper, from the runs written by run_simulations.py.

Usage:  python3 figures.py            (all figures -> figs/)
        python3 figures.py eps rho    (selected figures only)
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from model import geom, rolling_wd, contact, m, g, I, a1, a2, a3, b
OUT=os.path.join(HERE, 'figs') + os.sep
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.family':'serif','font.serif':['STIXGeneral','DejaVu Serif'],
 'mathtext.fontset':'stix','font.size':8,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,
 'legend.fontsize':6.5,'axes.linewidth':0.6,'xtick.major.width':0.5,'ytick.major.width':0.5,
 'xtick.major.size':2.5,'ytick.major.size':2.5,'xtick.direction':'in','ytick.direction':'in',
 'xtick.top':True,'ytick.right':True,'lines.linewidth':0.9,'legend.frameon':False,
 'savefig.dpi':400,'savefig.bbox':'tight','savefig.pad_inches':0.02})
COL=3.40; WIDE=7.0
K={'roll':'#000000','coul':'#0072B2','spin':'#D55E00','visc':'#CC79A7','w1':'#0072B2','w2':'#D55E00','w3':'#009E73','grey':'#8a8a8a'}
def L(p): D=np.load(p); return D['t'],D['Y']
def lab(ax,s,x=0.015,y=0.96): ax.text(x,y,s,transform=ax.transAxes,va='top',ha='left',fontsize=8)
C=lambda mu=0.15,eps=3e-5,rho=0.0,**k: dict(law='coulomb',mu=mu,eps=eps,rho=rho,**k)
def Nser(t,Y,P,step=1):
    out=[]
    for k in range(0,len(t),step):
        w,ang=Y[k,0:3],Y[k,3:6]
        if P['law']=='rolling':
            n,nd,r,rd,h=geom(ang,w); wd=rolling_wd(w,ang)
            aG=-np.cross(wd,r)-np.cross(w,np.cross(w,r))-np.cross(w,rd); out.append((m*aG+m*g*n)@n)
        else: out.append(contact(w,ang,Y[k,6:9],P)[2])
    return np.array(out)
def slip_energy(t,Y,P,step=2):
    S,E=[],[]
    for k in range(0,len(t),step):
        w,ang=Y[k,0:3],Y[k,3:6]; n,nd,r,rd,h=geom(ang,w)
        if P['law']=='rolling': v=-np.cross(w,r); S.append(0.0)
        else:
            v=Y[k,6:9]; vp=v+np.cross(w,r); vp=vp-(vp@n)*n; S.append(np.linalg.norm(vp))
        E.append(0.5*w@(I@w)+0.5*m*(v@v)+m*g*h)
    return t[::step],np.array(S),np.array(E)
def rev(t,w3,w30=-5):
    i=np.where((np.sign(w3)==-np.sign(w30))&(t>0.05))[0]; return t[i[0]] if len(i) else np.nan
def env(y,win):
    from scipy.ndimage import maximum_filter1d, minimum_filter1d
    return minimum_filter1d(y,win), maximum_filter1d(y,win)

R1 = R2 = os.path.join(HERE, 'runs') + os.sep
roll_n=L(R1+'rolling_nonpref.npz'); roll_p=L(R1+'rolling_pref.npz')
coul_n=L(R2+'coulomb_non.npz'); spin_n=L(R2+'spin05_non.npz'); visc_n=L(R1+'viscous_nonpref.npz')
PR=dict(law='rolling'); PV=dict(law='viscous',cv=500.0)

# ---------------- Fig 1: schematic -------------------------------------------
def fig_schematic():
    fig,ax=plt.subplots(figsize=(COL,1.85)); ax.set_aspect('equal'); ax.axis('off')
    th=np.radians(8.0); A,Bv,bb=1.65,0.62,0.24          # drawn proportions, not to scale
    u=np.linspace(np.pi,2*np.pi,240); Rm=np.array([[np.cos(th),-np.sin(th)],[np.sin(th),np.cos(th)]])
    pts=Rm@np.vstack([A*np.cos(u),Bv*np.sin(u)]); k=np.argmin(pts[1]); off=-pts[:,k:k+1]; pts=pts+off
    cap=Rm@np.vstack([[-A,A],[0,0]])+off
    ax.fill(np.r_[pts[0],cap[0][::-1]],np.r_[pts[1],cap[1][::-1]],color='#ecd9c2',ec='k',lw=0.7)
    ax.plot([-2.35,2.35],[0,0],color='k',lw=0.8)
    for x in np.arange(-2.3,2.35,0.17): ax.plot([x,x-0.11],[0,-0.11],color='k',lw=0.35)
    O=off[:,0].copy(); Gp=(Rm@np.array([0,-bb]))+off[:,0]; P=np.array([0.0,0.0])
    def arr(p,q,c='k',lw=0.8,ms=6): ax.add_patch(FancyArrowPatch(tuple(p),tuple(q),arrowstyle='-|>',mutation_scale=ms,color=c,lw=lw,shrinkA=0,shrinkB=0))
    e3=Rm@np.array([0,1.0]); e1=Rm@np.array([1.0,0])
    arr(Gp,Gp+0.70*e3); ax.text(*(Gp+0.74*e3+np.array([0.05,0.0])),r'$\mathbf{e}_3$',fontsize=8)
    arr(Gp,Gp+0.80*e1); ax.text(*(Gp+0.84*e1+np.array([0.03,-0.02])),r'$\mathbf{e}_1$',fontsize=8)
    arr((-2.1,0.06),(-2.1,0.85)); ax.text(-2.03,0.78,r'$\mathbf{E}_3$',fontsize=8)
    arr(Gp,P,c=K['coul'],lw=0.9); m_=(Gp+P)/2; ax.text(m_[0]+0.08,m_[1]-0.04,r'$\mathbf{r}_{P/G}$',fontsize=8,color=K['coul'])
    arr(P,(P[0],0.62),c=K['spin'],lw=0.9)
    ax.text(P[0]-0.34,0.50,r'$N\mathbf{n}$',fontsize=8,color=K['spin'])
    arr((P[0],-0.035),(P[0]-0.75,-0.035),c=K['spin'],lw=0.9); ax.text(P[0]-0.72,-0.24,r'$\mathbf{F}$',fontsize=8,color=K['spin'])
    ax.plot(*P,'o',ms=2.8,color='k'); ax.text(P[0]+0.06,-0.24,'$P$',fontsize=8)
    ax.plot(*O,'o',ms=2.6,color='k'); ax.text(O[0]-0.2,O[1]+0.05,'$O$',fontsize=8)
    ax.plot(*Gp,'o',ms=3.4,mfc='w',mec='k',mew=0.8); ax.text(Gp[0]-0.25,Gp[1]-0.12,'$G$',fontsize=8)
    ax.text((O[0]+Gp[0])/2+0.06,(O[1]+Gp[1])/2-0.02,'$b$',fontsize=7.5)
    ax.set_xlim(-2.35,2.35); ax.set_ylim(-0.32,1.2)
    fig.savefig(OUT+'fig1_schematic.pdf'); plt.close(fig)

# ---------------- Fig 3: rolling solution ------------------------------------
def fig_rolling():
    t,Y=roll_n; tp,Yp=roll_p
    fig,axs=plt.subplots(3,1,figsize=(COL,3.6),sharex=True,gridspec_kw=dict(hspace=0.08))
    ax=axs[0]
    ax.plot(t,Y[:,0],color=K['w1'],lw=0.45,label=r'$\omega_1$'); ax.plot(t,Y[:,1],color=K['w2'],lw=0.45,label=r'$\omega_2$')
    ax.plot(t,Y[:,2],color=K['w3'],lw=1.2,label=r'$\omega_3$'); ax.plot(tp,Yp[:,2],color=K['grey'],lw=0.9,ls='--',label=r'$\omega_3$, preferred')
    ax.set_ylabel(r'$\omega_i$ (rad s$^{-1}$)'); ax.set_ylim(-8,8); ax.legend(ncol=4,loc='lower left',bbox_to_anchor=(0.0,1.0),handlelength=1.4,columnspacing=0.9,borderaxespad=0.2)
    lab(ax,'(a)',y=0.13)
    ax=axs[1]; ax.plot(t,np.degrees(Y[:,4]),color='k',lw=0.5); ax.set_ylabel(r'$\theta$ (deg)'); lab(ax,'(b)')
    ax=axs[2]; sel=(t>0)&(t<8); Nn=Nser(t,Y,PR); ax.plot(t,Nn/(m*g),color='k',lw=0.45); ax.axhline(1,color=K['grey'],lw=0.6,ls='--')
    ax.set_ylabel(r'$N/mg$'); ax.set_xlabel(r'$t$ (s)'); ax.set_xlim(0,8); lab(ax,'(c)')
    fig.savefig(OUT+'fig3_rolling.pdf'); plt.close(fig)
    return Nn.min()/(m*g), Nn.max()/(m*g)

# ---------------- Fig 4: fixed normal reaction -------------------------------
def fig_fixN():
    fig,axs=plt.subplots(2,1,figsize=(COL,3.0),gridspec_kw=dict(hspace=0.42))
    ax=axs[0]; t,Y=roll_n
    ax.plot(t,Y[:,2],color='#bdbdbd',lw=3.2,label='rolling',solid_capstyle='butt')
    t1,Y1=L(R2+'mu_1.0.npz'); ax.plot(t1,Y1[:,2],color=K['coul'],lw=0.8,ls='-',label=r'slipping, $N$ solved')
    t2,Y2=L(R2+'fixN_mu_1.0.npz'); ax.plot(t2,Y2[:,2],color=K['spin'],lw=1.0,ls='--',label=r'slipping, $N\equiv mg$')
    ax.set_xlim(0,6); ax.set_ylim(-5.5,5.5); ax.set_xlabel(r'$t$ (s)'); ax.set_ylabel(r'$\omega_3$ (rad s$^{-1}$)')
    ax.legend(loc='lower right',handlelength=2.0); lab(ax,'(a)',x=0.02,y=0.95)
    ax=axs[1]; mus=np.array([0.15,0.2,0.5,1.0,2.0]); te=[]; tf=[]
    for mu in mus:
        te.append(rev(*(lambda d:(d[0],d[1][:,2]))(L(R2+('coulomb_non' if mu==0.15 else f'mu_{mu}')+'.npz'))))
        tf.append(rev(*(lambda d:(d[0],d[1][:,2]))(L(R2+f'fixN_mu_{mu}.npz'))))
    tr=rev(t,Y[:,2])
    ax.axhline(tr,color=K['roll'],lw=0.7,ls=':'); ax.text(1.2,tr+0.06,'rolling',fontsize=7)
    ax.plot(mus,te,'o-',color=K['coul'],ms=3.2,mfc='w',label=r'$N$ solved')
    ax.plot(mus,tf,'s--',color=K['spin'],ms=3.0,mfc='w',label=r'$N\equiv mg$')
    ax.set_xscale('log'); ax.set_xticks([0.15,0.2,0.5,1,2]); ax.set_xticklabels(['0.15','0.2','0.5','1','2']); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel(r'$\mu$'); ax.set_ylabel(r'$t_{\rm rev}$ (s)'); ax.legend(loc='upper right'); lab(ax,'(b)',x=0.40,y=0.95)
    fig.savefig(OUT+'fig4_fixN.pdf'); plt.close(fig)
    return dict(zip(mus,zip(te,tf))), tr

# ---------------- Fig 5: four contact laws -----------------------------------
def fig_compare():
    runs=[('rolling',roll_n,PR,K['roll']),(r'Coulomb, $\mu=0.15$',coul_n,C(),K['coul']),
          (r'Coulomb + drilling, $\rho=0.5$ mm',spin_n,C(rho=0.5e-3),K['spin']),(r'viscous, $c=500$ kg s$^{-1}$',visc_n,PV,K['visc'])]
    fig,axs=plt.subplots(2,2,figsize=(WIDE,3.7),gridspec_kw=dict(hspace=0.30,wspace=0.24))
    out={}
    for name,(t,Y),P,c in runs:
        axs[0,0].plot(t,Y[:,2],color=c,lw=1.0,label=name)
        lo,hi=env(np.degrees(Y[:,4]),int(0.12*(len(t)/t[-1])))
        axs[0,1].plot(t,hi,color=c,lw=0.9)
        ts,S,E=slip_energy(t,Y,P,step=4)
        if P['law']!='rolling': axs[1,0].plot(ts,S*1e3,color=c,lw=0.6)
        axs[1,1].plot(ts,E/E[0],color=c,lw=1.0)
        out[name]=dict(trev=rev(t,Y[:,2]),w3=Y[-1,2],thmax=np.degrees(abs(Y[:,4]).max()),dE=(E[-1]-E[0])/E[0],mono=bool(np.all(np.diff(E)<1e-9)),slipmax=S.max())
    a=axs[0,0]; a.set_ylabel(r'$\omega_3$ (rad s$^{-1}$)'); a.set_ylim(-5.6,5.6); lab(a,'(a)')
    h,l=a.get_legend_handles_labels(); lg=fig.legend(h,l,loc='lower center',ncol=4,bbox_to_anchor=(0.5,0.985),handlelength=2.2,columnspacing=1.6,fontsize=7)
    [ln.set_linewidth(2.2) for ln in lg.get_lines()]
    a=axs[0,1]; a.set_ylabel(r'$\theta$ envelope (deg)'); a.set_yscale('log'); a.set_ylim(1e-3,10); lab(a,'(b)')
    a=axs[1,0]; a.set_ylabel(r'$|\mathbf{v}_P|$ (mm s$^{-1}$)'); a.set_xlabel(r'$t$ (s)'); lab(a,'(c)')
    a=axs[1,1]; a.set_ylabel(r'$E/E_0$'); a.set_xlabel(r'$t$ (s)'); lab(a,'(d)',x=0.93,y=0.95)
    for a in axs.flat: a.set_xlim(0,8)
    fig.savefig(OUT+'fig5_compare.pdf'); plt.close(fig)
    return out

# ---------------- Fig 2: regularisation convergence --------------------------
def fig_eps():
    E=[(1e-2,R1+'eps_0.01'),(3e-3,R1+'eps_0.003'),(1e-3,R1+'eps_0.001'),(3e-4,R1+'eps_0.0003'),
       (1e-4,R1+'eps_0.0001'),(3e-5,R2+'coulomb_non'),(1e-5,R2+'eps_1e-05')]
    eps,tr,w3=[],[],[]
    for e,f in E:
        t,Y=L(f+'.npz'); eps.append(e); tr.append(rev(t,Y[:,2])); w3.append(Y[-1,2])
    eps=np.array(eps); tr=np.array(tr); w3=np.array(w3)
    fig,axs=plt.subplots(2,1,figsize=(COL,2.7),sharex=True,gridspec_kw=dict(hspace=0.1))
    axs[0].plot(eps,tr,'o-',color=K['coul'],ms=3.2,mfc='w'); axs[0].set_ylabel(r'$t_{\rm rev}$ (s)'); lab(axs[0],'(a)',x=0.9)
    axs[1].plot(eps,w3,'s-',color=K['coul'],ms=3.0,mfc='w'); axs[1].axhline(0,color=K['grey'],lw=0.5,ls=':')
    axs[1].set_ylabel(r'$\omega_3(8\,{\rm s})$ (rad s$^{-1}$)'); axs[1].set_xlabel(r'$\varepsilon$ (m s$^{-1}$)'); axs[1].set_xscale('log'); lab(axs[1],'(b)',x=0.9)
    for a in axs: a.axvspan(1e-5*0.8,3e-5*1.25,color='#e8eef5',zorder=0)
    fig.savefig(OUT+'fig2_eps.pdf'); plt.close(fig)
    return list(zip(eps,tr,w3))

# ---------------- Fig 6: viscous sweep ---------------------------------------
def fig_visc():
    fig,ax=plt.subplots(figsize=(COL,1.9))
    t,Y=roll_n; ax.plot(t,Y[:,2],color=K['roll'],lw=1.0,label='rolling')
    cs=[(50,R1+'cv_50'),(200,R1+'cv_200'),(500,R1+'viscous_nonpref'),(2000,R1+'cv_2000')]
    shades=['#e8c6dc','#d99bc4','#CC79A7','#7a2f5f']; res={}
    for (c,f),sh in zip(cs,shades):
        t2,Y2=L(f+'.npz'); ax.plot(t2,Y2[:,2],color=sh,lw=1.0,label=fr'$c={c}$'); res[c]=(rev(t2,Y2[:,2]),Y2[-1,2])
    for c,f in [(1.5,R1+'cv_1.5'),(10,R1+'cv_10')]:
        t2,Y2=L(f+'.npz'); res[c]=(rev(t2,Y2[:,2]),Y2[-1,2])
    ax.set_xlim(0,8); ax.set_ylim(-5.6,5.6); ax.set_xlabel(r'$t$ (s)'); ax.set_ylabel(r'$\omega_3$ (rad s$^{-1}$)')
    ax.legend(loc='upper left',handlelength=1.5,ncol=1,bbox_to_anchor=(0.0,1.0))
    fig.savefig(OUT+'fig6_viscous.pdf'); plt.close(fig); return res

# ---------------- Fig 7: drilling sweep --------------------------------------
def fig_rho():
    rhos=[0,0.1,0.2,0.5,1.0,2.0]; wp,wn,th,tr=[],[],[],[]
    for r in rhos:
        if r==0: tn,Yn=L(R2+'coulomb_non.npz'); tp,Yp=L(R2+'coulomb_pref.npz')
        elif r==0.5: tn,Yn=L(R2+'spin05_non.npz'); tp,Yp=L(R2+'rho_0.5_pref.npz')
        else: tn,Yn=L(R2+f'rho_{r}_non.npz'); tp,Yp=L(R2+f'rho_{r}_pref.npz')
        wp.append(Yp[-1,2]); wn.append(Yn[-1,2]); th.append(np.degrees(abs(Yn[:,4]).max())); tr.append(rev(tn,Yn[:,2]))
    fig,axs=plt.subplots(2,1,figsize=(COL,2.8),sharex=True,gridspec_kw=dict(hspace=0.1))
    for a in axs: a.axvspan(0.16,0.52,color='#fbe3d3',zorder=0)
    axs[0].plot(rhos,wp,'o-',color=K['coul'],ms=3.2,mfc='w',label='preferred')
    axs[0].plot(rhos,wn,'s-',color=K['spin'],ms=3.0,mfc='w',label='non-preferred')
    axs[0].set_ylabel(r'$\omega_3(8\,{\rm s})$ (rad s$^{-1}$)'); axs[0].legend(loc='upper right'); lab(axs[0],'(a)',x=0.55,y=0.92)
    axs[1].plot(rhos,th,'^-',color='k',ms=3.2,mfc='w'); axs[1].set_ylabel(r'$\theta_{\max}$ (deg)'); axs[1].set_xlabel(r'$\rho$ (mm)')
    axs[1].text(0.34,0.25,'Hertz\nestimate',ha='center',fontsize=6.5,color='#9a4b17'); lab(axs[1],'(b)',x=0.9,y=0.92)
    fig.savefig(OUT+'fig7_rho.pdf'); plt.close(fig)
    return list(zip(rhos,tr,wn,wp,th))

if __name__=='__main__':
    which=sys.argv[1:] or ['schem','roll','fixN','cmp','eps','visc','rho']
    if 'schem' in which: fig_schematic(); print('schematic ok')
    if 'roll' in which: print('rolling N/mg range', fig_rolling())
    if 'fixN' in which: print('fixN', fig_fixN())
    if 'cmp' in which:
        for k,v in fig_compare().items(): print('cmp',k,v)
    if 'eps' in which: print('eps', fig_eps())
    if 'visc' in which: print('visc', fig_visc())
    if 'rho' in which: print('rho', fig_rho())

