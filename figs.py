import numpy as np, json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from testes import *
plt.rcParams.update({"font.family":"DejaVu Serif","font.size":10,"axes.spines.top":False,"axes.spines.right":False,
                     "axes.edgecolor":"#8a8980","axes.grid":True,"grid.color":"#e1e0d9","grid.linewidth":0.6,"figure.dpi":200})
C1,C2,C3,C4="#2a78d6","#eb6834","#1baf7a","#6250d6"; INK="#3b3b38"
g=np.random.default_rng(3)
def Dinf(N): return (N-1)/N-stats.norm.cdf(-1/np.sqrt(N))
def Lnull(N,B):
    z=g.standard_normal((B,N)); u=np.sort(stats.norm.cdf((z-z.mean(1,keepdims=True))/z.std(1,ddof=1,keepdims=True)),1)
    i=np.arange(1,N+1); return np.maximum((i/N-u).max(1),(u-(i-1)/N).max(1))
# ---- dados figura 1
mc=json.load(open("mc_peq.json"))
Ns=[r["N"] for r in mc]; cs={a:[r[k] for r in mc] for a,k in [(0.10,"c10"),(0.05,"c05"),(0.01,"c01")]}
for N in [50,75,100,150,200,300,500,1000]:
    L=Lnull(N,40000 if N<=300 else 15000); Ns.append(N)
    for a in cs: cs[a].append(np.quantile(L,1-a))
Ns=np.array(Ns); Nfine=np.unique(np.round(np.geomspace(4,1000,300)).astype(int))
# ---- exemplos do usuário
A=[1,2,2,2,3]; B=[1.70,1.75,1.80,1.82,1.90,1.88,1.77,1.79,1.65,1.96]
null={N:Lnull(N,200000) for N in [5,6,10,11]}
def pL(x): return (null[len(x)]>=LILL(x)-1e-12).mean()
ex=[]
for nome,base,ad in [("A",A,[None,1,2,3,4,5,6,7,8,9,10]),("B",B,[None,2,2.1,2.2,2.3,2.4,2.5,2.6])]:
    for a in ad:
        x=np.array(base+([] if a is None else [a]),float)
        ex.append(dict(amostra=nome,add=a,D=LILL(x),pL=pL(x),pSW=stats.shapiro(x).pvalue,
                       media=x.mean(),dp=x.std(ddof=1)))
cA=np.quantile(null[6],0.95); cB=np.quantile(null[11],0.95)
def zstar(base,c):
    b=np.array(base,float); m,s=b.mean(),b.std(ddof=1); zs=np.arange(0,30,0.001)
    D=np.array([LILL(np.append(b,m+z*s)) for z in zs])>c
    last_ok=np.where(~D)[0].max(); return zs[last_ok+1], m+zs[last_ok+1]*s, m, s
zA=zstar(A,cA); zB=zstar(B,cB)
# ---- exemplo de aglomerado
xa=np.array([-2,-1.5,-1.5]+[-1]*7+[-0.5]*3+[0]*4+[0.5]*4+[1,2,2],float)
R24,R25=R_test(24),R_test(25); T24,T25=T_test(24),T_test(25)
c24=np.quantile(R24.null,.95); c25=np.quantile(R25.null,.95)
agl=[]
for nome,y,Rt,cc,Tt in [("orig",xa,R24,c24,T24),("moda",np.append(xa,-1.0),R25,c25,T25),("perto",np.append(xa,-0.99),R25,c25,T25)]:
    agl.append(dict(caso=nome,R=Rt.stat(y),cR=cc,pR=(Rt.null>=Rt.stat(y)-1e-12).mean(),
                    T=Tt.stat(y),cT=np.quantile(Tt.null,.95),pT=(Tt.null>=Tt.stat(y)-1e-12).mean(),
                    L=LILL(y),pSW=stats.shapiro(y).pvalue))
json.dump(dict(ex=ex,cA=cA,cB=cB,zA=zA,zB=zB,agl=agl,c_table={str(int(N)):{str(a):float(cs[a][i]) for a in cs} for i,N in enumerate(Ns)}),
          open("exemplos.json","w"),default=float,indent=1)
print(json.dumps(dict(cA=cA,cB=cB,zA=zA,zB=zB),default=float)); [print(e) for e in ex]; [print(a) for a in agl]

LAB={"pt":dict(f1x="Tamanho da amostra ampliada, $N=n+1$",f1y="Valor da estatística",dinf="$D_\\infty(n)$ (limite com $z\\to\\infty$)",
               crit="valor crítico $c_{N,\\alpha}$",f2x="Distância padronizada do elemento acrescentado, $z=(x-\\bar{x}_n)/s_n$",
               f2y="$D_{n+1}$ de Lilliefors",sA="Amostra A ($n=5$)",sB="Amostra B ($n=10$)",crit5="valor crítico (5%)",
               lim="limite $D_\\infty$",rej="região de rejeição",f3x="$x$",f3y="Probabilidade acumulada",fe="FDE com a cópia de $-1$",
               fn="normal ajustada (estimadores de $R$)",reg="região de avaliação $I(K)$",salto="salto de $F_n$ em $x=-1$"),
     "en":dict(f1x="Enlarged sample size, $N=n+1$",f1y="Statistic value",dinf="$D_\\infty(n)$ (limit as $z\\to\\infty$)",
               crit="critical value $c_{N,\\alpha}$",f2x="Standardized distance of the added element, $z=(x-\\bar{x}_n)/s_n$",
               f2y="Lilliefors $D_{n+1}$",sA="Sample A ($n=5$)",sB="Sample B ($n=10$)",crit5="critical value (5%)",
               lim="limit $D_\\infty$",rej="rejection region",f3x="$x$",f3y="Cumulative probability",fe="ECDF with the copy of $-1$",
               fn="fitted normal ($R$ estimators)",reg="evaluation region $I(K)$",salto="jump of $F_n$ at $x=-1$")}
from matplotlib.ticker import FuncFormatter
def virg(ax,lang):
    if lang!="pt": return
    f=FuncFormatter(lambda v,p: (f"{v:g}").replace(".",","))
    ax.yaxis.set_major_formatter(f)
    if ax.get_xscale()!="log": ax.xaxis.set_major_formatter(f)
for lang,T in LAB.items():
    # Figura 1
    fig,ax=plt.subplots(figsize=(6.3,3.6))
    ax.plot(Nfine,[Dinf(N) for N in Nfine],color=INK,lw=2,label=T["dinf"])
    for (a,col,ls) in [(0.10,C3,"-"),(0.05,C1,"-"),(0.01,C4,"-")]:
        ax.plot(Ns,cs[a],color=col,lw=2,ls=ls,label=f'{T["crit"]}, $\\alpha={a:.2f}$'.replace(".","{,}") if lang=="pt" else f'{T["crit"]}, $\\alpha={a:.2f}$')
    ax.set_xscale("log"); ax.set_xlabel(T["f1x"]); ax.set_ylabel(T["f1y"]); ax.set_ylim(0,0.6)
    ax.set_xticks([4,10,30,100,300,1000]); ax.set_xticklabels(["4","10","30","100","300","1000"])
    virg(ax,lang); ax.legend(frameon=False,fontsize=8.5,loc="lower left"); fig.tight_layout(); fig.savefig(f"fig1_{lang}.png"); plt.close(fig)
    # Figura 2
    fig,axs=plt.subplots(1,2,figsize=(6.6,3.2),sharey=True)
    for ax,base,c,tit in [(axs[0],A,cA,T["sA"]),(axs[1],B,cB,T["sB"])]:
        b=np.array(base,float); m,s=b.mean(),b.std(ddof=1); zs=np.linspace(0,15,1501)
        D=[LILL(np.append(b,m+z*s)) for z in zs]
        ax.axhspan(c,0.6,color=C2,alpha=0.08,lw=0); ax.text(14.8,0.585,T["rej"],ha="right",va="top",fontsize=8,color=C2)
        ax.plot(zs,D,color=C1,lw=2); ax.axhline(c,color=C2,lw=1.2,ls="--",label=T["crit5"])
        ax.axhline(Dinf(len(b)+1),color=INK,lw=1,ls=":",label=T["lim"])
        ax.set_title(tit,fontsize=10); ax.set_xlabel("$z$"); ax.set_ylim(0,0.6)
    [virg(a,lang) for a in axs]; axs[0].set_ylabel(T["f2y"]); axs[1].legend(frameon=False,fontsize=8,loc="lower right")
    fig.tight_layout(); fig.savefig(f"fig2_{lang}.png"); plt.close(fig)
    # Figura 3
    y=np.sort(np.append(xa,-1.0)); N=len(y); a=janela_sd(y,R25.m); K=y[a:a+R25.m]
    mu,sg=K.mean(),R25.c_sigma*K.std(ddof=1); lo,hi=y[a],y[a+R25.m-1]
    fig,ax=plt.subplots(figsize=(6.3,3.6))
    ax.axvspan(lo,hi,color=C3,alpha=0.08,lw=0,label=T["reg"])
    u=np.unique(y); Fv=np.searchsorted(y,u,side="right")/N
    ax.step(np.concatenate([[-2.6],u,[2.6]]),np.concatenate([[0],Fv,[1]]),where="post",color=C1,lw=2,label=T["fe"])
    xx=np.linspace(-2.6,2.6,400); ax.plot(xx,stats.norm.cdf((xx-mu)/sg),color=C2,lw=2,label=T["fn"])
    Fl=np.searchsorted(y,-1.0,side="left")/N; Fr=np.searchsorted(y,-1.0,side="right")/N
    ax.annotate("",xy=(-1.0,Fr),xytext=(-1.0,Fl),arrowprops=dict(arrowstyle="<->",color=INK,lw=1.2))
    ax.text(-0.93,(Fl+Fr)/2,T["salto"]+"\n$J=8/25$",fontsize=8.5,color=INK,va="center")
    ax.set_xlabel(T["f3x"]); ax.set_ylabel(T["f3y"]); ax.set_xlim(-2.6,2.6); ax.set_ylim(0,1.02)
    virg(ax,lang); ax.legend(frameon=False,fontsize=8.5,loc="lower right"); fig.tight_layout(); fig.savefig(f"fig3_{lang}.png"); plt.close(fig)
