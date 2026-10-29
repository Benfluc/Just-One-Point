import numpy as np, json
from scipy import stats
rng=np.random.default_rng(7)
def Dinf(N): return (N-1)/N - stats.norm.cdf(-1/np.sqrt(N))
def Lnull(N,B):
    out=[]
    for _ in range(B//100000):
        z=rng.standard_normal((100000,N))
        u=np.sort(stats.norm.cdf((z-z.mean(1,keepdims=True))/z.std(1,ddof=1,keepdims=True)),1)
        i=np.arange(1,N+1); out.append(np.maximum((i/N-u).max(1),(u-(i-1)/N).max(1)))
    return np.concatenate(out)
tab=[]
for N in range(4,40):
    L=Lnull(N,1_000_000); d=Dinf(N); k=int((L>=d).sum())
    ub=stats.beta.ppf(0.999,k+1,10**6-k)   # Clopper-Pearson unilateral 99,9%
    q=np.quantile(L,[0.90,0.95,0.99])
    tab.append(dict(N=N,Dinf=d,c10=q[0],c05=q[1],c01=q[2],k=k,ub=ub))
    print(N,f"Dinf={d:.4f} c10={q[0]:.4f} c05={q[1]:.4f} c01={q[2]:.4f} P(L>=Dinf)<= {ub:.2e} (k={k})")
json.dump(tab,open("mc_peq.json","w"))
