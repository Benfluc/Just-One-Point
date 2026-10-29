import numpy as np, json
from scipy import stats
from testes import *
g = np.random.default_rng(99)
zs = np.arange(0, 40.0001, 0.02)
out = {}
for N in [10, 25, 50, 100]:
    n = N - 1
    bn, bN = Bateria(n), Bateria(N)
    M = 400
    zL_first, zL_last, zSW, ks_viravel, ks_banda = [], [], [], 0, 0
    far = {k: [0, 0] for k in ["KS", "Lilliefors", "SW", "R", "T"]}
    cnt = 0
    while cnt < M:
        x = g.standard_normal(n)
        rb = bn.rejeita(x)
        # outlier extremo (z = 1e6) para cada teste, entre as amostras que ele não rejeita
        y = np.append(x, x.mean() + 1e6 * x.std(ddof=1)); ry = bN.rejeita(y)
        for k in far:
            if not rb[k]: far[k][0] += 1; far[k][1] += ry[k]
        # KS clássico: existe x que vira?  (máximo de D_{n+1} sobre x)
        if not rb["KS"]:
            xs = np.sort(x); cand = np.concatenate([[-1e9, 1e9], xs - 1e-9, xs + 1e-9, np.linspace(-5, 5, 401)])
            dmax = max(KS(np.append(x, c)) for c in cand)
            ks_viravel += dmax > bN.c["KS"]
            Dn = KS(x); ks_banda += Dn > bN.c["KS"] - (1 - bN.c["KS"]) / n
        if rb["Lilliefors"]: continue
        cnt += 1
        m, s = x.mean(), x.std(ddof=1)
        D = np.array([LILL(np.append(x, m + z * s)) for z in zs]) > bN.c["Lilliefors"]
        zL_first.append(zs[np.argmax(D)] if D.any() else np.inf)
        zL_last.append(zs[np.where(~D)[0].max() + 1] if D[-1] else np.inf)
        if not rb["SW"]:
            for z in np.arange(0, 60, 0.05):
                if stats.shapiro(np.append(x, m + z * s)).pvalue < 0.05: zSW.append(z); break
            else: zSW.append(np.inf)
    q = lambda v: [float(np.quantile(v, p)) for p in (0.05, 0.5, 0.95)]
    out[N] = dict(zL_first=q(zL_first), zL_last=q(zL_last), zL_finito=float(np.mean(np.isfinite(zL_last))),
                  zSW=q(zSW), far={k: v[1] / v[0] for k, v in far.items()}, far_counts={k: v for k, v in far.items()},
                  ks_viravel=None, c=bN.c)
    nks = far["KS"][0]
    out[N]["ks_viravel"] = ks_viravel / nks; out[N]["ks_banda"] = ks_banda / nks; out[N]["ks_counts"] = [ks_viravel, ks_banda, nks]
    print(N, json.dumps(out[N], default=float))
json.dump(out, open("exp_a.json", "w"), default=float)
