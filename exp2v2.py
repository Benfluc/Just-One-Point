"""Experimento 2 revisto: mesmas amostras para os cinco testes, calibração nominal e sob arredondamento,
intervalos de Wilson, taxas incondicionais e condicionais."""
import numpy as np, json, sys
from scipy import stats
from testes import KS, LILL, Bateria, R_test, T_test
N = int(sys.argv[1]); SEED = 1000 + N
g = np.random.default_rng(SEED); n = N - 1
TESTES = ["KS", "Lilliefors", "SW", "R", "T"]
bn, bN = Bateria(n), Bateria(N)
def estat(x, b):
    return dict(KS=KS(x), Lilliefors=LILL(x), SW=stats.shapiro(x).statistic, SWp=stats.shapiro(x).pvalue,
                R=b.R.stat(x), T=b.T.stat(x))
def rej_nom(s, b):
    return dict(KS=s["KS"] > b.c["KS"], Lilliefors=s["Lilliefors"] > b.c["Lilliefors"], SW=s["SWp"] < 0.05,
                R=s["R"] > b.c["R"], T=s["T"] > b.c["T"])
def calib(h, size, b, B=20000):
    v = {k: [] for k in ["KS", "Lilliefors", "SW", "R", "T"]}
    for _ in range(B):
        s = estat(np.round(g.standard_normal(size) / h) * h, b)
        for k in v: v[k].append(s[k])
    c = {k: float(np.quantile(v[k], 0.95)) for k in ["KS", "Lilliefors", "R", "T"]}
    c["SW"] = float(np.quantile(v["SW"], 0.05))
    return c
def rej_cal(s, c):
    return dict(KS=s["KS"] > c["KS"], Lilliefors=s["Lilliefors"] > c["Lilliefors"], SW=s["SW"] < c["SW"],
                R=s["R"] > c["R"], T=s["T"] > c["T"])
out = {}
for h in [0.25, 0.5, 1.0]:
    cn, cN = calib(h, n, bn), calib(h, N, bN)
    M = 2000; rec = []
    for _ in range(M):
        x = np.round(g.standard_normal(n) / h) * h
        vals, cnt = np.unique(x, return_counts=True); moda = vals[np.argmax(cnt)]
        adds = dict(moda=moda, perto=moda + h / 100, nova=g.standard_normal())
        s0 = estat(x, bn); sa = {a: estat(np.append(x, v), bN) for a, v in adds.items()}
        r = dict(nom0=rej_nom(s0, bn), cal0=rej_cal(s0, cn),
                 nom={a: rej_nom(sa[a], bN) for a in adds}, cal={a: rej_cal(sa[a], cN) for a in adds},
                 moda_med=bool(moda == np.median(np.append(x, moda))))
        rec.append(r)
    res = {}
    for cal in ["nom", "cal"]:
        base = {k: float(np.mean([r[cal + "0"][k] for r in rec])) for k in TESTES}
        depois = {a: {k: float(np.mean([r[cal][a][k] for r in rec])) for k in TESTES} for a in ["moda", "perto", "nova"]}
        cond = {}
        for k in TESTES:
            nr = [r for r in rec if not r[cal + "0"][k]]
            cond[k] = dict(nr=len(nr), **{a: int(sum(r[cal][a][k] for r in nr)) for a in ["moda", "perto", "nova"]})
        comum = [r for r in rec if not any(r[cal + "0"][k] for k in TESTES)]
        cond_comum = dict(nr=len(comum), **{k: {a: int(sum(r[cal][a][k] for r in comum)) for a in ["moda", "perto", "nova"]} for k in TESTES})
        # T: separação por moda == mediana
        nrT = [r for r in rec if not r[cal + "0"]["T"]]
        tsep = {str(flag): dict(nr=len([r for r in nrT if r["moda_med"] == flag]),
                                moda=int(sum(r[cal]["moda"]["T"] for r in nrT if r["moda_med"] == flag))) for flag in [True, False]}
        res[cal] = dict(base=base, depois=depois, cond=cond, comum=cond_comum, T_moda_mediana=tsep)
    res["crit_cal"] = dict(n=cn, N=cN)
    out[str(h)] = res
    print(N, h, json.dumps({c: {k: (res[c]["base"][k], res[c]["cond"][k]) for k in TESTES} for c in ["nom", "cal"]}), flush=True)
json.dump(dict(N=N, seed=SEED, M=2000, B_calib=20000, res=out), open(f"exp2v2_{N}.json", "w"), indent=1)
