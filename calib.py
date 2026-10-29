"""Valores críticos com erro de Monte Carlo e tamanho empírico dos cinco testes."""
import numpy as np, json
from scipy import stats
from testes import Bateria, KS
g = np.random.default_rng(3000); out = {}
for N in [9, 10, 24, 25, 49, 50, 99, 100]:
    b = Bateria(N); row = {}
    for k in ["Lilliefors", "R", "T"]:
        v = b.nulls[k]; B = len(v); vs = np.sort(v)
        # IC 95% de distribuição livre para o quantil 0,95 (estatísticas de ordem)
        lo = int(np.floor(B * 0.95 - 1.96 * np.sqrt(B * 0.95 * 0.05))); hi = int(np.ceil(B * 0.95 + 1.96 * np.sqrt(B * 0.95 * 0.05)))
        row[k] = dict(c=float(b.c[k]), lo=float(vs[lo]), hi=float(vs[hi]), B=B)
    row["KS"] = dict(c=float(b.c["KS"]), exato=True)
    rej = {k: 0 for k in ["KS", "Lilliefors", "SW", "R", "T"]}; M = 5000
    for _ in range(M):
        for k, v in b.rejeita(g.standard_normal(N)).items(): rej[k] += int(v)
    row["tamanho"] = {k: v / M for k, v in rej.items()}; row["M"] = M
    row["c_sigma"] = dict(R=float(b.R.c_sigma), T=float(b.T.c_sigma), i0=b.T.i0)
    out[N] = row; print(N, json.dumps(row))
json.dump(out, open("calib.json", "w"), indent=1)
