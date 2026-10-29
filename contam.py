"""Experimento 3: contaminação probabilística em amostras normais contínuas (mesmas amostras para todos os testes)."""
import numpy as np, json, sys
from scipy import stats
from testes import KS, LILL, Bateria
N = int(sys.argv[1]); SEED = 2000 + N; g = np.random.default_rng(SEED); n = N - 1
TESTES = ["KS", "Lilliefors", "SW", "R", "T"]
bn, bN = Bateria(n), Bateria(N)
fontes = {"N(0,1)": lambda: g.standard_normal(), "N(0,3^2)": lambda: 3 * g.standard_normal(),
          "N(0,5^2)": lambda: 5 * g.standard_normal(), "t_3": lambda: g.standard_t(3)}
M = 3000; rec = []
for _ in range(M):
    x = g.standard_normal(n); r0 = bn.rejeita(x)
    ra = {f: bN.rejeita(np.append(x, fn())) for f, fn in fontes.items()}
    rec.append((r0, ra))
res = {"cond": {}, "comum": {}}
for k in TESTES:
    nr = [r for r in rec if not r[0][k]]
    res["cond"][k] = dict(nr=len(nr), **{f: int(sum(r[1][f][k] for r in nr)) for f in fontes})
comum = [r for r in rec if not any(r[0].values())]
res["comum"] = dict(nr=len(comum), **{k: {f: int(sum(r[1][f][k] for r in comum)) for f in fontes} for k in TESTES})
json.dump(dict(N=N, seed=SEED, M=M, res=res), open(f"contam_{N}.json", "w"), indent=1)
print(N, json.dumps(res["cond"]))
