"""Proposição assintótica: zbar / N^{1/4} -> C, com C definido pelo processo-limite de Lilliefors (Durbin, 1973)."""
import numpy as np, json
from scipy import stats
Phi, phi = stats.norm.cdf, stats.norm.pdf
g = np.random.default_rng(4000)
def Lnull(N, B):
    out = []
    for _ in range(max(1, B // max(1, 2_000_000 // N))):
        b = max(1, 2_000_000 // N)
        z = g.standard_normal((b, N)); u = np.sort(Phi((z - z.mean(1, keepdims=True)) / z.std(1, ddof=1, keepdims=True)), 1)
        i = np.arange(1, N + 1); out.append(np.maximum((i / N - u).max(1), (u - (i - 1) / N).max(1)))
    return np.concatenate(out)[:B]
def L_grid(zstd, zs):
    """L da amostra ampliada para cada z da grade (vetorizado). zstd: pontos padronizados da amostra original."""
    n = len(zstd); N = n + 1
    rho = np.sqrt((n - 1) / n + zs ** 2 / N)
    Y = np.concatenate([(zstd[None, :] - zs[:, None] / N) / rho[:, None], (zs * n / N / rho)[:, None]], 1)
    U = np.sort(Phi(Y), 1); i = np.arange(1, N + 1)
    return np.maximum((i / N - U).max(1), (U - (i - 1) / N).max(1))
res = {}
for N in [25, 50, 100, 200, 400, 1000]:
    n = N - 1; cN = np.quantile(Lnull(N, 40000), 0.95); cn = np.quantile(Lnull(n, 40000), 0.95)
    zs = np.arange(0, 6 * N ** 0.25 + 6, 0.01); vals = []
    while len(vals) < 300:
        x = g.standard_normal(n); zstd = (x - x.mean()) / x.std(ddof=1)
        u = np.sort(Phi(zstd)); i = np.arange(1, n + 1)
        if max((i / n - u).max(), (u - (i - 1) / n).max()) > cn: continue
        D = L_grid(zstd, zs) > cN
        if not D[-1]: vals.append(np.inf); continue
        vals.append(zs[np.where(~D)[0].max() + 1])
    v = np.array(vals) / N ** 0.25
    res[N] = dict(q=[float(np.quantile(v, p)) for p in (0.05, 0.25, 0.5, 0.75, 0.95)], sqrtN_c=float(np.sqrt(N) * cN))
    print(N, res[N], flush=True)
# distribuição-limite de C aproximada com o processo de Lilliefors em N grande
NB = 20000; k = np.sqrt(NB) * np.quantile(Lnull(NB, 3000), 0.95); print("k_alpha ~", k)
u_grid = np.linspace(-5, 5, 4001); cs = np.arange(0, 6, 0.002); Cs = []
drift = (cs[:, None] ** 2 / 2) * (u_grid * phi(u_grid))[None, :]
while len(Cs) < 1500:
    x = g.standard_normal(NB); zstd = np.sort((x - x.mean()) / x.std(ddof=1))
    Fn = np.searchsorted(zstd, u_grid, side="right") / NB
    Z = np.sqrt(NB) * (Fn - Phi(u_grid))
    if np.abs(Z).max() > k: continue
    S = np.abs(Z[None, :] + drift).max(1) > k
    Cs.append(cs[np.where(~S)[0].max() + 1] if S[-1] else np.inf)
Cs = np.array(Cs)
res["limite"] = dict(k_alpha=float(k), q=[float(np.quantile(Cs, p)) for p in (0.05, 0.25, 0.5, 0.75, 0.95)],
                     deterministico=float(np.sqrt(2 * k / phi(1))))
print("C:", res["limite"])
json.dump(res, open("zbar.json", "w"), indent=1)
