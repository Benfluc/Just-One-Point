"""Verificações determinísticas da prova revista do Teorema 1 e da Proposição sobre Shapiro-Wilk."""
import numpy as np, json, warnings
from scipy import stats, integrate
warnings.filterwarnings("ignore")
Phi = stats.norm.cdf; phi0, phi1 = stats.norm.pdf(0), stats.norm.pdf(1)
def g(N, k): return (N - k) / N - Phi(-np.sqrt((N - 1) * k / (N * (N - k))))
def Dinf(N): return g(N, 1)
out = {}
# (1) Lema do supremo: max_k g = g(1) para N <= 34
marg = []
for N in range(4, 41):
    ks = np.arange(1, N); v = g(N, ks)
    marg.append(dict(N=N, kstar=int(ks[np.argmax(v)]), sup=float(v.max()), Dinf=float(Dinf(N)),
                     folga=float(Dinf(N) - v[1:].max()) if N > 2 else None))
out["supremo"] = marg
print("N com k*=1:", [m["N"] for m in marg if m["kstar"] == 1][-1], "| menor folga em 4..34:",
      min(m["folga"] for m in marg if m["N"] <= 34))
# outras cotas do lema: casos v>=0 e i=N
print("checa 1/2-1/N < Dinf e Phi(-1/sqrtN) < Dinf para 4..34:",
      all(0.5 - 1 / N < Dinf(N) and Phi(-1 / np.sqrt(N)) < Dinf(N) for N in range(4, 35)))
# (2) cota analítica com limiar d: P(L>=d) <= 2exp(-2Nt^2) + pi_N(t;d)
def pi(N, t, d):
    k = N - 1
    def f(r):
        lim = r * (d - t - phi1 * abs(np.log(r))) / phi0
        pm = 1.0 if lim <= 0 else 2 * stats.norm.sf(np.sqrt(N) * lim)
        return pm * stats.chi2.pdf(k * r * r, k) * 2 * k * r
    v, err = integrate.quad(f, 1e-8, 20, limit=1000, points=[np.exp(-(d - t) / phi1), 1.0, np.exp((d - t) / phi1)], epsabs=1e-12)
    return v + stats.chi2.cdf(k * 1e-16, k) + stats.chi2.sf(k * 400, k), err
def cota(N, d):
    best = (9, None)
    for t in np.linspace(0.15, 0.45, 61):
        if np.exp(-2 * N * t * t) > 0.5: continue          # região de validade de Massart
        p, e = pi(N, t, d); b = 2 * np.exp(-2 * N * t * t) + p
        if b < best[0]: best = (b, float(t))
    return best
def crude(N, t, d):
    dd = d - t; c = dd / (2 * phi1); kap = np.exp(-c) * dd / (2 * phi0); k = N - 1
    eu = np.exp(2 * c) - 1; el = 1 - np.exp(-2 * c)
    return (2 * np.exp(-2 * N * t * t) + 2 * stats.norm.sf(kap * np.sqrt(N))
            + np.exp(-k * (eu - np.log1p(eu)) / 2) + np.exp(-k * (-el - np.log1p(-el)) / 2))
# N1: cota grosseira (decrescente em N) com limiar 1/2 abaixo de 0,01
best = None
for t in np.linspace(0.12, 0.4, 57):
    N = 10
    while crude(N, t, 0.5) >= 0.01: N += 1
    if best is None or N < best[0]: best = (N, float(t))
N1, t1 = best; print("N1 (alpha=0,01) =", N1, "t =", t1)
tab = []
for N in range(35, N1 + 1):
    d = float(Dinf(N)); b, t = cota(N, d); tab.append(dict(N=N, Dinf=d, cota=b, t=t))
print("max cota em 35..N1:", max(r["cota"] for r in tab), "em N=", max(tab, key=lambda r: r["cota"])["N"])
print("primeiras:", [(r["N"], round(r["cota"], 4)) for r in tab[:6]])
out["cotas"] = tab; out["N1"] = N1; out["t1"] = t1
for a in [0.10, 0.05]:
    Nmin = next(N for N in range(9, 200) if all(cota(M, float(Dinf(M)))[0] < a for M in range(N, min(N + 15, N1 + 1))))
    out[f"N0_{a}"] = Nmin; print("alpha", a, "cota < alpha a partir de N =", Nmin)
# (3) Shapiro-Wilk: W_inf = W da configuração limite (N-1 pontos empatados e 1 afastado)
sw = []
for N in range(4, 5001):
    r = stats.shapiro(np.r_[np.zeros(N - 1), 1.0]); sw.append((N, float(r.statistic), float(r.pvalue)))
pmax = max(sw, key=lambda s: s[2]); print("SW: maior p-valor no limite", pmax)
out["sw"] = dict(pmax=pmax, amostra={str(N): (w, p) for N, w, p in sw if N in (4, 5, 6, 10, 25, 50, 100, 1000, 5000)})
# confere W_inf contra amostras com ponto muito afastado
rng = np.random.default_rng(1)
for N in [4, 10, 100]:
    b = rng.standard_normal(N - 1); print(N, stats.shapiro(np.append(b, 1e7)).statistic, dict(sw=[s for s in sw if s[0] == N]))
json.dump(out, open("prova_v2.json", "w"), indent=1)
