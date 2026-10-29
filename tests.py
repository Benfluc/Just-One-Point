"""Estatísticas de teste usadas no artigo: KS clássico, Lilliefors, R (Terán-García e Pérez-Fernández, 2024)
e T (Iturrate-Bobes, Pérez-Fernández e De Baets, 2026). Valores críticos por Monte Carlo."""
import numpy as np
from scipy import stats
rng = np.random.default_rng(1234)
Phi = stats.norm.cdf

def ks_sup(xs, G):
    """sup_x |F_n(x)-G(x)| para xs ordenado (aceita empates); G = valores da fda nos pontos."""
    n = len(xs); i = np.arange(1, n + 1)
    return max((i / n - G).max(), (G - (i - 1) / n).max())

def KS(x):  # F0 = N(0,1) totalmente especificada
    xs = np.sort(x); return ks_sup(xs, Phi(xs))

def LILL(x):
    xs = np.sort(x); return ks_sup(xs, Phi((xs - xs.mean()) / xs.std(ddof=1)))

def sup_restrito(xs, mu, sig, a, b):
    """sup_{x in [xs[a], xs[b]]} |F_n(x) - Phi((x-mu)/sig)|, com empates."""
    n = len(xs); lo, hi = xs[a], xs[b]
    u = np.unique(xs[(xs >= lo) & (xs <= hi)])
    G = Phi((u - mu) / sig)
    Fd = np.searchsorted(xs, u, side="right") / n
    Fe = np.searchsorted(xs, u, side="left") / n
    d = np.abs(Fd - G).max()
    m = u > lo
    if m.any(): d = max(d, np.abs(Fe[m] - G[m]).max())
    return d

def janela_sd(xs, m):
    c1 = np.concatenate([[0], np.cumsum(xs)]); c2 = np.concatenate([[0], np.cumsum(xs ** 2)])
    s1 = c1[m:] - c1[:-m]; s2 = c2[m:] - c2[:-m]
    return int(np.argmin((s2 - s1 ** 2 / m)))

class R_test:
    def __init__(self, N, r=0.1, B=50_000):
        self.N, self.r, self.m = N, r, int(np.floor((1 - r) * N))
        Z = np.sort(rng.standard_normal((B, N)), 1); m = self.m
        c1 = np.concatenate([np.zeros((B, 1)), np.cumsum(Z, 1)], 1); c2 = np.concatenate([np.zeros((B, 1)), np.cumsum(Z ** 2, 1)], 1)
        s1 = c1[:, m:] - c1[:, :-m]; s2 = c2[:, m:] - c2[:, :-m]
        a = np.argmin(s2 - s1 ** 2 / m, 1)
        idx = a[:, None] + np.arange(m)
        K = np.take_along_axis(Z, idx, 1)
        self.c_sigma = 1 / K.std(1, ddof=1).mean()
        mu = K.mean(1); sg = self.c_sigma * K.std(1, ddof=1)
        # sem empates sob H0: F_n(x_(j)) = j/N
        j = np.arange(1, N + 1)[None, :]
        G = Phi((Z - mu[:, None]) / sg[:, None])
        inside = (j - 1 >= a[:, None]) & (j - 1 <= (a + m - 1)[:, None])
        dd = np.where(inside, np.abs(j / N - G), 0)
        de = np.where(inside & (j - 1 > a[:, None]), np.abs((j - 1) / N - G), 0)
        self.null = np.maximum(dd.max(1), de.max(1))
    def stat(self, x):
        xs = np.sort(x); a = janela_sd(xs, self.m); K = xs[a:a + self.m]
        return sup_restrito(xs, K.mean(), self.c_sigma * K.std(ddof=1), a, a + self.m - 1)

class T_test:
    """OWA: v = estatística de ordem i0, w = mediana."""
    def __init__(self, N, frac=0.075, B=50_000):
        self.N = N; self.i0 = max(1, int(np.ceil(frac * N)))
        Z = np.sort(rng.standard_normal((B, N)), 1); a, b = self.i0 - 1, N - self.i0
        self.c_sigma = 1 / (Z[:, b] - Z[:, a]).mean()
        mu = np.median(Z, 1); sg = self.c_sigma * (Z[:, b] - Z[:, a])
        j = np.arange(1, N + 1)[None, :]
        G = Phi((Z - mu[:, None]) / sg[:, None])
        inside = (j - 1 >= a) & (j - 1 <= b)
        dd = np.where(inside, np.abs(j / N - G), 0)
        de = np.where(inside & (j - 1 > a), np.abs((j - 1) / N - G), 0)
        self.null = np.maximum(dd.max(1), de.max(1))
    def stat(self, x):
        xs = np.sort(x); a, b = self.i0 - 1, self.N - self.i0
        return sup_restrito(xs, np.median(xs), self.c_sigma * (xs[b] - xs[a]), a, b)

def nulos_simples(N, B=100_000):
    Z = rng.standard_normal((B, N)); Zs = np.sort(Z, 1); i = np.arange(1, N + 1)
    U = Phi(Zs); ks = np.maximum((i / N - U).max(1), (U - (i - 1) / N).max(1))
    U = np.sort(Phi((Z - Z.mean(1, keepdims=True)) / Z.std(1, ddof=1, keepdims=True)), 1)
    li = np.maximum((i / N - U).max(1), (U - (i - 1) / N).max(1))
    return ks, li

class Bateria:
    """Os cinco testes para um tamanho N, com valores críticos de nível alpha."""
    def __init__(self, N, alpha=0.05):
        self.N = N; ks, li = nulos_simples(N)
        self.R = R_test(N); self.T = T_test(N)
        q = lambda v: np.quantile(v, 1 - alpha)
        self.c = dict(KS=stats.kstwo.ppf(1 - alpha, N), Lilliefors=q(li), R=q(self.R.null), T=q(self.T.null))
        self.nulls = dict(Lilliefors=li, R=self.R.null, T=self.T.null)
        self.alpha = alpha
    def rejeita(self, x):
        return dict(KS=KS(x) > self.c["KS"], Lilliefors=LILL(x) > self.c["Lilliefors"],
                    SW=stats.shapiro(x).pvalue < self.alpha,
                    R=self.R.stat(x) > self.c["R"], T=self.T.stat(x) > self.c["T"])
