"""Gera as tabelas do artigo (PT e EN) a partir dos arquivos JSON. Uso: python tabelas.py"""
import json, math
T = ["Lilliefors", "R", "T", "SW", "KS"]
NOME = {"pt": {"Lilliefors": "Lilliefors", "R": "$R$", "T": "$T$", "SW": "Shapiro-Wilk", "KS": "KS clássico"},
        "en": {"Lilliefors": "Lilliefors", "R": "$R$", "T": "$T$", "SW": "Shapiro-Wilk", "KS": "Classical KS"}}
def wilson(k, n, z=1.96):
    if n == 0: return (float("nan"),) * 3
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return 100 * p, 100 * (c - h), 100 * (c + h)
def f(x, lang, nd=1):
    s = f"{x:.{nd}f}"; return s.replace(".", ",") if lang == "pt" else s
def ci(k, n, lang):
    if n < 30: return "—"
    p, lo, hi = wilson(k, n); sep = "; " if lang == "pt" else "; "
    return f"{f(p, lang)} [{f(lo, lang)}{sep}{f(hi, lang)}]"
def tab_exp2(lang, cal):
    H = {"pt": "| $N$ | $h$ | Teste | Rej. antes | Rej. após cópia | $n_{\\mathrm{nr}}$ | Inv. cópia da moda | Inv. moda $+\\,h/100$ | Inv. obs. nova |",
         "en": "| $N$ | $h$ | Test | Rej. before | Rej. after copy | $n_{\\mathrm{nr}}$ | Flips, copy of mode | Flips, mode $+\\,h/100$ | Flips, new obs. |"}[lang]
    L = [H, "|:----:|:----:|:------------|:-----:|:------:|:-----:|:---------:|:---------:|:---------:|"]
    for N in [25, 50, 100]:
        d = json.load(open(f"exp2v2_{N}.json"))["res"]
        for h in ["0.25", "0.5", "1.0"]:
            r = d[h][cal]; first = True
            if sum(r["cond"][k]["nr"] < 30 for k in T) >= 2:
                continue
            for k in T:
                c = r["cond"][k]; nr = c["nr"]
                hh = f(float(h), lang, 2 if h == "0.25" else 1)
                L.append(f"| {N if first else ''} | {hh if first else ''} | {NOME[lang][k]} | {f(100*r['base'][k], lang)} | {f(100*r['depois']['moda'][k], lang)} | {nr} | {ci(c['moda'], nr, lang)} | {ci(c['perto'], nr, lang)} | {ci(c['nova'], nr, lang)} |")
                first = False
    return "\n".join(L)
def tab_comum(lang, cal):
    H = {"pt": "| $N$ | $h$ | $n_{\\mathrm{comum}}$ | Lilliefors | $R$ | $T$ | Shapiro-Wilk | KS clássico |",
         "en": "| $N$ | $h$ | $n_{\\mathrm{common}}$ | Lilliefors | $R$ | $T$ | Shapiro-Wilk | Classical KS |"}[lang]
    L = [H, "|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|"]
    for N in [25, 50, 100]:
        d = json.load(open(f"exp2v2_{N}.json"))["res"]
        for h in ["0.25", "0.5", "1.0"]:
            c = d[h][cal]["comum"]; nr = c["nr"]
            if nr < 30: continue
            hh = f(float(h), lang, 2 if h == "0.25" else 1)
            L.append(f"| {N} | {hh} | {nr} | " + " | ".join(f"{f(100*c[k]['moda']/nr, lang)} ({f(100*c[k]['nova']/nr, lang)})" for k in T) + " |")
    return "\n".join(L)
def tab_contam(lang):
    F = ["N(0,1)", "N(0,3^2)", "N(0,5^2)", "t_3"]
    H = {"pt": "| $N$ | Teste | $n_{\\mathrm{nr}}$ | $N(0,1)$ | $N(0,3^2)$ | $N(0,5^2)$ | $t_3$ |",
         "en": "| $N$ | Test | $n_{\\mathrm{nr}}$ | $N(0,1)$ | $N(0,3^2)$ | $N(0,5^2)$ | $t_3$ |"}[lang]
    L = [H, "|:----:|:------------|:-----:|:---------:|:---------:|:---------:|:---------:|"]
    for N in [25, 50, 100]:
        d = json.load(open(f"contam_{N}.json"))["res"]["cond"]; first = True
        for k in T:
            c = d[k]; L.append(f"| {N if first else ''} | {NOME[lang][k]} | {c['nr']} | " + " | ".join(ci(c[x], c["nr"], lang) for x in F) + " |"); first = False
    return "\n".join(L)
def tab_calib(lang):
    d = json.load(open("calib.json"))
    H = {"pt": "| $N$ | $c$ KS (exato) | $c$ Lilliefors [IC 95%] | $c$ $R$ [IC 95%] | $c$ $T$ [IC 95%] | $i_0$ | Tamanho (%): KS / L / SW / $R$ / $T$ |",
         "en": "| $N$ | $c$ KS (exact) | $c$ Lilliefors [95% CI] | $c$ $R$ [95% CI] | $c$ $T$ [95% CI] | $i_0$ | Size (%): KS / L / SW / $R$ / $T$ |"}[lang]
    L = [H, "|:----:|:------:|:------------:|:------------:|:------------:|:---:|:------------:|"]
    for N, r in d.items():
        g = lambda k: f"{f(r[k]['c'], lang, 4)} [{f(r[k]['lo'], lang, 4)}; {f(r[k]['hi'], lang, 4)}]"
        tm = " / ".join(f(100 * r["tamanho"][k], lang) for k in ["KS", "Lilliefors", "SW", "R", "T"])
        L.append(f"| {N} | {f(r['KS']['c'], lang, 4)} | {g('Lilliefors')} | {g('R')} | {g('T')} | {r['c_sigma']['i0']} | {tm} |")
    return "\n".join(L)
def tab_A1(lang):
    d = json.load(open("prova_v2.json"))
    sel = [r for r in d["cotas"] if r["N"] in (35, 36, 38, 40, 50, 75, 100, 125, 153)]
    H = {"pt": "| $N$ | $D_\\infty(N-1)$ | $t$ ótimo | Cota para $P_{H_0}(L_N\\ge D_\\infty)$ |",
         "en": "| $N$ | $D_\\infty(N-1)$ | Optimal $t$ | Bound on $P_{H_0}(L_N\\ge D_\\infty)$ |"}[lang]
    L = [H, "|:---:|:---:|:---:|:---:|"]
    for r in sel:
        b = r["cota"]; m, e = f"{b:.1e}".split("e"); bs = f(b, lang, 4) if b >= 1e-4 else "$" + (m.replace(".", "{,}") if lang == "pt" else m) + "\\times10^{" + str(int(e)) + "}$"
        L.append(f"| {r['N']} | {f(r['Dinf'], lang, 4)} | {f(r['t'], lang, 3)} | {bs} |")
    return "\n".join(L)
def tab_zbar(lang):
    d = json.load(open("zbar.json"))
    H = {"pt": "| $N$ | P5 | P25 | Mediana | P75 | P95 |", "en": "| $N$ | P5 | P25 | Median | P75 | P95 |"}[lang]
    L = [H, "|:---:|:---:|:---:|:---:|:---:|:---:|"]
    for N in ["25", "50", "100", "200", "400", "1000"]:
        L.append(f"| {N} | " + " | ".join(f(q, lang, 2) for q in d[N]["q"]) + " |")
    lab = {"pt": "Limite $C$", "en": "Limit $C$"}[lang]
    L.append(f"| {lab} | " + " | ".join(f(q, lang, 2) for q in d["limite"]["q"]) + " |")
    return "\n".join(L)
if __name__ == "__main__":
    out = {}
    for lang in ["pt", "en"]:
        out[f"exp2_nom_{lang}"] = tab_exp2(lang, "nom"); out[f"exp2_cal_{lang}"] = tab_exp2(lang, "cal")
        out[f"comum_cal_{lang}"] = tab_comum(lang, "cal"); out[f"comum_nom_{lang}"] = tab_comum(lang, "nom")
        out[f"contam_{lang}"] = tab_contam(lang); out[f"calib_{lang}"] = tab_calib(lang); out[f"A1_{lang}"] = tab_A1(lang)
        try: out[f"zbar_{lang}"] = tab_zbar(lang)
        except FileNotFoundError: pass
    json.dump(out, open("tabelas.json", "w"), ensure_ascii=False, indent=1)
    for k in ["exp2_nom_pt", "exp2_cal_pt", "comum_cal_pt", "comum_nom_pt", "calib_pt", "A1_pt"]: print("###", k); print(out[k])

def tab_exp2_compacta(lang):
    H = {"pt": "| $N$ | $h$ | Teste | Nom.: rej. antes | Nom.: inv. moda | Nom.: inv. nova | Cal.: rej. antes | Cal.: inv. moda | Cal.: inv. nova |",
         "en": "| $N$ | $h$ | Test | Nom.: rej. before | Nom.: flips, mode | Nom.: flips, new | Cal.: rej. before | Cal.: flips, mode | Cal.: flips, new |"}[lang]
    L = [H, "|:----:|:----:|:------------|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|"]
    for N in [25, 50, 100]:
        d = json.load(open(f"exp2v2_{N}.json"))["res"]
        for h in ["0.25", "0.5", "1.0"]:
            first = True
            for k in T:
                cells = []
                for cal in ["nom", "cal"]:
                    r = d[h][cal]; c = r["cond"][k]; nr = c["nr"]
                    cells += [f(100 * r["base"][k], lang)] + ([f(100 * c[a] / nr, lang) for a in ["moda", "nova"]] if nr >= 30 else ["—", "—"])
                hh = f(float(h), lang, 2 if h == "0.25" else 1)
                L.append(f"| {N if first else ''} | {hh if first else ''} | {NOME[lang][k]} | " + " | ".join(cells) + " |"); first = False
    return "\n".join(L)
def tab_extremo(lang):
    d = json.load(open("exp_a.json"))
    H = {"pt": "| $N=n+1$ | KS: existe $x$ que inverte | KS: satisfaz (6) | KS, $z=10^6$ | Lilliefors, $z=10^6$ | Shapiro-Wilk, $z=10^6$ | $R$, $z=10^6$ | $T$, $z=10^6$ |",
         "en": "| $N=n+1$ | KS: some $x$ flips | KS: satisfies (6) | KS, $z=10^6$ | Lilliefors, $z=10^6$ | Shapiro-Wilk, $z=10^6$ | $R$, $z=10^6$ | $T$, $z=10^6$ |"}[lang]
    L = [H, "|:-------:|:---------:|:---------:|:---------:|:---------:|:---------:|:---------:|:---------:|"]
    for N in ["10", "25", "50", "100"]:
        v = d[N]; kv, kb, nks = v["ks_counts"]; fc = v["far_counts"]
        cells = [ci(kv, nks, lang), ci(kb, nks, lang)] + [ci(fc[k][1], fc[k][0], lang) for k in ["KS", "Lilliefors", "SW", "R", "T"]]
        L.append(f"| {N} | " + " | ".join(cells) + " |")
    return "\n".join(L)
if __name__ == "__main__":
    for lang in ["pt", "en"]:
        out[f"exp2c_{lang}"] = tab_exp2_compacta(lang)
        try: out[f"extremo_{lang}"] = tab_extremo(lang)
        except (KeyError, FileNotFoundError): pass
    json.dump(out, open("tabelas.json", "w"), ensure_ascii=False, indent=1)
    print(out["exp2c_pt"])
