# Basta um ponto — código de reprodução

Código que reproduz todos os resultados numéricos, tabelas e figuras do artigo
"Basta um ponto: resistência do teste de Lilliefors e de suas versões robustas à inclusão de uma única observação".

## Ambiente
Python 3.13, NumPy 2.5, SciPy 1.18, Matplotlib 3.11 (ver `requirements.txt`).
Todas as simulações usam geradores `numpy.random.default_rng` com sementes fixas declaradas em cada script.

## Arquivos
| Script | Conteúdo | Saída |
|---|---|---|
| `testes.py` | Estatísticas KS, Lilliefors, R (Terán-García e Pérez-Fernández, 2024) e T (Iturrate-Bobes et al., 2026); calibração por Monte Carlo | — |
| `prova_v2.py` | Verificações determinísticas da prova do Teorema 1 (Lema 2, cotas para N ≥ 35) e da Proposição 2 (Shapiro-Wilk) | `prova_v2.json` |
| `calib.py` | Valores críticos com IC 95% e tamanho empírico dos cinco testes | `calib.json` |
| `exp_a.py` | Experimento 1: observação na cauda | `exp_a.json` |
| `exp2v2.py N` | Experimento 2: aglomerados por arredondamento (N = 25, 50, 100) | `exp2v2_N.json` |
| `contam.py N` | Experimento 3: contaminação probabilística (N = 25, 50, 100) | `contam_N.json` |
| `zbar.py` | Experimento 4: comportamento assintótico do limiar | `zbar.json` |
| `mc_peq.py` | Valores críticos de Lilliefors para 4 ≤ N ≤ 39 com 10^6 réplicas (Tabela 1 e Figura 1) | `mc_peq.json` |
| `figs.py` | Exemplos ilustrativos e Figuras 1–3 (PT e EN) | `exemplos.json`, `fig*.png` |
| `tabelas.py` | Gera as tabelas do artigo a partir dos arquivos JSON | `tabelas.md` |

Os resultados usados no artigo estão em `resultados/` e as figuras em `figuras/`. Os scripts gravam os arquivos JSON no diretório corrente.

## Execução
```
pip install -r requirements.txt
bash run_all.sh
```
O tempo total é de algumas horas em um computador com dois núcleos; os experimentos 2 e 3 podem ser executados em paralelo.
