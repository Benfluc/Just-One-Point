#!/bin/bash
set -e
python3 prova_v2.py
python3 calib.py
python3 exp_a.py
for N in 25 50 100; do python3 exp2v2.py $N; python3 contam.py $N; done
python3 mc_peq.py
python3 zbar.py
python3 figs.py
python3 tabelas.py
