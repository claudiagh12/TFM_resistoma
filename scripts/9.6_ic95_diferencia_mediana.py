#!/usr/bin/env python3 
 
import pandas as pd 
import numpy as np 
 
BASE = "/home/cgamboah/TFM_resistoma" 
 
N_BOOT = 10000 
SEED = 20260927 
 
rng = np.random.default_rng(SEED) 
 
for umbral in [1, 5, 10]: 
 
    archivo = ( 
        f"{BASE}/definitivo/matrices/profundidad_{umbral}" 
        f"/riqueza_por_muestra.tsv" 
    ) 
 
    datos = pd.read_csv(archivo, sep="\t") 
 
    colonizado = datos.loc[ 
        datos["grupo"] == "Colonizado", "riqueza" 
    ].to_numpy() 
 
    no_colonizado = datos.loc[ 
        datos["grupo"] == "No colonizado", "riqueza" 
    ].to_numpy() 
 
    diferencia_observada = ( 
        np.median(colonizado) - 
        np.median(no_colonizado) 
    ) 
 
    diferencias = np.empty(N_BOOT) 
 
    for i in range(N_BOOT): 
        muestra_col = rng.choice( 
            colonizado, 
            size=len(colonizado), 
            replace=True 
        ) 
 
        muestra_no = rng.choice( 
            no_colonizado, 
            size=len(no_colonizado), 
            replace=True 
        ) 
 
        diferencias[i] = ( 
            np.median(muestra_col) - 
            np.median(muestra_no) 
        ) 
 
    ic_inf, ic_sup = np.percentile( 
        diferencias, 
        [2.5, 97.5] 
    ) 
 
    print("\n" + "=" * 60) 
    print(f"DEPTH >= {umbral}") 
    print("=" * 60) 
    print(f"n colonizado = {len(colonizado)}") 
    print(f"n no colonizado = {len(no_colonizado)}") 
    print(f"Diferencia observada = {diferencia_observada:.2f}") 
    print(f"IC95% bootstrap = {ic_inf:.2f} - {ic_sup:.2f}") 
 
print("\nAnálisis completado.") 
print(f"Bootstrap: {N_BOOT} remuestreos") 
print(f"Semilla: {SEED}") 
