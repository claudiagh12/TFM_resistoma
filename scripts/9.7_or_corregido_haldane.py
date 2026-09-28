#!/usr/bin/env python3

import pandas as pd

BASE = "/home/cgamboah/TFM_resistoma"

archivo = (
    f"{BASE}/definitivo/resultados/"
    "aro_significativos_fdr05.tsv"
)

datos = pd.read_csv(archivo, sep="\t")

resultados = []

for _, fila in datos.iterrows():

    a = fila["positivos_colonizado"]
    b = fila["total_colonizado"] - a
    c = fila["positivos_no_colonizado"]
    d = fila["total_no_colonizado"] - c

    # Corrección de Haldane-Anscombe
    or_corregido = (
        (a + 0.5) * (d + 0.5)
        / ((b + 0.5) * (c + 0.5))
    )

    resultados.append({
        "ARO": fila["ARO"],
        "positivos_colonizado": a,
        "total_colonizado": fila["total_colonizado"],
        "positivos_no_colonizado": c,
        "total_no_colonizado": fila["total_no_colonizado"],
        "Fisher_p": fila["p"],
        "FDR": fila["FDR"],
        "OR_Haldane_Anscombe": or_corregido
    })

resultado = pd.DataFrame(resultados)

print("=" * 80)
print("OR CORREGIDO MEDIANTE HALDANE-ANScombe")
print("=" * 80)

print(
    resultado.to_string(index=False)
)

salida = (
    f"{BASE}/definitivo/resultados/"
    "or_corregido_haldane_12_aro.tsv"
)

resultado.to_csv(
    salida,
    sep="\t",
    index=False
)

print("\nArchivo guardado:")
print(salida)
