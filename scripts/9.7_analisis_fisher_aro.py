#!/usr/bin/env python3

import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


BASE = "/home/cgamboah/TFM_resistoma"

archivo = (
    f"{BASE}/definitivo/matrices/profundidad_1/"
    "matriz_presencia_ausencia.tsv"
)

datos = pd.read_csv(archivo, sep="\t")

# Columnas que no corresponden a ARO
columnas_muestra = [
    "muestra",
    "run",
    "grupo"
]

aro = [
    c for c in datos.columns
    if c not in columnas_muestra
]

print("=" * 70)
print("ANÁLISIS DE PRESENCIA/AUSENCIA DE LOS ARO")
print("=" * 70)

print(f"Muestras: {len(datos)}")
print(f"ARO analizados: {len(aro)}")

# Separar grupos
colonizado = datos[datos["grupo"] == "Colonizado"]
no_colonizado = datos[datos["grupo"] == "No colonizado"]

resultados = []

for nombre in aro:

    positivos_col = int(colonizado[nombre].sum())
    negativos_col = len(colonizado) - positivos_col

    positivos_no = int(no_colonizado[nombre].sum())
    negativos_no = len(no_colonizado) - positivos_no

    tabla = [
        [positivos_col, negativos_col],
        [positivos_no, negativos_no]
    ]

    odds_ratio, p = fisher_exact(
        tabla,
        alternative="two-sided"
    )

    resultados.append({
        "ARO": nombre,
        "positivos_colonizado": positivos_col,
        "total_colonizado": len(colonizado),
        "prevalencia_colonizado": (
            positivos_col / len(colonizado) * 100
        ),
        "positivos_no_colonizado": positivos_no,
        "total_no_colonizado": len(no_colonizado),
        "prevalencia_no_colonizado": (
            positivos_no / len(no_colonizado) * 100
        ),
        "p": p
    })


resultados = pd.DataFrame(resultados)

# Corrección Benjamini-Hochberg
resultados["FDR"] = multipletests(
    resultados["p"],
    method="fdr_bh"
)[1]

# Ordenar por FDR
resultados = resultados.sort_values(
    "FDR"
).reset_index(drop=True)

# Significativos
significativos = resultados[
    resultados["FDR"] < 0.05
].copy()

print("\n" + "=" * 70)
print("RESULTADOS")
print("=" * 70)

print(
    f"ARO con FDR < 0,05: {len(significativos)}"
)

print("\nDeterminantes significativos:")

if len(significativos) > 0:

    print(
        significativos[
            [
                "ARO",
                "positivos_colonizado",
                "prevalencia_colonizado",
                "positivos_no_colonizado",
                "prevalencia_no_colonizado",
                "p",
                "FDR"
            ]
        ].to_string(index=False)
    )


# Guardar todos los ARO
salida_todos = (
    f"{BASE}/definitivo/resultados/"
    "analisis_fisher_266_aro.tsv"
)

resultados.to_csv(
    salida_todos,
    sep="\t",
    index=False
)

# Guardar solo significativos
salida_sig = (
    f"{BASE}/definitivo/resultados/"
    "aro_significativos_fdr05.tsv"
)

significativos.to_csv(
    salida_sig,
    sep="\t",
    index=False
)

print("\nArchivos guardados:")
print(salida_todos)
print(salida_sig)
