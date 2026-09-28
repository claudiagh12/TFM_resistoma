#!/usr/bin/env python3

import pandas as pd
from itertools import combinations

BASE = "/home/cgamboah/TFM_resistoma"

# Entradas definitivas
ARO_FILE = f"{BASE}/definitivo/resultados/aro_significativos_fdr05.tsv"
MATRIX_FILE = f"{BASE}/definitivo/matrices/profundidad_1/matriz_presencia_ausencia.tsv"

# Salidas
OUT = f"{BASE}/definitivo/resultados"

# 1. Leer los 12 ARO significativos
aros = pd.read_csv(ARO_FILE, sep="\t")["ARO"].tolist()

# 2. Leer matriz de presencia/ausencia
datos = pd.read_csv(MATRIX_FILE, sep="\t")

# 3. Seleccionar únicamente los 12 ARO
matriz = datos[aros]

# 4. Número de ARO por muestra
resumen = pd.DataFrame({
    "run": datos["run"],
    "n_significant_ARO": matriz.sum(axis=1)
})

resumen.to_csv(
    f"{OUT}/coocurrencia_12_aro_resumen.tsv",
    sep="\t",
    index=False
)

# 5. Coocurrencia y Jaccard
resultados = []

for aro1, aro2 in combinations(aros, 2):

    a = matriz[aro1].sum()
    b = matriz[aro2].sum()
    juntos = ((matriz[aro1] == 1) & (matriz[aro2] == 1)).sum()

    union = a + b - juntos

    jaccard = juntos / union if union > 0 else 0

    resultados.append([
        aro1,
        aro2,
        juntos,
        a,
        b,
        jaccard
    ])

pares = pd.DataFrame(
    resultados,
    columns=[
        "ARO_1",
        "ARO_2",
        "cooccurrence_n",
        "ARO_1_n",
        "ARO_2_n",
        "jaccard_index"
    ]
)

pares.to_csv(
    f"{OUT}/coocurrencia_12_aro_pares.tsv",
    sep="\t",
    index=False
)

# 6. Matriz de coocurrencia
cooc = matriz.T.dot(matriz)

cooc.to_csv(
    f"{OUT}/coocurrencia_12_aro_matriz.tsv",
    sep="\t"
)

print("Análisis de coocurrencia completado.")
print(f"ARO analizados: {len(aros)}")
print(f"Muestras: {len(datos)}")
print(f"Parejas: {len(pares)}")
