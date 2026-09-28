#!/usr/bin/env python3

import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


BASE = "/home/cgamboah/TFM_resistoma"

# Matriz definitiva del análisis principal
archivo = (
    f"{BASE}/definitivo/matrices/profundidad_1/"
    "matriz_presencia_ausencia.tsv"
)

# Resultado del análisis de sensibilidad
salida = (
    f"{BASE}/definitivo/resultados/"
    "analisis_sensibilidad_11_qpcr_excluidas.tsv"
)

# ---------------------------------------------------------
# 1. Cargar matriz definitiva
# ---------------------------------------------------------

datos = pd.read_csv(archivo, sep="\t")

print("=" * 70)
print("ANÁLISIS DE SENSIBILIDAD DE 9.7")
print("Exclusión de 11 muestras cultivo negativo/qPCR positivo")
print("=" * 70)

print(f"Muestras iniciales: {len(datos)}")

# ---------------------------------------------------------
# 2. Muestras excluidas
# ---------------------------------------------------------

excluir = {
    "ERR9883616",
    "ERR9880783",
    "ERR9880715",
    "ERR9880724",
    "ERR9880722",
    "ERR9880786",
    "ERR9880789",
    "ERR9880795",
    "ERR9882318",
    "ERR9882526",
    "ERR9882649"
}

# Verificación
encontradas = set(datos["run"]) & excluir
faltantes = excluir - set(datos["run"])

if len(encontradas) != 11 or faltantes:
    raise SystemExit(
        "ERROR: no se encontraron las 11 muestras esperadas "
        "en la matriz definitiva."
    )

# ---------------------------------------------------------
# 3. Excluir las 11 muestras
# ---------------------------------------------------------

sensibilidad = datos[
    ~datos["run"].isin(excluir)
].copy()

print(f"Muestras excluidas: {len(excluir)}")
print(f"Muestras para sensibilidad: {len(sensibilidad)}")

print("\nGrupos después de la exclusión:")
print(sensibilidad["grupo"].value_counts())

# Verificación del tamaño esperado
n_colonizado = (
    sensibilidad["grupo"] == "Colonizado"
).sum()

n_no_colonizado = (
    sensibilidad["grupo"] == "No colonizado"
).sum()

if n_colonizado != 52 or n_no_colonizado != 40:
    raise SystemExit(
        "ERROR: los tamaños de grupo no coinciden con "
        "52 colonizados y 40 no colonizados."
    )

# ---------------------------------------------------------
# 4. Identificar ARO
# ---------------------------------------------------------

columnas_muestra = [
    "muestra",
    "run",
    "grupo"
]

aros = [
    c for c in sensibilidad.columns
    if c not in columnas_muestra
]

print(f"\nARO analizados: {len(aros)}")

if len(aros) != 266:
    raise SystemExit(
        f"ERROR: se esperaban 266 ARO y se encontraron {len(aros)}."
    )

# ---------------------------------------------------------
# 5. Fisher para cada ARO
# ---------------------------------------------------------

resultados = []

colonizado = sensibilidad[
    sensibilidad["grupo"] == "Colonizado"
]

no_colonizado = sensibilidad[
    sensibilidad["grupo"] == "No colonizado"
]

for nombre in aros:

    positivos_col = int(
        colonizado[nombre].sum()
    )

    negativos_col = (
        len(colonizado) - positivos_col
    )

    positivos_no = int(
        no_colonizado[nombre].sum()
    )

    negativos_no = (
        len(no_colonizado) - positivos_no
    )

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
        "odds_ratio_fisher": odds_ratio,
        "p": p
    })

# ---------------------------------------------------------
# 6. Corrección Benjamini-Hochberg
# ---------------------------------------------------------

resultados = pd.DataFrame(resultados)

resultados["FDR"] = multipletests(
    resultados["p"],
    method="fdr_bh"
)[1]

# Ordenar por FDR
resultados = resultados.sort_values(
    "FDR"
).reset_index(drop=True)

# ---------------------------------------------------------
# 7. Mostrar resultados significativos
# ---------------------------------------------------------

significativos = resultados[
    resultados["FDR"] < 0.05
].copy()

print("\n" + "=" * 70)
print("RESULTADOS DE SENSIBILIDAD")
print("=" * 70)

print(
    f"ARO con FDR < 0,05: {len(significativos)}"
)

if len(significativos) > 0:

    print(
        significativos[
            [
                "ARO",
                "positivos_colonizado",
                "total_colonizado",
                "positivos_no_colonizado",
                "total_no_colonizado",
                "p",
                "FDR"
            ]
        ].to_string(index=False)
    )

# ---------------------------------------------------------
# 8. Guardar resultados
# ---------------------------------------------------------

resultados.to_csv(
    salida,
    sep="\t",
    index=False
)

print("\nArchivo guardado:")
print(salida)
