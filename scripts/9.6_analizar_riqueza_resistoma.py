#!/usr/bin/env python3

import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu


BASE = "/home/cgamboah/TFM_resistoma"

umbrales = [1, 5, 10]

for umbral in umbrales:

    archivo = (
        f"{BASE}/definitivo/matrices/profundidad_{umbral}"
        f"/matriz_presencia_ausencia.tsv"
    )

    datos = pd.read_csv(archivo, sep="\t")

    # Las primeras columnas contienen información de la muestra
    columnas_grupo = ["muestra", "run", "grupo"]

    aro = [
        columna for columna in datos.columns
        if columna not in columnas_grupo
    ]

    # Riqueza = número de ARO presentes en cada muestra
    datos["riqueza"] = datos[aro].sum(axis=1)

    colonizado = datos.loc[
        datos["grupo"] == "Colonizado", "riqueza"
    ]

    no_colonizado = datos.loc[
        datos["grupo"] == "No colonizado", "riqueza"
    ]

    resultado = mannwhitneyu(
        colonizado,
        no_colonizado,
        alternative="two-sided"
    )

    # Tamaño de efecto r = |Z| / sqrt(N)
    n_total = len(colonizado) + len(no_colonizado)

    # Conversión de U a Z
    media_u = len(colonizado) * len(no_colonizado) / 2
    desviacion_u = np.sqrt(
        len(colonizado)
        * len(no_colonizado)
        * (n_total + 1)
        / 12
    )

    z = (resultado.statistic - media_u) / desviacion_u
    r = abs(z) / np.sqrt(n_total)

    diferencia_mediana = (
        colonizado.median() - no_colonizado.median()
    )

    print("\n" + "=" * 60)
    print(f"DEPTH >= {umbral}")
    print("=" * 60)

    print("\nColonizado:")
    print(f"  n = {len(colonizado)}")
    print(f"  mínimo = {colonizado.min()}")
    print(f"  máximo = {colonizado.max()}")
    print(f"  media = {colonizado.mean():.2f}")
    print(f"  mediana = {colonizado.median():.2f}")

    print("\nNo colonizado:")
    print(f"  n = {len(no_colonizado)}")
    print(f"  mínimo = {no_colonizado.min()}")
    print(f"  máximo = {no_colonizado.max()}")
    print(f"  media = {no_colonizado.mean():.2f}")
    print(f"  mediana = {no_colonizado.median():.2f}")

    print("\nComparación Mann–Whitney:")
    print(f"  U = {resultado.statistic:.2f}")
    print(f"  p = {resultado.pvalue:.6g}")
    print(f"  Z = {z:.4f}")
    print(f"  r = {r:.4f}")
    print(f"  diferencia de medianas = {diferencia_mediana:.2f}")

    # Guardar riqueza por muestra
    salida = (
        f"{BASE}/definitivo/matrices/profundidad_{umbral}"
        f"/riqueza_por_muestra.tsv"
    )

    datos[
        ["muestra", "run", "grupo", "riqueza"]
    ].to_csv(
        salida,
        sep="\t",
        index=False
    )

print("\nAnálisis completado.")
