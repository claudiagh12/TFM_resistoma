import csv
import numpy as np
from scipy.stats import spearmanr
from pathlib import Path

BASE = Path.home() / "TFM_resistoma"

# Profundidad de secuenciación
profundidad = {}

with open(
    BASE / "definitivo" / "matrices" / "profundidad_secuenciacion.tsv",
    encoding="utf-8"
) as archivo:

    lector = csv.DictReader(archivo, delimiter="\t")

    for fila in lector:
        profundidad[fila["run"]] = int(fila["pares_despues"])


for nivel in [1, 5, 10]:

    archivo = (
        BASE / "definitivo" / "matrices"
        / f"profundidad_{nivel}"
        / "matriz_presencia_ausencia.tsv"
    )

    profundidades = []
    riquezas = []

    with open(archivo, encoding="utf-8") as entrada:

        lector = csv.DictReader(entrada, delimiter="\t")

        for fila in lector:

            run = fila["run"]

            riqueza = sum(
                int(valor)
                for nombre, valor in fila.items()
                if nombre not in {"muestra", "run", "grupo"}
            )

            profundidades.append(profundidad[run])
            riquezas.append(riqueza)

    rho, p = spearmanr(profundidades, riquezas)

    print()
    print(f"Depth >= {nivel}")
    print(f"  n = {len(riquezas)}")
    print(f"  rho de Spearman = {rho:.4f}")
    print(f"  p = {p:.6g}")
