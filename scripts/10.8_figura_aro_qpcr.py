#!/usr/bin/env python3

import os
import pandas as pd
import matplotlib.pyplot as plt


INPUT = "definitivo/resultados/comparacion_12_aro_blast_qpcr.tsv"
OUT_DIR = "definitivo/figuras"

os.makedirs(OUT_DIR, exist_ok=True)


# Leer resultados definitivos
df = pd.read_csv(INPUT, sep="\t")


# Solo muestras colonizadas
df = df[df["grupo"] == "Colonizado"].copy()


# Convertir qPCR a número
qpcr = [c for c in df.columns if "qPCR" in c][0]

df = df.dropna(subset=["numero_12_aro", qpcr])


# Figura
fig, ax = plt.subplots(figsize=(8, 6))

ax.scatter(
    df["numero_12_aro"],
    df[qpcr],
    s=45,
    alpha=0.7
)

ax.set_yscale("log")

ax.set_xlabel("Número de determinantes ARO detectados")
ax.set_ylabel("Carga de K. pneumoniae (copias/ng ADN)")

ax.set_title(
    "Relación entre determinantes ARO y carga de K. pneumoniae"
)

ax.set_xticks(range(13))

ax.grid(
    True,
    linestyle="--",
    alpha=0.3
)

plt.tight_layout()


# Guardar
fig.savefig(
    f"{OUT_DIR}/10.8_Figura_aro_qpcr.png",
    dpi=300,
    bbox_inches="tight"
)

fig.savefig(
    f"{OUT_DIR}/10.8_Figura_aro_qpcr.pdf",
    bbox_inches="tight"
)

plt.close()

print("Figura 4 generada correctamente.")
print(f"Muestras colonizadas: {len(df)}")
