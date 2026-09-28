#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# Archivo de entrada
INPUT = "definitivo/matrices/profundidad_1/matriz_presencia_ausencia.tsv"

# Carpeta de salida
OUT_DIR = "definitivo/figuras"
os.makedirs(OUT_DIR, exist_ok=True)


# Leer matriz
df = pd.read_csv(INPUT, sep="\t")

# Identificar columnas ARO
metadata = ["muestra", "run", "grupo"]
aro_columns = [c for c in df.columns if c not in metadata]

# Calcular riqueza
df["n_ARO"] = df[aro_columns].sum(axis=1)

# Separar grupos
colonizado = df.loc[df["grupo"] == "Colonizado", "n_ARO"]
no_colonizado = df.loc[df["grupo"] == "No colonizado", "n_ARO"]


# Figura
fig, ax = plt.subplots(figsize=(7, 6))

data = [colonizado, no_colonizado]
labels = [
    f"Colonizado\n(n = {len(colonizado)})",
    f"No colonizado\n(n = {len(no_colonizado)})"
]

bp = ax.boxplot(
    data,
    labels=labels,
    widths=0.55,
    patch_artist=True,
    showfliers=False
)

# Colores
bp["boxes"][0].set_facecolor("#D97706")
bp["boxes"][1].set_facecolor("#A5A5A5")

for box in bp["boxes"]:
    box.set_alpha(0.75)

# Puntos individuales
rng = np.random.default_rng(12345)

for i, values in enumerate(data, start=1):
    jitter = rng.uniform(-0.12, 0.12, len(values))
    ax.scatter(
        i + jitter,
        values,
        s=25,
        color="black",
        alpha=0.65,
        zorder=3
    )

# Etiquetas
ax.set_xlabel("Estado de colonización")
ax.set_ylabel("Número de determinantes ARO por muestra")
ax.set_title("Riqueza del resistoma según el estado de colonización")

ax.yaxis.grid(True, linestyle="--", alpha=0.3)
ax.set_axisbelow(True)

plt.tight_layout()

# Guardar
fig.savefig(
    f"{OUT_DIR}/10.3_Figura_riqueza_resistoma.png",
    dpi=300,
    bbox_inches="tight"
)

fig.savefig(
    f"{OUT_DIR}/10.3_Figura_riqueza_resistoma.pdf",
    bbox_inches="tight"
)

plt.close()

print("Figura 1 generada correctamente.")
