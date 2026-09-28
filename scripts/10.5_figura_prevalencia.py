#!/usr/bin/env python3

import os
import pandas as pd
import matplotlib.pyplot as plt


INPUT = "definitivo/resultados/aro_significativos_fdr05.tsv"
OUT_DIR = "definitivo/figuras"

os.makedirs(OUT_DIR, exist_ok=True)


# Leer resultados definitivos
df = pd.read_csv(INPUT, sep="\t")


# Nombres de los ARO
nombres = {
    "ARO:3003922": "oqxA",
    "ARO:3003923": "oqxB",
    "ARO:3004041": "acrA",
    "ARO:3004111": "FosA6",
    "ARO:3004122": "OmpK37",
    "ARO:3004580": "KpnE",
    "ARO:3004583": "KpnF",
    "ARO:3004588": "KpnG",
    "ARO:3005047": "eptB",
    "ARO:3005053": "ArnT",
    "ARO:3005059": "LptD",
    "ARO:3007682": "MdtQ"
}

df["nombre"] = df["ARO"].map(nombres)

# Ordenar de menor a mayor prevalencia
df = df.sort_values("prevalencia_colonizado")


# Figura
fig, ax = plt.subplots(figsize=(9, 6))

barras = ax.bar(
    df["nombre"],
    df["prevalencia_colonizado"]
)


# Porcentaje sobre cada barra
for barra, porcentaje in zip(
    barras,
    df["prevalencia_colonizado"]
):
    ax.text(
        barra.get_x() + barra.get_width() / 2,
        porcentaje + 0.8,
        f"{porcentaje:.1f}%",
        ha="center",
        va="bottom",
        fontsize=9
    )


ax.set_xlabel("Determinante ARO")
ax.set_ylabel("Prevalencia en muestras colonizadas (%)")

ax.set_title(
    "Prevalencia de los 12 determinantes ARO en muestras colonizadas"
)

ax.set_ylim(0, 45)

plt.xticks(rotation=45, ha="right")

plt.tight_layout()


# Guardar
fig.savefig(
    f"{OUT_DIR}/10.5_Figura_prevalencia_12_aro.png",
    dpi=300,
    bbox_inches="tight"
)

fig.savefig(
    f"{OUT_DIR}/10.5_Figura_prevalencia_12_aro.pdf",
    bbox_inches="tight"
)

plt.close()

print("Figura 2 generada correctamente.")
