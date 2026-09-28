#!/usr/bin/env python3

import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist


INPUT = "definitivo/matrices/profundidad_1/matriz_presencia_ausencia.tsv"
ARO_FILE = "definitivo/resultados/aro_significativos_fdr05.tsv"
OUT_DIR = "definitivo/figuras"

os.makedirs(OUT_DIR, exist_ok=True)


# Datos
df = pd.read_csv(INPUT, sep="\t")
aros = pd.read_csv(ARO_FILE, sep="\t")["ARO"].tolist()

df = df[df["grupo"] == "Colonizado"].copy()

matriz = df[aros].astype(float).values


# Clustering
link_muestras = linkage(
    pdist(matriz, metric="euclidean"),
    method="average"
)

link_aros = linkage(
    pdist(matriz.T, metric="euclidean"),
    method="average"
)

orden_muestras = dendrogram(
    link_muestras,
    no_plot=True
)["leaves"]

orden_aros = dendrogram(
    link_aros,
    no_plot=True
)["leaves"]


matriz = matriz[orden_muestras][:, orden_aros]


# Nombres ARO
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

labels_aros = [nombres[a] for a in [aros[i] for i in orden_aros]]


# Figura
fig = plt.figure(figsize=(12, 14))

ax1 = fig.add_axes([0.10, 0.78, 0.70, 0.15])
dendrogram(link_muestras, ax=ax1, no_labels=True)
ax1.axis("off")

ax2 = fig.add_axes([0.82, 0.10, 0.15, 0.68])
dendrogram(link_aros, ax=ax2, orientation="right", no_labels=True)
ax2.axis("off")

ax3 = fig.add_axes([0.10, 0.10, 0.70, 0.68])

im = ax3.imshow(
    matriz,
    aspect="auto",
    interpolation="nearest"
)

ax3.set_xticks(range(len(labels_aros)))
ax3.set_xticklabels(labels_aros, rotation=45, ha="right")

ax3.set_yticks(range(len(df)))
ax3.set_yticklabels(
    df.iloc[orden_muestras]["muestra"],
    fontsize=7
)

ax3.set_xlabel("Determinantes ARO")
ax3.set_ylabel("Muestras colonizadas")

fig.colorbar(
    im,
    ax=ax3,
    orientation="horizontal",
    pad=0.08,
    label="Presencia / ausencia"
)

fig.suptitle(
    "Determinantes ARO en muestras colonizadas",
    fontsize=14
)


# Guardar
plt.savefig(
    f"{OUT_DIR}/10.10_Figura_heatmap_12_aro.png",
    dpi=300,
    bbox_inches="tight"
)

plt.savefig(
    f"{OUT_DIR}/10.10_Figura_heatmap_12_aro.pdf",
    bbox_inches="tight"
)

plt.close()

print("Figura 5 generada correctamente.")
