#!/usr/bin/env python3

import os
import pandas as pd
import matplotlib.pyplot as plt


INPUT = "definitivo/resultados/caracterizacion_funcional_12_aro.tsv"
OUT_DIR = "definitivo/figuras"

os.makedirs(OUT_DIR, exist_ok=True)


# Leer resultados definitivos
df = pd.read_csv(INPUT, sep="\t")


# Contar ARO por mecanismo
mecanismos = df["Resistance_mechanism"].value_counts()
mecanismos.index = mecanismos.index.map({
    "antibiotic efflux": "Eflujo de antibióticos",
    "antibiotic inactivation": "Inactivación de antibióticos",
    "antibiotic target alteration": "Alteración de la diana",
    "reduced permeability to antibiotic": "Permeabilidad reducida"
})


# Figura
fig, ax = plt.subplots(figsize=(8, 5))

barras = ax.bar(
    mecanismos.index,
    mecanismos.values
)


# Número de ARO sobre cada barra
for barra, valor in zip(barras, mecanismos.values):
    ax.text(
        barra.get_x() + barra.get_width() / 2,
        valor + 0.1,
        str(valor),
        ha="center",
        va="bottom"
    )


ax.set_xlabel("Mecanismo de resistencia")
ax.set_ylabel("Número de determinantes ARO")

ax.set_title(
    "Mecanismos de resistencia de los 12 determinantes ARO"
)

plt.xticks(rotation=25, ha="right")

plt.tight_layout()


# Guardar
fig.savefig(
    f"{OUT_DIR}/10.6_Figura_mecanismos.png",
    dpi=300,
    bbox_inches="tight"
)

fig.savefig(
    f"{OUT_DIR}/10.6_Figura_mecanismos.pdf",
    bbox_inches="tight"
)

plt.close()

print("Figura 3 generada correctamente.")
print(mecanismos)
