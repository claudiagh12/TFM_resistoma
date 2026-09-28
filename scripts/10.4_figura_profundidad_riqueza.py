import csv
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from scipy.stats import spearmanr
from pathlib import Path

BASE = Path.home() / "TFM_resistoma"

profundidad_file = BASE / "definitivo/matrices/profundidad_secuenciacion.tsv"
matriz_file = BASE / "definitivo/matrices/profundidad_1/matriz_presencia_ausencia.tsv"
salida = BASE / "definitivo/figuras"

salida.mkdir(parents=True, exist_ok=True)

# Profundidad de secuenciación
profundidad = {}

with open(profundidad_file, encoding="utf-8") as archivo:
    for fila in csv.DictReader(archivo, delimiter="\t"):
        profundidad[fila["run"]] = int(fila["pares_despues"])

# Riqueza y grupo
x = []
y = []
grupos = []

with open(matriz_file, encoding="utf-8") as archivo:
    for fila in csv.DictReader(archivo, delimiter="\t"):

        riqueza = sum(
            int(valor)
            for nombre, valor in fila.items()
            if nombre not in {"muestra", "run", "grupo"}
        )

        x.append(profundidad[fila["run"]])
        y.append(riqueza)
        grupos.append(fila["grupo"])

# Correlación
rho, p = spearmanr(x, y)

# Figura
fig, ax = plt.subplots(figsize=(8, 6))

for grupo, marcador in [
    ("Colonizado", "o"),
    ("No colonizado", "s")
]:

    datos_x = [x[i] for i in range(len(x)) if grupos[i] == grupo]
    datos_y = [y[i] for i in range(len(y)) if grupos[i] == grupo]

    ax.scatter(
        datos_x,
        datos_y,
        marker=marcador,
        label=grupo,
        alpha=0.75
    )

# Eje X en millones de pares
ax.xaxis.set_major_formatter(
    FuncFormatter(lambda valor, posicion: f"{valor / 1e6:.0f}")
)

ax.set_xlabel(
    "Profundidad de secuenciación (millones de pares retenidos)"
)

ax.set_ylabel("Riqueza del resistoma (número de ARO)")

ax.set_title(
    "Relación entre profundidad de secuenciación y riqueza del resistoma"
)

ax.text(
    0.05,
    0.95,
    f"Spearman ρ = {rho:.4f}\np = {p:.6g}",
    transform=ax.transAxes,
    va="top"
)

ax.legend()
ax.grid(alpha=0.2)

fig.tight_layout()

fig.savefig(
    salida / "10.4_Figura_profundidad_riqueza.png",
    dpi=300,
    bbox_inches="tight"
)

fig.savefig(
    salida / "10.4_Figura_profundidad_riqueza.pdf",
    bbox_inches="tight"
)

plt.close()

print("Figura 10.4 generada correctamente.")
print(f"n = {len(y)}")
print(f"Spearman rho = {rho:.4f}")
print(f"p = {p:.6g}")
