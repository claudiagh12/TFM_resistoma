import csv
import numpy as np
from scipy.stats import mannwhitneyu

ARCHIVO = "definitivo/matrices/profundidad_secuenciacion.tsv"

datos = []

with open(ARCHIVO, encoding="utf-8") as archivo:
    lector = csv.DictReader(archivo, delimiter="\t")

    for fila in lector:
        datos.append({
            "grupo": fila["grupo"],
            "pares": int(fila["pares_despues"])
        })


colonizado = np.array([
    x["pares"] for x in datos
    if x["grupo"] == "Colonizado"
])

no_colonizado = np.array([
    x["pares"] for x in datos
    if x["grupo"] == "No colonizado"
])


def resumen(x):
    return {
        "n": len(x),
        "mediana": np.median(x),
        "min": np.min(x),
        "max": np.max(x)
    }


c = resumen(colonizado)
nc = resumen(no_colonizado)

resultado = mannwhitneyu(
    colonizado,
    no_colonizado,
    alternative="two-sided"
)

z = (resultado.statistic - (len(colonizado) * len(no_colonizado) / 2)) / np.sqrt(
    len(colonizado) * len(no_colonizado) * (len(colonizado) + len(no_colonizado) + 1) / 12
)

r = abs(z) / np.sqrt(len(datos))


print("PROFUNDIDAD DE SECUENCIACIÓN POST-fastp")
print()
print("Colonizado:")
print(f"  n = {c['n']}")
print(f"  mediana = {c['mediana']:,.0f} pares")
print(f"  rango = {c['min']:,.0f} – {c['max']:,.0f}")

print()
print("No colonizado:")
print(f"  n = {nc['n']}")
print(f"  mediana = {nc['mediana']:,.0f} pares")
print(f"  rango = {nc['min']:,.0f} – {nc['max']:,.0f}")

print()
print("Mann–Whitney:")
print(f"  U = {resultado.statistic:.1f}")
print(f"  p = {resultado.pvalue:.6g}")
print(f"  r = {r:.4f}")
