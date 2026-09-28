import csv
import sys
from pathlib import Path

BASE = Path.home() / "TFM_resistoma"
KMA = BASE / "definitivo" / "kma"
META = BASE / "metadata" / "TFM_PRJEB52877_metadata_maestro.tsv"

IDENTIDAD_MINIMA = 90.0
COBERTURA_MINIMA = 80.0

if len(sys.argv) != 2 or sys.argv[1] not in {"1", "5", "10"}:
    print("Uso: python 9.4_matrices_resistoma.py {1|5|10}")
    sys.exit(1)

DEPTH_MINIMO = float(sys.argv[1])
SALIDA = BASE / "definitivo" / "matrices" / f"profundidad_{sys.argv[1]}"
SALIDA.mkdir(parents=True, exist_ok=True)


# 1. Leer metadata

metadata = {}

with open(META, encoding="utf-8-sig", newline="") as archivo:
    lector = csv.reader(archivo, delimiter="\t", quotechar='"')
    next(lector)

    for fila in lector:
        if len(fila) >= 7:
            metadata[fila[6].strip()] = {
                "muestra": fila[0].strip(),
                "grupo": fila[2].strip()
            }


# 2. Leer KMA y aplicar filtros

registros = []
valores = {}
aros = set()
runs = set()

for archivo in sorted(KMA.glob("*.res")):

    run = archivo.stem

    if run not in metadata:
        continue

    runs.add(run)

    with open(archivo, encoding="utf-8") as entrada:
        lector = csv.DictReader(entrada, delimiter="\t")

        for fila in lector:

            if (
                float(fila["Template_Identity"]) < IDENTIDAD_MINIMA
                or float(fila["Template_Coverage"]) < COBERTURA_MINIMA
                or float(fila["Depth"]) < DEPTH_MINIMO
            ):
                continue

            partes = fila["#Template"].split("|")

            if not partes or not partes[0].startswith("ARO:"):
                continue

            aro = partes[0]

            nombre = ""
            for parte in partes:
                if parte.startswith("Name:"):
                    nombre = parte.replace("Name:", "")
                    break

            registro = {
                "muestra": metadata[run]["muestra"],
                "run": run,
                "grupo": metadata[run]["grupo"],
                "aro": aro,
                "nombre": nombre,
                "depth": float(fila["Depth"]),
                "identidad": float(fila["Template_Identity"]),
                "cobertura": float(fila["Template_Coverage"])
            }

            registros.append(registro)
            valores[(run, aro)] = registro["depth"]
            aros.add(aro)


print(f"Depth mínimo: {DEPTH_MINIMO:g}")
print(f"RUNs en metadata: {len(metadata)}")
print(f"RUNs procesados: {len(runs)}")
print(f"Detecciones: {len(registros)}")
print(f"ARO únicos: {len(aros)}")


# 3. Tabla larga

with open(SALIDA / "resistoma_longitudinal.tsv", "w", encoding="utf-8", newline="") as salida:

    campos = [
        "muestra", "run", "grupo", "aro", "nombre",
        "depth", "identidad", "cobertura"
    ]

    escritor = csv.DictWriter(
        salida,
        fieldnames=campos,
        delimiter="\t",
        lineterminator="\n"
    )

    escritor.writeheader()
    escritor.writerows(registros)


# 4. Matriz Depth

aros = sorted(aros)

with open(SALIDA / "matriz_profundidad.tsv", "w", encoding="utf-8", newline="") as salida:

    escritor = csv.writer(salida, delimiter="\t", lineterminator="\n")
    escritor.writerow(["muestra", "run", "grupo"] + aros)

    for run in sorted(metadata):
        escritor.writerow(
            [
                metadata[run]["muestra"],
                run,
                metadata[run]["grupo"]
            ]
            + [valores.get((run, aro), 0) for aro in aros]
        )


# 5. Matriz presencia/ausencia

with open(SALIDA / "matriz_presencia_ausencia.tsv", "w", encoding="utf-8", newline="") as salida:

    escritor = csv.writer(salida, delimiter="\t", lineterminator="\n")
    escritor.writerow(["muestra", "run", "grupo"] + aros)

    for run in sorted(metadata):
        escritor.writerow(
            [
                metadata[run]["muestra"],
                run,
                metadata[run]["grupo"]
            ]
            + [1 if (run, aro) in valores else 0 for aro in aros]
        )


print(f"Salida: {SALIDA}")
