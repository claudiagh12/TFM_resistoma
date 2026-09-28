import csv
import json
from pathlib import Path

BASE = Path.home() / "TFM_resistoma"
FASTP = BASE / "definitivo" / "control_calidad" / "fastp"
META = BASE / "metadata" / "TFM_PRJEB52877_metadata_maestro.tsv"
SALIDA = BASE / "definitivo" / "matrices" / "profundidad_secuenciacion.tsv"


# Leer metadata
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


resultados = []

for archivo in sorted(FASTP.glob("*_fastp.json")):

    run = archivo.name.replace("_fastp.json", "")

    if run not in metadata:
        continue

    with open(archivo, encoding="utf-8") as entrada:
        datos = json.load(entrada)

    reads_antes_r1 = datos["read1_before_filtering"]["total_reads"]
    reads_antes_r2 = datos["read2_before_filtering"]["total_reads"]

    reads_despues_r1 = datos["read1_after_filtering"]["total_reads"]
    reads_despues_r2 = datos["read2_after_filtering"]["total_reads"]

    pares_antes = min(reads_antes_r1, reads_antes_r2)
    pares_despues = min(reads_despues_r1, reads_despues_r2)

    resultados.append({
        "muestra": metadata[run]["muestra"],
        "run": run,
        "grupo": metadata[run]["grupo"],
        "reads_antes_r1": reads_antes_r1,
        "reads_antes_r2": reads_antes_r2,
        "reads_despues_r1": reads_despues_r1,
        "reads_despues_r2": reads_despues_r2,
        "pares_antes": pares_antes,
        "pares_despues": pares_despues,
        "retencion_porcentaje": (
            pares_despues / pares_antes
        ) * 100
    })


resultados.sort(key=lambda x: x["run"])


with open(SALIDA, "w", encoding="utf-8", newline="") as salida:

    campos = [
        "muestra",
        "run",
        "grupo",
        "reads_antes_r1",
        "reads_antes_r2",
        "reads_despues_r1",
        "reads_despues_r2",
        "pares_antes",
        "pares_despues",
        "retencion_porcentaje"
    ]

    escritor = csv.DictWriter(
        salida,
        fieldnames=campos,
        delimiter="\t",
        lineterminator="\n"
    )

    escritor.writeheader()
    escritor.writerows(resultados)


print(f"RUNs procesados: {len(resultados)}")
print(f"Archivo creado: {SALIDA}")
