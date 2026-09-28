from pathlib import Path
import csv
import re

BASE = Path.home() / "TFM_resistoma"

# Fuentes
ARO_SIGNIFICATIVOS = BASE / "definitivo/resultados/aro_significativos_fdr05.tsv"
CARD = BASE / "historico/resultados_antiguos/resistome/CARD_annotation_significant_12.tsv"

# Salida
OUT = BASE / "definitivo/resultados/caracterizacion_funcional_12_aro.tsv"

# Campos funcionales que se conservarán
CARD_FIELDS = [
    "ARO",
    "Name",
    "CARD_short_name",
    "ARO_description",
    "AMR_gene_family",
    "Antibiotics",
    "Drug_classes",
    "Resistance_mechanism",
    "Model_type",
]


def read_significant_aros(path):
    ar_os = []

    with path.open(encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)

        for row in reader:
            if row and row[0].startswith("ARO:"):
                ar_os.append(row[0])

    return ar_os


def read_card(path):
    records = {}

    with path.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:
            aro = row["ARO"].strip()

            if aro:
                records[aro] = {
                    field: row.get(field, "").strip()
                    for field in CARD_FIELDS
                }

    return records


def main():
    significant = read_significant_aros(ARO_SIGNIFICATIVOS)
    card = read_card(CARD)

    significant_set = set(significant)
    card_set = set(card)

    missing_card = significant_set - card_set
    extra_card = card_set - significant_set

    print("===== VALIDACIÓN 9.8 =====")
    print(f"ARO significativos: {len(significant_set)}")
    print(f"ARO anotados en CARD: {len(card_set)}")

    if missing_card:
        print("\nERROR: ARO significativos sin anotación CARD:")
        for aro in sorted(missing_card):
            print(aro)
        raise SystemExit(1)

    if extra_card:
        print("\nAVISO: ARO presentes en CARD pero no en resultados definitivos:")
        for aro in sorted(extra_card):
            print(aro)

    # Orden: el mismo orden de los resultados definitivos
    records = [card[aro] for aro in significant if aro in card]

    OUT.parent.mkdir(parents=True, exist_ok=True)

    with OUT.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=CARD_FIELDS,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(records)

    print(f"\nRegistros funcionales generados: {len(records)}")
    print(f"Archivo: {OUT}")

    # Comprobación final
    if len(records) != len(significant_set):
        raise SystemExit(
            "ERROR: el número de registros generados no coincide "
            "con los ARO significativos."
        )

    print("\n===== VALIDACIÓN FINAL =====")
    print("OK: los 12 ARO tienen anotación funcional CARD.")
    print("OK: no se incorporaron variables estadísticas.")
    print("OK: no se modificó el archivo CARD histórico.")


if __name__ == "__main__":
    main()
