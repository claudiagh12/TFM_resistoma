#!/usr/bin/env bash

set -u
set -o pipefail

# ============================================================
# TFM RESISTOMA
# PROCESAMIENTO DEFINITIVO: QC + KMA + CARD
#
# FASTQ -> FastQC -> fastp -> FastQC -> KMA
#
# No se realiza filtrado de hospedador en este pipeline.
# El filtrado con Bowtie2 se evaluó previamente como piloto.
# ============================================================

BASE="$HOME/TFM_resistoma"

META="$BASE/metadata/TFM_PRJEB52877_metadata_maestro.tsv"

WORK="/tmp/TFM_resistoma_kma_work"

QC_RAW="$BASE/qc/fastqc_raw_final"
QC_CLEAN="$BASE/qc/fastqc_clean_final"
FASTP_OUT="$BASE/qc/fastp_final"

KMA_OUT="$BASE/results/kma_qc_final"
LOG_DIR="$BASE/logs/kma_final"

KMA="$HOME/miniconda3/envs/tfm_resistoma/bin/kma"
DB="$BASE/localDB/bwt/card_reference/kma"

FASTQC="$(which fastqc)"
FASTP="$(which fastp)"
WGET="$(which wget)"

THREADS=2

mkdir -p \
    "$WORK" \
    "$QC_RAW" \
    "$QC_CLEAN" \
    "$FASTP_OUT" \
    "$KMA_OUT" \
    "$LOG_DIR"

LOG="$LOG_DIR/kma_batch.log"

echo "============================================================" | tee -a "$LOG"
echo "TFM RESISTOMA - PROCESAMIENTO DEFINITIVO QC + KMA" | tee -a "$LOG"
echo "Fecha: $(date)" | tee -a "$LOG"
echo "KMA: $KMA" | tee -a "$LOG"
echo "Threads: $THREADS" | tee -a "$LOG"
echo "============================================================" | tee -a "$LOG"

# ------------------------------------------------------------
# Comprobación de herramientas
# ------------------------------------------------------------

for TOOL in "$FASTQC" "$FASTP" "$WGET" "$KMA"; do
    if [ ! -x "$TOOL" ]; then
        echo "ERROR: herramienta no encontrada o no ejecutable: $TOOL" | tee -a "$LOG"
        exit 1
    fi
done

if [ ! -f "$META" ]; then
    echo "ERROR: metadata no encontrada: $META" | tee -a "$LOG"
    exit 1
fi

if [ ! -d "$DB" ]; then
    echo "ERROR: base KMA no encontrada: $DB" | tee -a "$LOG"
    exit 1
fi

echo "FastQC: $("$FASTQC" --version)" | tee -a "$LOG"
echo "fastp: $("$FASTP" --version)" | tee -a "$LOG"
echo "KMA: $("$KMA" -v 2>&1 | head -n 1)" | tee -a "$LOG"

# ------------------------------------------------------------
# Obtener RUNs desde metadata
# ------------------------------------------------------------

mapfile -t RUNS < <(
    python - "$META" <<'PY'
import csv
import sys

metadata = sys.argv[1]

with open(metadata, newline="", encoding="utf-8") as f:
    reader = csv.reader(f, delimiter="\t", quotechar='"')
    header = next(reader)

    run_idx = header.index("run_accession")

    for row in reader:
        if len(row) > run_idx and row[run_idx].strip():
            print(row[run_idx].strip())
PY
)

echo "Número de RUN a procesar: ${#RUNS[@]}" | tee -a "$LOG"

# ------------------------------------------------------------
# Procesamiento
# ------------------------------------------------------------

for RUN in "${RUNS[@]}"; do

    echo "" | tee -a "$LOG"
    echo "============================================================" | tee -a "$LOG"
    echo "Procesando RUN: $RUN" | tee -a "$LOG"
    echo "Inicio: $(date)" | tee -a "$LOG"
    echo "============================================================" | tee -a "$LOG"

    SAMPLE_DIR="$WORK/$RUN"
    RAW_DIR="$SAMPLE_DIR/raw"
    CLEAN_DIR="$SAMPLE_DIR/clean"

    mkdir -p "$RAW_DIR" "$CLEAN_DIR"

    # --------------------------------------------------------
    # Si KMA ya terminó correctamente, no repetir muestra
    # --------------------------------------------------------

    if [ -s "$KMA_OUT/$RUN.res" ]; then
        echo "KMA ya existe para $RUN. Se omite esta muestra." | tee -a "$LOG"
        continue
    fi

    # --------------------------------------------------------
    # Obtener URLs FASTQ desde metadata
    # --------------------------------------------------------

    URLS=$(python - "$META" "$RUN" <<'PY'
import csv
import sys

metadata = sys.argv[1]
target = sys.argv[2]

with open(metadata, newline="", encoding="utf-8") as f:
    reader = csv.reader(f, delimiter="\t", quotechar='"')
    header = next(reader)

    run_idx = header.index("run_accession")
    ftp_idx = header.index("fastq_ftp")

    for row in reader:
        if len(row) > run_idx and row[run_idx].strip() == target:
            print(row[ftp_idx].strip())
            break
PY
)

    if [ -z "$URLS" ]; then
        echo "ERROR: no se encontraron FASTQ para $RUN" | tee -a "$LOG"
        continue
    fi

    # --------------------------------------------------------
    # Descargar FASTQ
    # --------------------------------------------------------

    echo "Descargando FASTQ..." | tee -a "$LOG"

    R1_URL=$(echo "$URLS" | cut -d ";" -f 1)
    R2_URL=$(echo "$URLS" | cut -d ";" -f 2)
    R1_URL="https://$R1_URL"
    R2_URL="https://$R2_URL"

    R1="$RAW_DIR/${RUN}_1.fastq.gz"
    R2="$RAW_DIR/${RUN}_2.fastq.gz"

    "$WGET" -c -O "$R1" "$R1_URL" 2>&1 | tee "$LOG_DIR/${RUN}_download_R1.log"

    "$WGET" -c -O "$R2" "$R2_URL" 2>&1 | tee "$LOG_DIR/${RUN}_download_R2.log"

    if [ ! -s "$R1" ] || [ ! -s "$R2" ]; then
        echo "ERROR: FASTQ incompletos para $RUN" | tee -a "$LOG"
        rm -rf "$SAMPLE_DIR"
        continue
    fi

    echo "FASTQ descargados correctamente." | tee -a "$LOG"

    # --------------------------------------------------------
    # Integridad gzip
    # --------------------------------------------------------

    echo "Comprobando integridad gzip..." | tee -a "$LOG"

    if ! gzip -t "$R1" "$R2"; then
        echo "ERROR: fallo de integridad gzip en $RUN" | tee -a "$LOG"
        rm -rf "$SAMPLE_DIR"
        continue
    fi

    echo "Integridad gzip OK." | tee -a "$LOG"

    # --------------------------------------------------------
    # FastQC inicial
    # --------------------------------------------------------

    echo "FastQC inicial..." | tee -a "$LOG"

    "$FASTQC" \
        --threads "$THREADS" \
        --outdir "$QC_RAW" \
        "$R1" "$R2" \
        2>&1 | tee "$LOG_DIR/${RUN}_fastqc_raw.log"

    # --------------------------------------------------------
    # fastp
    # --------------------------------------------------------

    echo "Ejecutando fastp..." | tee -a "$LOG"

    CLEAN_R1="$CLEAN_DIR/${RUN}_1.fastq.gz"
    CLEAN_R2="$CLEAN_DIR/${RUN}_2.fastq.gz"

    "$FASTP" \
        -i "$R1" \
        -I "$R2" \
        -o "$CLEAN_R1" \
        -O "$CLEAN_R2" \
        --detect_adapter_for_pe \
        --cut_front \
        --cut_tail \
        --cut_mean_quality 20 \
        --length_required 50 \
        --thread "$THREADS" \
        -h "$FASTP_OUT/${RUN}_fastp.html" \
        -j "$FASTP_OUT/${RUN}_fastp.json" \
        2>&1 | tee "$LOG_DIR/${RUN}_fastp.log"

    if [ ! -s "$CLEAN_R1" ] || [ ! -s "$CLEAN_R2" ]; then
        echo "ERROR: fastp no produjo FASTQ para $RUN" | tee -a "$LOG"
        rm -rf "$SAMPLE_DIR"
        continue
    fi

    # --------------------------------------------------------
    # FastQC posterior
    # --------------------------------------------------------

    echo "FastQC posterior a fastp..." | tee -a "$LOG"

    "$FASTQC" \
        --threads "$THREADS" \
        --outdir "$QC_CLEAN" \
        "$CLEAN_R1" "$CLEAN_R2" \
        2>&1 | tee "$LOG_DIR/${RUN}_fastqc_clean.log"

    # --------------------------------------------------------
    # Contar pares después de fastp
    # --------------------------------------------------------

    READS_R1=$(gzip -cd "$CLEAN_R1" | wc -l)
    PAIRS=$((READS_R1 / 4))

    echo "Pares después de fastp: $PAIRS" | tee -a "$LOG"

    # --------------------------------------------------------
    # KMA + CARD
    # --------------------------------------------------------

    echo "Ejecutando KMA + CARD..." | tee -a "$LOG"

    "$KMA" \
        -mem_mode \
        -ex_mode \
        -1t1 \
        -vcf \
        -ipe "$CLEAN_R1" "$CLEAN_R2" \
        -t "$THREADS" \
        -t_db "$DB" \
        -o "$KMA_OUT/$RUN" \
        2>&1 | tee "$LOG_DIR/${RUN}_kma.log"

    KMA_STATUS=${PIPESTATUS[0]}

    if [ "$KMA_STATUS" -ne 0 ]; then
        echo "ERROR: KMA falló para $RUN" | tee -a "$LOG"
        rm -rf "$SAMPLE_DIR"
        continue
    fi

    if [ ! -s "$KMA_OUT/$RUN.res" ]; then
        echo "ERROR: KMA terminó pero no produjo .res para $RUN" | tee -a "$LOG"
        rm -rf "$SAMPLE_DIR"
        continue
    fi

    KMA_HITS=$(tail -n +2 "$KMA_OUT/$RUN.res" | wc -l)

    echo "KMA terminado correctamente." | tee -a "$LOG"
    echo "Detecciones reportadas por KMA: $KMA_HITS" | tee -a "$LOG"

    # --------------------------------------------------------
    # Resumen de procesamiento
    # --------------------------------------------------------

    cat > "$KMA_OUT/${RUN}_processing_summary.tsv" <<EOF
RUN	Pairs_post_fastp	KMA_hits
$RUN	$PAIRS	$KMA_HITS
EOF

    # --------------------------------------------------------
    # Eliminar FASTQ temporales
    # --------------------------------------------------------

    rm -rf "$SAMPLE_DIR"

    echo "Muestra $RUN completada correctamente." | tee -a "$LOG"
    echo "Fin: $(date)" | tee -a "$LOG"

done

echo "" | tee -a "$LOG"
echo "============================================================" | tee -a "$LOG"
echo "PROCESAMIENTO DEFINITIVO TERMINADO" | tee -a "$LOG"
echo "Fin: $(date)" | tee -a "$LOG"
echo "============================================================" | tee -a "$LOG"
