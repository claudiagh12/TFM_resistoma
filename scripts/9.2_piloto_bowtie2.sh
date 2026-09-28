#!/bin/bash

set -u
set -o pipefail

# ============================================================
# TFM RESISTOMA
# Pipeline QC + eliminación de hospedador + KMA/CARD
#
# Flujo:
#
# ENA FASTQ
#    ↓
# Integridad gzip
#    ↓
# FastQC inicial
#    ↓
# fastp
#    ↓
# FastQC posterior
#    ↓
# Bowtie2 + GRCh38
#    ↓
# FASTQ no humanos
#    ↓
# KMA + CARD
#    ↓
# .res completo
#
# IMPORTANTE:
# No se aplica aquí un umbral de profundidad.
#
# Los análisis:
#   Depth >= 1
#   Depth >= 5
#   Depth >= 10
#
# se realizarán posteriormente mediante un script independiente
# de análisis de sensibilidad.
#
# PRUEBA PILOTO:
#   ERR9880676
#   ERR9880677
# ============================================================


# ============================================================
# CONFIGURACIÓN
# ============================================================

BASE="$HOME/TFM_resistoma"

META="$BASE/metadata/TFM_PRJEB52877_metadata_maestro.tsv"

WORK="/tmp/TFM_resistoma_qc_work"

QC_RAW="$BASE/qc/fastqc_raw_final"
QC_CLEAN="$BASE/qc/fastqc_clean_final"
FASTP_OUT="$BASE/qc/fastp_final"

HOST_OUT="$BASE/host_filter/final"

KMA_OUT="$BASE/results/kma_qc"

LOG_DIR="$BASE/logs/qc_kma"

LOG="$LOG_DIR/qc_host_kma_batch.log"


# ============================================================
# HERRAMIENTAS
# ============================================================

KMA="$HOME/miniconda3/envs/tfm_resistoma/bin/kma"

DB="$BASE/localDB/bwt/card_reference/kma"

HOST_INDEX="$BASE/host_filter/grch38/GRCh38_noalt_as"

FASTQC="$(which fastqc)"
FASTP="$(which fastp)"
BOWTIE2="$(which bowtie2)"

THREADS=4


# ============================================================
# MUESTRAS PILOTO
# ============================================================

RUNS=(
    "ERR9880676"
    "ERR9880677"
)


# ============================================================
# CREAR DIRECTORIOS
# ============================================================

mkdir -p \
    "$WORK" \
    "$QC_RAW" \
    "$QC_CLEAN" \
    "$FASTP_OUT" \
    "$HOST_OUT" \
    "$KMA_OUT" \
    "$LOG_DIR"


# ============================================================
# LOG PRINCIPAL
# ============================================================

echo "============================================================" | tee -a "$LOG"
echo " TFM RESISTOMA - QC + HOST FILTER + KMA" | tee -a "$LOG"
echo " PRUEBA PILOTO" | tee -a "$LOG"
echo " Inicio: $(date)" | tee -a "$LOG"
echo "============================================================" | tee -a "$LOG"
echo "" | tee -a "$LOG"


# ============================================================
# COMPROBAR HERRAMIENTAS
# ============================================================

echo "Comprobando herramientas..." | tee -a "$LOG"

if [ ! -x "$KMA" ]; then
    echo "ERROR: KMA no encontrado en:" | tee -a "$LOG"
    echo "$KMA" | tee -a "$LOG"
    exit 1
fi

if [ ! -x "$FASTQC" ]; then
    echo "ERROR: FastQC no encontrado." | tee -a "$LOG"
    exit 1
fi

if [ ! -x "$FASTP" ]; then
    echo "ERROR: fastp no encontrado." | tee -a "$LOG"
    exit 1
fi

if [ ! -x "$BOWTIE2" ]; then
    echo "ERROR: Bowtie2 no encontrado." | tee -a "$LOG"
    exit 1
fi

if [ ! -s "$DB" ] && [ ! -s "${DB}.index" ]; then
    echo "ADVERTENCIA: comprobar manualmente la base KMA:" | tee -a "$LOG"
    echo "$DB" | tee -a "$LOG"
fi


echo "" | tee -a "$LOG"
echo "Versiones:" | tee -a "$LOG"

"$FASTQC" --version | tee -a "$LOG"

"$FASTP" --version | tee -a "$LOG"

"$BOWTIE2" --version 2>&1 | head -4 | tee -a "$LOG"

"$KMA" -v | tee -a "$LOG"

echo "" | tee -a "$LOG"


# ============================================================
# COMPROBAR ÍNDICE HUMANO GRCh38
# ============================================================

echo "Comprobando índice GRCh38..." | tee -a "$LOG"

if [ ! -s "${HOST_INDEX}.1.bt2" ] || \
   [ ! -s "${HOST_INDEX}.2.bt2" ] || \
   [ ! -s "${HOST_INDEX}.3.bt2" ] || \
   [ ! -s "${HOST_INDEX}.4.bt2" ] || \
   [ ! -s "${HOST_INDEX}.rev.1.bt2" ] || \
   [ ! -s "${HOST_INDEX}.rev.2.bt2" ]; then

    echo "ERROR: índice GRCh38 incompleto." | tee -a "$LOG"
    exit 1
fi

echo "Índice GRCh38 OK." | tee -a "$LOG"
echo "" | tee -a "$LOG"


# ============================================================
# PROCESAMIENTO DE CADA RUN
# ============================================================

for RUN in "${RUNS[@]}"; do

    echo "============================================================" | tee -a "$LOG"
    echo "Procesando RUN: $RUN" | tee -a "$LOG"
    echo "============================================================" | tee -a "$LOG"


    # --------------------------------------------------------
    # DIRECTORIOS TEMPORALES
    # --------------------------------------------------------

    SAMPLE_DIR="$WORK/$RUN"

    RAW_DIR="$SAMPLE_DIR/raw"
    CLEAN_DIR="$SAMPLE_DIR/clean"
    HOST_DIR="$SAMPLE_DIR/host"

    mkdir -p \
        "$RAW_DIR" \
        "$CLEAN_DIR" \
        "$HOST_DIR"


    # --------------------------------------------------------
    # ARCHIVOS FASTQ
    # --------------------------------------------------------

    R1="$RAW_DIR/${RUN}_1.fastq.gz"
    R2="$RAW_DIR/${RUN}_2.fastq.gz"

    C1="$CLEAN_DIR/${RUN}_1.fastq.gz"
    C2="$CLEAN_DIR/${RUN}_2.fastq.gz"

    NONHUMAN_PREFIX="$HOST_DIR/${RUN}_nonhuman"

    NONHUMAN_R1="${NONHUMAN_PREFIX}.1"
    NONHUMAN_R2="${NONHUMAN_PREFIX}.2"


    # --------------------------------------------------------
    # OBTENER METADATA
    # --------------------------------------------------------

    LINE=$(awk -F'\t' -v run="$RUN" '$7 == run {print; exit}' "$META")

    if [ -z "$LINE" ]; then

        echo "ERROR: RUN $RUN no encontrado en metadata." | tee -a "$LOG"

        rm -rf "$SAMPLE_DIR"

        continue

    fi


    FASTQ_URLS=$(echo "$LINE" | cut -f17)

    URL1=$(echo "$FASTQ_URLS" | cut -d';' -f1)
    URL2=$(echo "$FASTQ_URLS" | cut -d';' -f2)


    # ========================================================
    # DESCARGA R1
    # ========================================================

    if [ ! -s "$R1" ]; then

        echo "Descargando R1..." | tee -a "$LOG"

        wget \
          --tries=5 \
          --timeout=60 \
          --continue \
          --show-progress \
          "https://$URL1" \
          -O "$R1.part"

        if [ $? -ne 0 ]; then

            echo "ERROR descargando R1 -> $RUN" | tee -a "$LOG"

            rm -rf "$SAMPLE_DIR"

            continue

        fi

        mv "$R1.part" "$R1"

    fi


    # ========================================================
    # DESCARGA R2
    # ========================================================

    if [ ! -s "$R2" ]; then

        echo "Descargando R2..." | tee -a "$LOG"

        wget \
          --tries=5 \
          --timeout=60 \
          --continue \
          --show-progress \
          "https://$URL2" \
          -O "$R2.part"

        if [ $? -ne 0 ]; then

            echo "ERROR descargando R2 -> $RUN" | tee -a "$LOG"

            rm -rf "$SAMPLE_DIR"

            continue

        fi

        mv "$R2.part" "$R2"

    fi


    # ========================================================
    # COMPROBAR FASTQ
    # ========================================================

    if [ ! -s "$R1" ] || [ ! -s "$R2" ]; then

        echo "ERROR: FASTQ incompletos -> $RUN" | tee -a "$LOG"

        rm -rf "$SAMPLE_DIR"

        continue

    fi

    echo "FASTQ descargados correctamente." | tee -a "$LOG"


    # ========================================================
    # INTEGRIDAD GZIP
    # ========================================================

    echo "Comprobando integridad gzip..." | tee -a "$LOG"

    if ! gzip -t "$R1"; then

        echo "ERROR: R1 corrupto -> $RUN" | tee -a "$LOG"

        rm -rf "$SAMPLE_DIR"

        continue

    fi


    if ! gzip -t "$R2"; then

        echo "ERROR: R2 corrupto -> $RUN" | tee -a "$LOG"

        rm -rf "$SAMPLE_DIR"

        continue

    fi


    echo "Integridad gzip OK." | tee -a "$LOG"


    # ========================================================
    # FASTQC INICIAL
    # ========================================================

    echo "FastQC inicial..." | tee -a "$LOG"

    "$FASTQC" \
      -t "$THREADS" \
      -o "$QC_RAW" \
      "$R1" "$R2" \
      >> "$LOG" 2>&1

    if [ $? -ne 0 ]; then

        echo "ERROR: FastQC inicial falló -> $RUN" | tee -a "$LOG"

        rm -rf "$SAMPLE_DIR"

        continue

    fi


    # ========================================================
    # FASTP
    # ========================================================

    echo "Ejecutando fastp..." | tee -a "$LOG"

    "$FASTP" \
      -i "$R1" \
      -I "$R2" \
      -o "$C1" \
      -O "$C2" \
      -h "$FASTP_OUT/${RUN}_fastp.html" \
      -j "$FASTP_OUT/${RUN}_fastp.json" \
      --detect_adapter_for_pe \
      --cut_front \
      --cut_tail \
      --cut_mean_quality 20 \
      --length_required 50 \
      --thread "$THREADS" \
      >> "$LOG" 2>&1

    if [ $? -ne 0 ]; then

        echo "ERROR: fastp falló -> $RUN" | tee -a "$LOG"

        rm -rf "$SAMPLE_DIR"

        continue

    fi


    # ========================================================
    # COMPROBAR FASTQ POST-FASTP
    # ========================================================

    if [ ! -s "$C1" ] || [ ! -s "$C2" ]; then

        echo "ERROR: FASTQ post-fastp inexistentes -> $RUN" | tee -a "$LOG"

        rm -rf "$SAMPLE_DIR"

        continue

    fi


    if ! gzip -t "$C1" || ! gzip -t "$C2"; then

        echo "ERROR: FASTQ post-fastp corruptos -> $RUN" | tee -a "$LOG"

        rm -rf "$SAMPLE_DIR"

        continue

    fi


    # ========================================================
    # FASTQC POST-FASTP
    # ========================================================

    echo "FastQC posterior a fastp..." | tee -a "$LOG"

    "$FASTQC" \
      -t "$THREADS" \
      -o "$QC_CLEAN" \
      "$C1" "$C2" \
      >> "$LOG" 2>&1

    if [ $? -ne 0 ]; then

        echo "ERROR: FastQC posterior falló -> $RUN" | tee -a "$LOG"

        rm -rf "$SAMPLE_DIR"

        continue

    fi


    # ========================================================
    # CONTEO DE PARES POST-FASTP
    # ========================================================

    CLEAN_PAIRS=$(zcat "$C1" | awk 'END {print NR/4}')

    echo "Pares después de fastp: $CLEAN_PAIRS" | tee -a "$LOG"


    # ========================================================
    # ELIMINACIÓN DE HOSPEDADOR
    # ========================================================

    echo "Eliminación de lecturas humanas con Bowtie2..." | tee -a "$LOG"


    "$BOWTIE2" \
      -x "$HOST_INDEX" \
      -1 "$C1" \
      -2 "$C2" \
      --very-sensitive \
      --threads "$THREADS" \
      --un-conc-gz "$NONHUMAN_PREFIX" \
      -S /dev/null \
      2> "$LOG_DIR/${RUN}_bowtie2.log"


    BOWTIE_STATUS=$?


    if [ $BOWTIE_STATUS -ne 0 ]; then

        echo "ERROR: Bowtie2 falló -> $RUN" | tee -a "$LOG"

        rm -rf "$SAMPLE_DIR"

        continue

    fi


    # ========================================================
    # COMPROBAR FASTQ NO HUMANOS
    # ========================================================

    if [ ! -s "$NONHUMAN_R1" ] || [ ! -s "$NONHUMAN_R2" ]; then

        echo "ERROR: FASTQ no humanos inexistentes -> $RUN" | tee -a "$LOG"

        rm -rf "$SAMPLE_DIR"

        continue

    fi


    if ! gzip -t "$NONHUMAN_R1" || ! gzip -t "$NONHUMAN_R2"; then

        echo "ERROR: FASTQ no humanos corruptos -> $RUN" | tee -a "$LOG"

        rm -rf "$SAMPLE_DIR"

        continue

    fi


    NONHUMAN_PAIRS=$(zcat "$NONHUMAN_R1" | awk 'END {print NR/4}')

    echo "Pares no humanos: $NONHUMAN_PAIRS" | tee -a "$LOG"


    # ========================================================
    # KMA + CARD
    # ========================================================

    echo "Ejecutando KMA + CARD..." | tee -a "$LOG"


    "$KMA" \
      -mem_mode \
      -ex_mode \
      -1t1 \
      -vcf \
      -ipe "$NONHUMAN_R1" "$NONHUMAN_R2" \
      -t 2 \
      -t_db "$DB" \
      -o "$KMA_OUT/$RUN" \
      2>&1 | tee "$LOG_DIR/${RUN}_kma.log"


    KMA_STATUS=${PIPESTATUS[0]}


    if [ $KMA_STATUS -ne 0 ]; then

        echo "ERROR: KMA falló -> $RUN" | tee -a "$LOG"

        rm -rf "$SAMPLE_DIR"

        continue

    fi


    echo "KMA terminado correctamente." | tee -a "$LOG"


    # ========================================================
    # COMPROBAR RESULTADO KMA
    # ========================================================

    RES="$KMA_OUT/${RUN}.res"


    if [ ! -s "$RES" ]; then

        echo "ERROR: no se generó $RES" | tee -a "$LOG"

        rm -rf "$SAMPLE_DIR"

        continue

    fi


    # ========================================================
    # CONSERVAR .RES COMPLETO
    #
    # NO SE APLICA FILTRO DE DEPTH AQUÍ.
    #
    # Posteriormente se evaluarán:
    #   Depth >= 1
    #   Depth >= 5
    #   Depth >= 10
    #
    # manteniendo:
    #   Template_Identity >= 90 %
    #   Template_Coverage >= 80 %
    # ========================================================

    KMA_HITS=$(tail -n +2 "$RES" | wc -l)

    echo "Detecciones reportadas por KMA: $KMA_HITS" | tee -a "$LOG"


    # ========================================================
    # ELIMINAR ARCHIVOS AUXILIARES DE KMA
    #
    # Se conserva únicamente:
    #   .res
    # ========================================================

    rm -f \
      "$KMA_OUT/${RUN}.aln" \
      "$KMA_OUT/${RUN}.frag.gz" \
      "$KMA_OUT/${RUN}.fsa" \
      "$KMA_OUT/${RUN}.vcf.gz"


    # ========================================================
    # RESUMEN DE PROCESAMIENTO
    # ========================================================

    SUMMARY="$HOST_OUT/${RUN}_processing_summary.tsv"


    {
        echo -e "RUN\tPairs_post_fastp\tPairs_nonhuman\tKMA_hits"
        echo -e "$RUN\t$CLEAN_PAIRS\t$NONHUMAN_PAIRS\t$KMA_HITS"
    } > "$SUMMARY"


    # ========================================================
    # ELIMINAR FASTQ TEMPORALES
    # ========================================================

    echo "Eliminando FASTQ temporales..." | tee -a "$LOG"

    rm -rf "$SAMPLE_DIR"


    echo "Muestra $RUN completada correctamente." | tee -a "$LOG"
    echo "" | tee -a "$LOG"

done


# ============================================================
# FIN
# ============================================================

echo "============================================================" | tee -a "$LOG"
echo " PROCESAMIENTO PILOTO FINALIZADO" | tee -a "$LOG"
echo " Fin: $(date)" | tee -a "$LOG"
echo "============================================================" | tee -a "$LOG"
echo "" | tee -a "$LOG"

echo "Resultados KMA:" | tee -a "$LOG"

ls -lh "$KMA_OUT" | tee -a "$LOG"
