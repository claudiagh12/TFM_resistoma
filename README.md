# TFM – Resistoma intestinal asociado a *Klebsiella pneumoniae*

## Descripción

Este repositorio contiene los scripts y resultados correspondientes al análisis
bioinformático realizado en el Trabajo Fin de Máster (TFM) en Bioinformática.

El estudio caracteriza el resistoma intestinal en relación con la
colonización por *Klebsiella pneumoniae* mediante datos de metagenómica
shotgun.

## Dataset

El análisis se realizó a partir del proyecto público de ENA:

- **Proyecto:** PRJEB52877
- **Estudio:** ERP137622
- **Muestras:** 103
- **Colonizados:** 52
- **No colonizados:** 51
- **Tecnología:** Illumina NovaSeq
- **Secuenciación:** paired-end shotgun metagenomics

Los datos de secuenciación originales no se incluyen en este repositorio y
están disponibles públicamente en ENA mediante el proyecto PRJEB52877.

## Flujo de análisis

El análisis comprende las siguientes etapas:

1. **Control de calidad:** FastQC y fastp antes y después del procesamiento.
2. **Identificación de resistencia:** KMA utilizando CARD v4.0.1.
3. **Construcción del resistoma:** matrices para profundidades ≥1, ≥5 y ≥10.
4. **Análisis de riqueza y profundidad:** comparación entre grupos, correlación
   y regresión binomial negativa.
5. **Análisis de asociación:** prueba exacta de Fisher y corrección de
   Benjamini-Hochberg.
6. **Caracterización funcional:** mecanismos de resistencia y clases de
   antibióticos.
7. **Coocurrencia:** análisis descriptivo de los 12 determinantes asociados.
8. **Análisis de sensibilidad:** exclusión de las 11 muestras
   cultivo-negativas/qPCR-positivas.
9. **Visualización:** generación de las figuras finales del TFM.

De forma exploratoria, se realizó un análisis piloto mediante Bowtie2 frente a
GRCh38 para evaluar la posible presencia de secuencias de origen humano.
Este procedimiento no constituyó un filtrado sistemático de las 103 muestras.

Los criterios principales de detección mediante KMA fueron:

- identidad ≥ 90 %;
- cobertura ≥ 80 %;
- profundidad ≥ 1.

## Resultados principales

Se identificaron **266 determinantes ARO** en las 103 muestras.

La riqueza del resistoma fue:

- **Colonizados:** mediana 61,5 ARO/muestra.
- **No colonizados:** mediana 43 ARO/muestra.
- **Mann–Whitney U = 1818,5; p = 0,00117.**

Se identificaron **12 determinantes ARO asociados al estado de colonización**:

- oqxA
- oqxB
- *K. pneumoniae* acrA
- FosA6
- *K. pneumoniae* OmpK37
- *K. pneumoniae* KpnE
- *K. pneumoniae* KpnF
- *K. pneumoniae* KpnG
- eptB
- ArnT
- LptD
- MdtQ

## Figuras

Las figuras finales generadas mediante scripts incluyen:

- **Sección 10.3:** riqueza del resistoma según el estado de colonización.
- **Sección 10.4:** relación entre profundidad de secuenciación y riqueza del resistoma.
- **Sección 10.5:** prevalencia de los 12 determinantes ARO asociados.
- **Sección 10.6:** mecanismos de resistencia de los 12 determinantes ARO.
- **Sección 10.8:** relación entre el número de determinantes ARO y la carga por qPCR.
- **Sección 10.10:** patrón de presencia/ausencia de los 12 determinantes ARO
  en las muestras colonizadas mediante mapa de calor y agrupamiento jerárquico.

Las figuras se encuentran en:

```text
definitivo/figuras/
```

## Scripts

Los scripts se encuentran organizados según las etapas metodológicas del TFM:

```text
scripts/
├── 9.2-9.3_procesamiento_qc_kma.sh
├── 9.2_piloto_bowtie2.sh
├── 9.4_matrices_resistoma.py
├── 9.5_*.py
├── 9.6_*.py
├── 9.7_*.py
├── 9.8_caracterizacion_funcional_aro.py
├── 9.9_coocurrencia_12_aro.py
└── 10.*_figura*.py
```

Los scripts `10.*` corresponden a la generación de las figuras finales.

## Estructura del repositorio

```text
TFM_resistoma/
├── README.md
├── environment.yml
├── scripts/
└── definitivo/
    ├── control_calidad/
    ├── figuras/
    ├── kma/
    ├── matrices/
    └── resultados/
```

- `scripts/`: scripts utilizados en el análisis.
- `definitivo/control_calidad/`: resultados de FastQC y fastp.
- `definitivo/kma/`: resultados principales de KMA.
- `definitivo/matrices/`: matrices del resistoma.
- `definitivo/resultados/`: resultados estadísticos y tablas.
- `definitivo/figuras/`: figuras finales.

## Reproducibilidad

El análisis se realizó en un entorno Conda denominado `tfm_resistoma`.

Principales herramientas:

| Herramienta | Versión |
|---|---:|
| Python | 3.11.15 |
| pandas | 2.3.2 |
| NumPy | 2.0.1 |
| SciPy | 1.16.0 |
| Matplotlib | 3.10.5 |
| Statsmodels | 0.15.0 |
| FastQC | 0.12.1 |
| fastp | 1.3.6 |
| KMA | 1.6.13 |
| Bowtie2 | 2.5.5 |

Para crear el entorno:

```bash
conda env create -f environment.yml
conda activate tfm_resistoma
```

## Consideraciones de interpretación

La detección de un determinante de resistencia mediante metagenómica shotgun
indica su presencia en la comunidad microbiana analizada, pero no permite
establecer por sí sola que pertenezca a *Klebsiella pneumoniae*.

La coocurrencia de determinantes tampoco constituye evidencia directa de
transferencia horizontal o de localización conjunta en un mismo
microorganismo o elemento genético.

## Trabajo Fin de Máster

**Máster en Bioinformática**  
**Universidad Internacional de Valencia (VIU)**


**Autora:** Claudia Marcela Gamboa Hernández
