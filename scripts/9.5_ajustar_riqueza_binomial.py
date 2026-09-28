#!/usr/bin/env python3

import pandas as pd
import numpy as np
import statsmodels.formula.api as smf


BASE = "/home/cgamboah/TFM_resistoma"

for umbral in [1, 5, 10]:

    archivo_riqueza = (
        f"{BASE}/definitivo/matrices/profundidad_{umbral}"
        f"/riqueza_por_muestra.tsv"
    )

    archivo_profundidad = (
        f"{BASE}/definitivo/matrices/profundidad_secuenciacion.tsv"
    )

    riqueza = pd.read_csv(archivo_riqueza, sep="\t")
    profundidad = pd.read_csv(archivo_profundidad, sep="\t")

    datos = riqueza.merge(
        profundidad[["run", "pares_despues"]],
        on="run",
        how="inner"
    )

    datos["grupo"] = pd.Categorical(
        datos["grupo"],
        categories=["No colonizado", "Colonizado"]
    )

    datos["log_profundidad"] = np.log10(
        datos["pares_despues"]
    )

    modelo = smf.negativebinomial(
        "riqueza ~ C(grupo) + log_profundidad",
        data=datos
    ).fit(disp=False)

    nombre = "C(grupo)[T.Colonizado]"

    coef = modelo.params[nombre]
    p = modelo.pvalues[nombre]

    irr = np.exp(coef)

    ci_inf = np.exp(
        modelo.conf_int().loc[nombre, 0]
    )

    ci_sup = np.exp(
        modelo.conf_int().loc[nombre, 1]
    )

    print("\n" + "=" * 60)
    print(f"DEPTH >= {umbral}")
    print("=" * 60)

    print(f"Muestras analizadas: {len(datos)}")
    print(f"IRR grupo Colonizado: {irr:.3f}")
    print(f"IC95%: {ci_inf:.3f} - {ci_sup:.3f}")
    print(f"p: {p:.6g}")

    print("\nModelo binomial negativo:")
    print(modelo.summary())
