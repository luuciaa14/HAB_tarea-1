import argparse
import os
import sys
from typing import List

import pandas as pd

# Sinónimos de genes mitocondriales habituales:
MITO_SYNONYMS = {
    "ND1": "MT-ND1",
    "ND2": "MT-ND2",
    "ND3": "MT-ND3",
    "ND4": "MT-ND4",
    "ND4L": "MT-ND4L",
    "ND5": "MT-ND5",
    "ND6": "MT-ND6",
    "CYB": "MT-CYB",
    "CO1": "MT-CO1",
    "CO2": "MT-CO2",
    "CO3": "MT-CO3",
    "ATP6": "MT-ATP6",
    "ATP8": "MT-ATP8",
    "RNR1": "MT-RNR1",
    "RNR2": "MT-RNR2",
}

# FUNCIONES

# Leer el archivo de entrada, ignorando líneas que empiezan
# por # y las que están en blanco.
def read_gene_list(path: str) -> List[str]:
    genes: List[str] = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            genes.append(s)
    return genes

# Normalización heurística para genes habituales
def heuristics_fix_mito_symbols(genes: List[str]) -> List[str]:
    """
    Aplica normalización heurística para símbolos mitocondriales comunes.
    Si el símbolo no está en el diccionario, se deja tal cual.
    """
    fixed = []
    for g in genes:
        fixed.append(MITO_SYNONYMS.get(g.upper(), g))
    return fixed

def main():
    # Definición de argumentos
    # Archivo de entrada y prefijo de salida
    parser = argparse.ArgumentParser(
        description= "Preprocesado de genes: lectura + normalización mitocondrial."
    )
    parser.add_argument("-i", "--input", required=True, help="Ruta al archivo de genes.")
    parser.add_argument("-o", "--output", required=True, help="Prefijo de slaida, p. ej. results/analisis_cox_nd1_atp6.")
    
    # Procesamiento de argumentos
    args = parser.parse_args()

    # Preparar la ruta de salida
    out_prefix = args.output

    # Si no se indica directorio, "reuslts" por defecto
    # Si no existe directorio, se crea
    out_dir = os.path.dirname(out_prefix) if os.path.dirname(out_prefix) else "results"
    os.makedirs(out_dir, exist_ok=True)

    # 1. Leer la lista de genes del archivo de entrada
    try:
        genes_raw = read_gene_list(args.input)
    except FileNotFoundError:
        print(f"[ERROR] No se encuentra el archivo de entrada: {args.input}", file=sys.stderr)
        sys.exit(1)

    # Archivo vacío o sin genes válidos
    if not genes_raw:
        print("[ERROR] El archo de entrada no contiene genes válidos.", file=sys.stderr)
        sys.exit(1)

    # 2. Normalización heurística de genes mitocondriales
    genes_fixed = heuristics_fix_mito_symbols(genes_raw)

    # 3. Guardar resultados en un .tsv
    df = pd.DataFrame(
        {
            "gene_input": genes_raw,
            "gene_normalized": genes_fixed,
        }
    )
    
    out_clean = f"{out_prefix}_genes_clean_tsv"
    df.to_csv(out_clean, sep="\t", index=False)

    # 4. Mensaje final informativo
    print("Preprocesamiento completado.")
    print(f"- Entrada: {args.input}")
    print(f"- Salida (gens normalizados): {out_clean}")