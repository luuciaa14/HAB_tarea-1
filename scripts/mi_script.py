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
            parts = [g.strip() for g in s.replace(",", " ").split() if g.strip()]
            genes.extend(parts)
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
    parser.add_argument("-o", "--output", required=True, help="Prefijo de salida, p. ej. results/analisis_cox_nd1_atp6.")
    
    # Procesamiento de argumentos
    args = parser.parse_args()

    input_path = os.path.abspath(args.input)
    out_prefix = os.path.abspath(args.output)
    out_dir = os.path.dirname(out_prefix) if os.path.dirname(out_prefix) else os.path.abspath("results")

    os.makedirs(out_dir, exist_ok=True)
    
    print(f"[INFO] CWD: {os.getcwd()}")
    print(f"[INFO] Input : {input_path}")
    print(f"[INFO] Output prefix: {out_prefix}")
    print(f"[INFO] Output dir   : {out_dir}")

    # 1. Leer la lista de genes del archivo de entrada
    try:
        genes_raw = read_gene_list(input_path)
    except FileNotFoundError:
        print(f"[ERROR] No se encuentra el archivo de entrada: {args.input}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"[ERROR] Fallo leyendo el archivo de entrada: {e}", file=sys.stderr)
        sys.exit(1)

    # Archivo vacío o sin genes válidos
    if not genes_raw:
        print("[ERROR] El archivo de entrada no contiene genes válidos.", file=sys.stderr)
        sys.exit(1)

    # 2. Normalización heurística de genes mitocondriales
    genes_fixed = heuristics_fix_mito_symbols(genes_raw)

    # 3. Guardar resultados en un .tsv    
    out_clean = f"{out_prefix}_genes_clean.tsv"
    try:
        df = pd.DataFrame({"gene_input": genes_raw, "gene_normalized": genes_fixed})
        df.to_csv(out_clean, sep="\t", index=False)
    except Exception as e:
        print(f"[ERROR] No pude escribir el TSV de salida: {e}", file=sys.stderr)
        sys.exit(1)

    # 4. Mensaje final informativo
    print("Preprocesamiento completado.")
    print(f"- Entrada: {args.input}")
    print(f"- Salida (genes normalizados): {out_clean}")

if __name__ == "__main__":
    main()