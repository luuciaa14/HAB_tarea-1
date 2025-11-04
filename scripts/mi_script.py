import argparse
import os
import sys
from typing import List
import pandas as pd
import mygene
from gseapy import enrichr
from gprofiler import GProfiler

#Librerias de Enrichr

ENRICHR_LIBRARIES = [
    "GO_Biological_Process_2023",
    "GO_Molecular_Function_2023",
    "GO_Cellular_Component_2023",
    "KEGG_2021_Human",
    "Reactome_2022",
]

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
    fixed = []
    for g in genes:
        fixed.append(MITO_SYNONYMS.get(g.upper(), g))
    return fixed

# Función de mapeo
def map_genes_with_mygene(genes: List[str]) -> pd.DataFrame:
    mg = mygene.MyGeneInfo()
    try:
        df = mg.querymany(
            genes,
            scopes = "symbol,alias,names",
            fields = "symbols,name,entrezgene,ensembl.gene,taxid",
            species = "human",
            as_dataframe = True,
            returnall = False,
            verbose = False,
        )
    except Exception as e:
        print(f"[ERROR] Fallo consultando MyGene.info: {e}", file = sys.stderr)
        return pd.DataFrame()
    
    # Si viene vacío
    if df is None or df.empty:
        return pd.DataFrame()
    
    # Pasamos a un mejor formato
    df = df.reset_index().rename(columns={"query": "query_in"})

    def _extract_ensembl(val):
        if isinstance(val, dict) and "gene" in val:
            return val["gene"]
        if isinstance(val, list) and len(val) > 0:
            first = val[0]
            if isinstance(first, dict) and "gene" in first:
                return first["gene"]
        return None
    
    if "ensembl" in df.columns:
        df["ensembl_gene"] = df["ensembl"].apply(_extract_ensembl)
    else:
        df["ensembl_gene"] = None

    # Columnas clave
    wanted = ["query_in", "symbol", "name", "entrezgene", "ensembl_gene", "taxid", "notfound"]
    for col in wanted:
        if col not in df.columns:
            df[col] = None
    df = df[wanted].drop_duplicates()
    return df

def run_enrichr(gene_symbols: List[str], out_prefix: str) -> dict:
    print("[INFO] Ejecutando enriquecimiento con Enrichr")
    results = {}
    for lib in ENRICHR_LIBRARIES:
        try:
            enr = enrichr(
                gene_list=gene_symbols,
                gene_sets=[lib],
                organism="Human",
                cutoff=1.0,
            )

            if enr.results is None or enr.results.empty:
                print("[WARN] Enrichr no devolvió resultados para {lib}")
                continue

            df_res = enr.results.copy()

            rename_map = {
                "Term": "term",
                "Adjusted P-value": "adj_p",
                "P-value": "p_value",
                "Odds Ratio": "odds_ratio",
                "Combined Score": "combined_score",
                "Genes": "overlap_genes",
                "Overlap": "overlap",
            }

            for old, new in rename_map.items():
                if old in df_res.columns:
                    df_res.rename(columns={old: new}, inplace=True)

            if "adj_p" in df_res.columns:
                df_res = df_res.sort_values("adj_p", ascending=True)
            elif "p_value" in df_res.columns:
                df_res = df_res.sort_values("p_value", ascending=True)

            out_path = f"{out_prefix}_enrichr_{lib}.tsv"
            df_res.to_csv(out_path, sep="\t", index=False)
            print(f"[OK] Guardado enriquecimiento de {lib} en: {out_path}")

            results[lib] = df_res

        except Exception as e:
            print(f"[WARN] Fallo al consultar {lib} en Enrichr: {e}", file=sys.stderr)
            continue

def run_gprofiler(gene_symbols: List[str], out_prefix: str) -> pd.DataFrame:
    print("[INFO] Ejecutando enriquecimiento con g:Profiler...")
    try:
        gp = GProfiler(return_dataframe=True)
        res = gp.profile(
            organism="hsapiens",
            query=gene_symbols,
            user_threshold=1.0,
            all_results=True,
        )
    except Exception as e:
        print(f"[WARN] g:Profiler falló: {e}", file=sys.stderr)
        return pd.DataFrame()
    
    if res is None or res.empty:
        print("[WARN] g:Profiler no devolvió resultados.")
        return pd.DataFrame()

    # Ordenamos por p-valor si está
    if "p_value" in res.columns:
        res = res.sort_values("p_value", ascending=True)

    out_path = f"{out_prefix}_gprofiler.tsv"
    res.to_csv(out_path, sep="\t", index=False)
    print(f"[OK] Guardado enriquecimiento g:Profiler en: {out_path}")
    return res

def make_summary(enrichr_results: dict, gprof_df: pd.DataFrame, out_prefix: str) -> None:
    frames = []

    if enrichr_results is None:
        enrichr_results = {}

    # De Enrichr cogemos top 10 de cada librería
    for lib, df in enrichr_results.items():
        if df is None or df.empty:
            continue

        top = df.head(10).copy()
        top["source"] = f"Enrichr:{lib}"

        for col in ["term", "adj_p", "p_value", "combined_score", "odds_ratio", "overlap", "overlap_genes"]:
            if col not in top.columns:
                top[col] = None

        frames.append(top[["source", "term", "adj_p", "p_value", "combined_score", "odds_ratio", "overlap", "overlap_genes"]])

    # De g:Profiler cogemos top 10 por p_value
    if gprof_df is not None and not gprof_df.empty:
        gtmp = gprof_df.copy()
        gtmp["term"] = gtmp.get("term_name", gtmp.get("name", None))
        gtmp["source"] = "g:Profiler"

        for col in ["adj_p", "combined_score", "odds_ratio", "overlap", "overlap_genes"]:
            if col not in gtmp.columns:
                gtmp[col] = None

        gtmp = gtmp.sort_values("p_value", ascending=True).head(10)
        frames.append(gtmp[["source", "term", "p_value", "adj_p", "combined_score", "odds_ratio", "overlap", "overlap_genes"]])

    if not frames:
        print("[WARN] No hay resultados para resumir.")
        return

    summary_df = pd.concat(frames, ignore_index=True)
    out_path = f"{out_prefix}_resumen.tsv"
    summary_df.to_csv(out_path, sep="\t", index=False)
    print(f"[OK] Guardado resumen integrado en: {out_path}")

def main():
    # Definición de argumentos
    # Archivo de entrada y prefijo de salida
    parser = argparse.ArgumentParser(
        description= "Preprocesado de genes: lectura + normalización mitocondrial + mapeo con MyGene.info."
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
        print(f"[ERROR] No se encuentra el archivo de entrada: {input_path}", file=sys.stderr)
        sys.exit(1)

    # Archivo vacío o sin genes válidos
    if not genes_raw:
        print("[ERROR] El archivo de entrada no contiene genes válidos.", file=sys.stderr)
        sys.exit(1)

    # 2. Normalización heurística de genes mitocondriales
    genes_fixed = heuristics_fix_mito_symbols(genes_raw)

    # 3. Guardar resultados en un .tsv    
    clean_path = f"{out_prefix}_genes_clean.tsv"
    df_clean = pd.DataFrame({"gene_input": genes_raw, "gene_normalized": genes_fixed})
    df_clean.to_csv(clean_path, sep="\t", index=False)
    print(f"[OK] Guardado preprocesado en: {clean_path}")

    # 4. Mapeo con MyGene.info
    mapping_df = map_genes_with_mygene(genes_fixed)

    if mapping_df is None or mapping_df.empty:
        print("[WARN] No se obtuvieron resultados de mapeo. Revisa los nombres de los genes o la conexión.")
    else:
        mapping_path = f"{out_prefix}_mapping.tsv"
        mapping_df.to_csv(mapping_path, sep="\t", index=False)
        print(f"[OK] Guardado mapeo en: {mapping_path}")

    # 5. Enriquecimiento (usando los símbolos normalizados)
    enrichr_results = run_enrichr(genes_fixed, out_prefix)

    # 6. Enriquecimiento con g:Profiler
    gprof_df = run_gprofiler(genes_fixed, out_prefix)

    # 7. Resumen
    make_summary(enrichr_results, gprof_df, out_prefix)

    print("[OK] Versión final completada.")

if __name__ == "__main__":
    main()