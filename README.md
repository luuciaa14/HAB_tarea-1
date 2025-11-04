
# Tarea 1: Análisis Funcional 

## ¿Qué es un análisis funcional? 

El análisis funcional de genes es una etapa fundamental en bioinformática y biología molecular que busca interpretar el papel biológico de uno o varios genes dentro de la célula. 

Su objetivo principal es pasar de una lista de genes a una comprensión de los procesos biológicos, rutas metabólicas y funciones moleculares en los que esos genes participan.

En definitiva, el análisis funcional transforma datos genómicos en conocimiento biológico interpretable, ayudando a contextualizar los genes dentro de las vías metabólicas, procesos fisiológicos y mecanismos moleculares del organismo.

## Contexto biológico

En esta tarea se analizarán tres genes mitocondriales y nucleares que están estrechamente relacionados con la cadena respiratoria y con la producción de energía celular:
 * `COX4I2` Codifica una subunidad de la citocromo c oxidasa (complejo IV), la última enzima de la cadena de transporte de electrones. Participa en la reducción del oxígeno a agua y regula la eficiencia respiratoria en función de la disponibilidad de oxígeno.
 * `ND1`Forma parte del complejo I (NADH deshidrogenasa) del sistema mitocondrial. Es esencial para la transferencia de electrones desde el NADH al ubiquinona, paso inicial en la fosforilación oxidativa.
 * `ATP6`Codifica una subunidad del complejo V (ATP sintasa), responsable de la síntesis de ATP aprovechando el gradiente de protones generado por los complejos anteriores. 

 Estos tres genes desempeñan funciones clave en la fosforilación oxidativa y la producción de energía mitocondrial.

 ## Preprocesado y normalización de genes

 En primer lugar, hemos hecho un script que lee y limpia la lista de genes de entrada, aplicando una normalización heurística para los genes mitocondriales más comunes.

Para ello se han llevado a cabo los siguientes pasos:

1. **Lectura del archivo de entrada.**
    * Se ignoran las líneas vacías y comentarios (líneas que empiezan por #)
2. **Normalización heurística.**
    * Los genes mitocondriales se renombran con su forma estándar en humano (p. ej. ND1 -> MT-ND1, ATP6 -> MT-ATP6), según el diccionario `MITO_SYNONYMS` definido en el script.
3. **Exportación de resultados.**
    * Se genera un archivo .tsv con dos columnas; `gene_input` (el nombre original del gen) y `gene_normalized`(el símbolo normalizado).
    * El archivo se guarda automáticamente en la carpeta `results/`.

**Resultado:**

Después de ejecutar el comando:
```
python scripts/mi_script.py -i data/genes_input.txt -o results/analisis_cox_nd1_atp6
```
Se crea el archivo:
```
results/analisis_cox_nd1_atp6_genes_clean.tsv
```
Ejemplo de contenido:

| gene_input | gene_normalized |
|------------|----------|
COX4I2 | COX4I2
ND1 | MT-ND1
ATP6 | MT-ATP6

## Mapeo de genes con MyGene.info

Para poder hacer un análisis funcional es muy útil trabajar con identificadores estables y no solo con el símbolo que escribe el usuario. Algunos símbolos pueden tener sinónimos, mayúsculas/minúsculas o formas "no estándar" (como los genes mitocondriales).

En esta segunda etapa del script, después de normalizar los genes, se consulta el servicio público MyGene.info usando la librería `mygene`y obtiene, para cada gen:

* Símbolo oficial (`symbol`)
* Nombre del gen (`name`)
* Identificador Entrez (`entrezgene`)
* Identificador Ensembl (`ensembl_gene`)
* Si el gen no se pudo mapear

El resultado se guarda en un archivo TSV llamado, por ejemplo:

```
results/analisis_cox_nd1_atp6_mapping.tsv
```

Este archivo será la base para las etapas de enriquecimiento funcional.

## Enriquecimiento funcional (Enrichr)

Una vez los genes están normalizados y mapeados queremos saber en qué procesos aparecen estos genes más de lo esperado. Para averiguarlo, utilizamos un análisis de enriquecimiento funcional. Esta es una técnica bioinformática que identifica qué procesos biológicos o funciones están sobrerrepresentados en una lista de genes, proteínas u otras moléculas.

En este caso usamos **Enrichr** a través de la librería de Python `gseapy`. Enrichr es un servicio que tiene colecciones de genes ya agrupadas en categorías biológicas. Nosotros consultaremos:

* `GO_Biological_Process_2023`
* `GO_Molecular_Function_2023`
* `GO_Cellular_Component_2023`
* `KEGG_2021_Human`
* `Reactome_2022`

Para cada una de estas colecciones, el script:

1. Envía la lista de genes.
2. Recibe las categorías enriquecidas.
3. Las ordena por significación.
4. Guarda un archivo .tsv por colección en la carpeta `results`.

## Enriquecimiento con g:Profiler

A pesar de que Enrichr suele ser suficiente, las listas de genes suelen ser cortas o alguna de las colecciones no están disponibles. Para estos casos se usa la herramienta **g:Profiler**, usando la librería `gprofiler-official`. Esta herramienta permite anotar la lista de genes frente a varias bases de datos y devuelve una tabla con los términos significativos ordenados por el valor p. 

En esta etapa, el script:

1. Usa la lista de símbolos ya normalizados.
2. Llama a g:Profiler para el organismo humano.
3. Guarda el resultado en:

```
results/analisis_cox_nd1_atp6_gprofiler.tsv
```

## Resumen integrador de resultados

Para facilitar la interpretación, el script combina los resultados más relevantes de Enrichr y g:Profiler en un único archivo resumen:

```
results/analisis_cox_nd1_atp6_resumen.tsv
```

Este resumen incluye los términos más significativos de cada fuente, ordenados por su valor p o valor p ajustado, permitiendo comparar rápidamente las categorías comunes.

**Visualización del resumen:**

| source	| term	| p_value |
|-----------|-------|---------|
| g:Profiler	| respiratory chain complex |	1.7495192096134485e-05 |					
| g:Profiler 	| oxidative phosphorylation	| 0.00016400056416194064 |					
| g:Profiler	| Oxidative phosphorylation	| 0.00019077901673278882 |				
| g:Profiler	| proton transmembrane transporter activity	| 0.00025605478100053635 |				
| g:Profiler	| proton transmembrane transport	| 0.000363191140733216 |				
| g:Profiler	| aerobic respiration	| 0.0004175865333328184 |			
| g:Profiler	| Diabetic cardiomyopathy	| 0.0006585331746415463 |	
| g:Profiler	| cellular respiration	| 0.0007530184933910186 |
| g:Profiler	| Chemical carcinogenesis - reactive oxygen species	| 0.0008634939616625755 |			
| g:Profiler	| Thermogenesis	| 0.0009996041175237241 |	


## Conclusiones y referencias

El análisis confirma que los genes COX4I2, ND1 y ATP6 están asociados a procesos de fosforilación oxidativa y respiración mitocondrial, concordando con su función conocida en la cadena transportadora de electrones.

**Bases de datos y recursos utilizados:**

*https://mygene.info/*
*https://maayanlab.cloud/Enrichr/*
*https://biit.cs.ut.ee/gprofiler/*
*https://geneontology.org/*
*https://www.genome.jp/kegg/*
*https://reactome.org/*