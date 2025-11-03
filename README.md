
# Tarea 1: Análisis Funcional 

## ¿Qué es un análisis funcional? 

El análisis funcional de genes es una etapa fundamental en bioinformática y biología molecular que busca interpretar el papel biológico de uno o varios genes dentro de la célula. 

Su objetivo principañ es pasar de una lista de genes a una comprensión de los procesos biológicos, rutas metabólicas y funciones moleculares en los que esos genes participan.

En definitiva, el análisis funcional transforma datos genómicos en conocimiento biológico interpretables, ayudando a contextualizar los genes dentro de las vías metabólicas, procesos fisiológicos y mecanismos molesculares del organismo.

## Contexto biológico

En esta tarea se analizarán tres genes mitocondriales y nucleares que están estrechamente relacionados con la cadena respiratoria y con la producción de energía celular:
 * `COX4I2` Codifica una subunidad de la citocromo c oxidasa (complejo IV), la última enzima de la cadena de transporte de electrones. Participa en la reducción del oxígeno a agua y regula la eficiencia respiratoria en función de la disponibilidad de oxígeno.
 * `ND1`Forma parte del complejo I (NADH deshidrogenasa) del sistema mitocondrial. Es esencial para la transferencia de electrones desde el NADH al ubiquinona, paso inicial en la fosforilación oxidativa.
 * `ATP6`Codifica una subunidad del complejo V (ATP sintasa), responsable de la síntesis de ATP aprovechando el gradiente de protones generado por los complejos anteriores. 

 Estos tres genes desempeñan funciones clave en la fosforilación oxidativa y la producción de energía mitocondrial.

 ## Preprocesado y normalización de genes

 En un primer lugar, hemos hecho un script que lee y limpia la lista de genes de entrada, aplicando una normalización heurística para los genes mitocondriales más comunes.

Para ello se han llevado a cabo los siguientes pasos:

1. **Lectura del archivo de entrada.**
    * Se ignoran las líneas vacías y comentarios (líneas que empeizan por #)
2. **Normalización heurística.**
    * Los genes mitocondriales se renombran con su forma entándar en humano (p. ej. ND1 -> MT-ND1, ATP6 -> MT-ATP6), según el diccionario `MITO_SYNONYMS` definido en el script.
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
Ejemplo contenido:

| gene_input | gene_normalized |
|------------|----------|
COX4I2 | COX4I2
ND1 | MT-ND1
ATP6 | MT-ATP6