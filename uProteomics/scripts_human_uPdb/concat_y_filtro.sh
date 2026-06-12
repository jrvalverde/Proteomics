 #!/bin/bash

# 1. Comprobamos que el usuario ha pasado un argumento
if [ -z "$1" ]; then
    echo "Uso: ./concat_y_filtro
    .sh <directorio>"
    exit 1
fi

# 2. Comprobamos que existen archivos .fa en el directorio
FILES=("$INPUT_DIR"/*_nr.fa)

if [ ${#FILES[@]} -eq 0 ]; then
    echo "No hay archivos *_nr.fa en el directorio"
    exit 1
fi

# 3. Guardamos el directorio en una variable
INPUT_DIR="$1"

# 4. Creamos variable para el directorio output
OUTPUT_DIR="$INPUT_DIR/DB_concat_final"

# 5. Creamos el directorio output
mkdir -p "$OUTPUT_DIR"

# 6 Concatenar archivos del directorio
cat "$INPUT_DIR"/*_nr.fa > "$OUTPUT_DIR/DB_concat_final.fa"
 
# 7. Mensaje de aviso de unión archivos del directorio
echo "Concatenando todos los archivos _nr.fa del directorio ${INPUT_DIR}"
    
# 8. Contar secuencias archivo concatenado
SEQ_INICIO=$(grep -c '^>' $OUTPUT_DIR/DB_concat_final.fa)
    
# 9. CD-HIT
cd-hit -i "$OUTPUT_DIR/DB_concat_final.fa" \
-o "$OUTPUT_DIR/DB_concat_final_nr2.fa" \
-c 1.0 -n 5 -d 0 -M 3000 -T 25

# 10. Contar secuencias archivo output de cd-hit
SEQ_FINAL=$(grep -c '^>' $OUTPUT_DIR/DB_concat_final_nr2.fa)

# 11. Mensaje resumen de resultados del filtrado
echo "Procesado terminado"
echo "  Secuencias original: ${SEQ_INICIO}"
echo "  Secuencias final: ${SEQ_FINAL}"

# 12. Guardamos nombre y ruta del archivo resumen en una variable
RESUMEN="$OUTPUT_DIR/resumen_DB_concat_final.txt"

# 13. Guardamos info del mensaje resumen en archivo txt
echo -e "DB_concat_final_nr2.fa\nSecuencias archivo concatenado: ${SEQ_INICIO}\nSecuencias archivo filtrado (cd-hit): ${SEQ_FINAL}" >> "$RESUMEN"  
echo "Resumen resultados guardado en ${RESUMEN}"
    
    
    
    
