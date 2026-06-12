 #!/bin/bash

# 1. Comprobamos que el usuario ha pasado un argumento
if [ -z "$1" ]; then
    echo "Uso: ./hacer_DB_nr.sh <directorio>"
    exit 1
fi

# 2. Guardamos el directorio en una variable
INPUT_DIR="$1"

# 3. Creamos variable para el directorio output
OUTPUT_DIR="$INPUT_DIR/DB_nr"

# 4. Creamos el directorio output
mkdir -p "$OUTPUT_DIR"

# 5. Guardamos nombre y ruta del archivo resumen en una variable
RESUMEN="$OUTPUT_DIR/resumen_DB_nr.txt"

# 6. Creamos archivo resumen y añadimos encabezados
echo -e "Base_de_datos\tSecuencias_antes\tSecuencias_despues" > "$RESUMEN"

# 7. Obtener nombres unicos de las bases de datos (antes del primer guión) con bucle for
for DB in $(ls $INPUT_DIR/*.fa | xargs -n1 basename | cut -d'_' -f1 | sort | uniq)
do
    
    echo "Procesando base de datos: $DB"
    
    ARCHIVOS=$(ls $INPUT_DIR/${DB}_*.fa 2>/dev/null)
    NUM_ARCHIVOS=$(echo "$ARCHIVOS" | wc -w)
    
    # 7.1 Contar secuencias originales
    SEQ_ANTES=$(grep -c "^>" $INPUT_DIR/${DB}_*.fa)
    
    if  [ "$NUM_ARCHIVOS" -gt 1 ]; then
    	echo " Hay $NUM_ARCHIVOS archivos -> se concatenan y se hace CD-HIT"
	
	# 7.1.1 Concatenar archivos
	cat $INPUT_DIR/${DB}_*.fa > $OUTPUT_DIR/${DB}_concat.fa
	
	# 7.1.2 CD-HIT"
	cd-hit -i "$OUTPUT_DIR/${DB}_concat.fa" \
	-o "$OUTPUT_DIR/${DB}_nr.fa" \
	-c 1.0 -n 5 -d 0 -M 3000 -T 25
    SEQ_DESPUES=$(grep -c "^>" $OUTPUT_DIR/${DB}_nr.fa)
    
    rm $OUTPUT_DIR/${DB}_concat.fa
    
    else
        echo " Solo hay 1 archivo -> NO se hace CD-HIT"
	
	# 7.1.3 Copiar directamente
	cp $INPUT_DIR/${DB}_*.fa $OUTPUT_DIR/${DB}_nr.fa
	
	SEQ_DESPUES=$SEQ_ANTES
    fi
    
    echo -e "${DB}\t${SEQ_ANTES}\t${SEQ_DESPUES}" >> "$RESUMEN"
done
echo "Proceso terminado >_<"

	
	
	
	
