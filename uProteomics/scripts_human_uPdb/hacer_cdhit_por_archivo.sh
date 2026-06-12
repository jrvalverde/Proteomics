#!/bin/bash

# 1. Comprobamos que el usuario ha pasado un argumento
if [ -z "$1" ]; then
    echo "Uso: ./hacer_cdhit_por_archivo.sh <directorio>"
    exit 1
fi

# 2. Guardamos el directorio en una variable
DIR="$1"

# 3. Creamos la carpeta de cd-hits dentro del directorio (si no existe)
## -p evita que de error si carpeta ya existe.
mkdir -p "$DIR/cd-hits"

# 4. Creamos un archivo resumen de número de secuencias antes y después del procesado
resumen="$DIR/cd-hits/resumen.cd-hits.txt"
echo -e "Archivo\tNum_sec_originales\Num_sec_nr_redundantes" > "$resumen"

shopt -s nullglob # si no hay archivos .fa el bucle no se ejecutará

# 5. Recorremos todos los archivos .fa dentro del directorio
for file in "$DIR"/*.fa; do
    
    # 5.1 Extraemos el nombre del archivo sin ruta con basename
    base=$(basename "$file")
    
    # 5.2 Quitamos la extensión .fa
    name="${base%.fa}"
    
    # 5.3 Guardamos en una variable el nombre del archivo de salida
    output="$DIR/cd-hits/${name}_cdhit.fa"
    
    # 5.3 Ejecutamos el cd-hit en el archivo
    cd-hit -i "$file" \
    -o "$output" \
    -c 1.0 -n 5 -d 0 -l 8 -M 3000 -T 25 > /dev/null # ocultamos texto que sale del cd-hit
    
    # 5.4 Contar número de secuencias en archivo original y filtrado
    original=$(grep -c "^>" "$file")
    filtrado=$(grep -c "^>" "$output")
    
    echo "  Original: $original"
    echo "  No redundantes: $filtrado"
    echo "--------------------------------------"
    
    # 5.5 Guardamos estos resultados en el archivo resumen creado
    echo -e "$base\t$original\t$filtrado" >> "$resumen" 
done
