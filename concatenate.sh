#!/bin/bash

# Define el nombre del archivo de salida
output_file="cut_sky.txt"

# Encuentra y concatena todos los archivos que empiezan con "cut_sky" y terminan en ".txt" dentro de subdirectorios
find ./z* -type f -name "cut_sky*txt" -exec cat {} + > "$output_file"

echo "Archivos concatenados en $output_file"

