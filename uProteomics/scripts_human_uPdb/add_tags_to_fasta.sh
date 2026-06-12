#!/bin/bash

indir="$1"

if [[ -z "$indir" || ! -d "$indir" ]]; then
  echo "Usage: $0 fasta-folder"
  exit 1
fi

DB=$(basename "$indir")
outdir="${indir}_tagged"

mkdir -p "$outdir"

for infile in "$indir"/*.fa*; do
  # Saltar si no es archivo
  [[ -f "$infile" ]] || continue

  fname=$(basename "$infile")

  # Quitar extensiones múltiples: .fasta, .fa, .faa, .gz, .zip
  base="$fname"
  base="${base%.*}"

  outfile="$outdir/$base.fa"

  echo "Tagging: $fname → $outfile (DB=$DB)"

  if [[ "$infile" == *.gz ]]; then
    zcat "$infile"
  elif [[ "$infile" == *.zip ]]; then
    unzip -p "$infile"
  else
    cat "$infile"
  fi \
  | sed "/^>/ s/$/ DB=$DB/" > "$outfile"
done

echo "✅ Finished, results are in: $outdir"
