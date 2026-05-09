#!/usr/bin/env bash
set -euo pipefail

label="$1"
indir="$2"
outdir="$3"

mkdir -p "${outdir}"

for f in "${indir}"/*.fna "${indir}"/*.fa "${indir}"/*.fasta "${indir}"/*.fna.gz
do
  [ -e "$f" ] || continue

  base=$(basename "$f")
  base=${base%.gz}
  base=${base%.fna}
  base=${base%.fa}
  base=${base%.fasta}

  outfile="${outdir}/${label}__${base}.fna"

  if [[ "$f" == *.gz ]]; then
    gzip -cd "$f" > "${outfile}"
  else
    ln -sf "$(realpath "$f")" "${outfile}"
  fi
done
