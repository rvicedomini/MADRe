#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# MADRe GTDB Reference Builder
# Builds a combined FASTA: >taxid|accession
#
# Usage:
#   bash build_madre_gtdb.sh \
#       --tar gtdb_genomes_reps_r226.tar.gz \
#       --meta bac120_metadata_r226.tsv \
#       --out MADRe_reference \
#       --threads 32
###############################################################################

# ------------------------ Parse arguments ------------------------------------
TAR=""
META=""
OUTDIR="MADRe_reference"
THREADS=$(command -v nproc >/dev/null 2>&1 && nproc || echo 8)

while [[ $# -gt 0 ]]; do
  case "$1" in
    -t|--tar) TAR="$2"; shift 2 ;;
    -m|--meta) META="$2"; shift 2 ;;
    -o|--out) OUTDIR="$2"; shift 2 ;;
    -p|--threads) THREADS="$2"; shift 2 ;;
    -h|--help)
      echo "Usage: $0 --tar <gtdb_tar.gz> --meta <metadata.tsv> [--out <dir>] [--threads <n>]"
      exit 0
      ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

if [[ -z "$TAR" || -z "$META" ]]; then
  echo "Error: must provide both --tar and --meta" >&2
  exit 1
fi

# ------------------------- Setup ---------------------------------------------
OUT_FASTA="${OUTDIR}/$(basename "${TAR%.tar.gz}")_MADRe.fna"
SCRATCH="${OUTDIR}/_tmp"
# mkdir -p "$OUTDIR" "$SCRATCH"

# Ensure output and scratch directories exist
if [[ ! -d "$OUTDIR" ]]; then
  echo "Creating output directory: $OUTDIR"
  mkdir -p "$OUTDIR"
fi

if [[ ! -d "$SCRATCH" ]]; then
  echo "Creating temporary scratch directory: $SCRATCH"
  mkdir -p "$SCRATCH"
fi

echo "Tarball: $TAR"
echo "Metadata: $META"
echo "Output directory: $OUTDIR"
echo "Threads: $THREADS"
echo

# ------------------------- Check tools ---------------------------------------
for tool in tar awk find zcat; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "Error: '$tool' not found in PATH." >&2
    exit 1
  fi
done
PARALLEL=$(command -v parallel || echo "")
PIGZ=$(command -v pigz || echo "")
ZCAT=$(command -v zcat || command -v gzcat || echo "gzip -dc")

# ------------------------- Step 1: Extract tarball ---------------------------
echo "[1/6] Extracting genomes from tarball..."
if [[ -n "$PIGZ" ]]; then
  tar --use-compress-program="$PIGZ" -xf "$TAR" -C "$SCRATCH"
else
  tar -xzf "$TAR" -C "$SCRATCH"
fi

# ------------------------- Step 2: List all FASTAs ---------------------------
echo "[2/6] Finding .fna.gz files..."
find "$SCRATCH" -type f -name '*_genomic.fna.gz' > "$SCRATCH/files.list"
echo "  $(wc -l < "$SCRATCH/files.list") files found."

# ------------------------- Step 3: Build acc→taxid map -----------------------
echo "[3/6] Building accession to taxid map..."
NCBI_TAX_COL=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="ncbi_taxid"){print i; exit}}' "$META")
if [[ -z "$NCBI_TAX_COL" ]]; then
  echo "Error: column 'ncbi_taxid' not found in $META" >&2
  exit 1
fi
awk -F'\t' -vC="$NCBI_TAX_COL" 'NR>1 && $1!="" && $C!="" {print $1"\t"$C}' "$META" > "$SCRATCH/acc2taxid.tsv"

# ------------------------- Step 4: Make job table ----------------------------
echo "[4/6] Preparing job list..."
awk -vMAP="$SCRATCH/acc2taxid.tsv" -vOFS="\t" '
  BEGIN{
    while((getline < MAP)>0){if($1!=""&&$2!=""){m[$1]=$2}}
    close(MAP)
  }
  {
    f=$0
    if (match(f,/(GCA|GCF)_[0-9]+\.[0-9]+/,a)){
      acc=a[0]
      tax=m[acc]
      if (tax=="") tax=m["RS_"acc]
      if (tax=="") tax=m["GB_"acc]
      if (tax=="") tax="0"
      print f, acc, tax
    }
  }
' "$SCRATCH/files.list" > "$SCRATCH/jobs.tsv"

echo "  $(wc -l < "$SCRATCH/jobs.tsv") genomes mapped to taxids."

# ------------------------- Step 5: Rewrite headers ---------------------------
echo "[5/6] Rewriting headers in parallel..."
mkdir -p "$SCRATCH/chunks"

if [[ -n "$PARALLEL" ]]; then
  parallel --colsep '\t' -j "$THREADS" --halt soon,fail=1 '
    f={1}; acc={2}; tax={3};
    '"$ZCAT"' "$f" | awk -vA="$acc" -vT="$tax" '"'"'
      /^>/ {print ">|"T"|"A; next}
           {print}
    '"'"' > "'"$SCRATCH"'/chunks/{#}.fna"
  ' :::: "$SCRATCH/jobs.tsv"
else
  echo "GNU parallel not found; running sequentially."
  i=0
  while IFS=$'\t' read -r f acc tax; do
    i=$((i+1))
    $ZCAT "$f" | awk -vA="$acc" -vT="$tax" '/^>/{print ">|"T"|"A; next}{print}' \
      > "$SCRATCH/chunks/${i}.fna"
  done < "$SCRATCH/jobs.tsv"
fi

# ------------------------- Step 6: Merge & compress --------------------------
echo "[6/6] Merging and compressing..."

find "${SCRATCH}" -type f -name "*.fna" -print0 | xargs -0 cat | pigz -p ${THREADS} > "${OUTDIR}/madre_gtdb_r226.fna.gz"

# ------------------------- Cleanup ------------------------------------------
echo "Cleaning up temporary files..."
rm -rf "$SCRATCH"

echo
echo "Done. Output: ${OUTDIR}/madre_gtdb.fna.gz"
