#!/bin/bash
########## This code counts the E-MTAB-11188 FASTQs of each run with kallisto | bustools (exon-only) ##########
# kb "standard" workflow: only reads compatible with mature (spliced) transcripts are counted, as in the Cell Ranger v2
# matrices of GSE180298. Output: one unfiltered, Cell Ranger-style count matrix per run (E-MTAB-11188/kb_counts/RUN/
# cellranger), read by 3.3.3. The FASTQ MD5 checksums are verified against the ENA file report written by 3.3.1.
# Software (README): kb-python 0.30.2 (env/kb-python_requirements.txt), kallisto 0.52.0, bustools 0.45.1 and the human
# "standard" index of pachterlab/kallisto-transcriptome-indices v1 (index.idx, t2g.txt).
# Usage : bash 3.3.2.scRNAseq_EMTAB11188_kb_count.sh RUN [RUN ...]   (KB, KALLISTO, BUSTOOLS and INDEX_DIR set the tool and
#         index locations; EMTAB11188_DIR or DATA_DIR change the data folder, as in 0.config.R)

set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")" && pwd)"
DATASET_DIR="${EMTAB11188_DIR:-${DATA_DIR:-$REPO_DIR/data}/E-MTAB-11188}"
FASTQ_DIR="${FASTQ_DIR:-$DATASET_DIR/fastq}"
MANIFEST="$DATASET_DIR/ENA_manifest.tsv"
WORK_DIR="${WORK_DIR:-$DATASET_DIR/kb_work}"  # large intermediates (BUS files)
OUT_DIR="$DATASET_DIR/kb_counts"
INDEX_DIR="${INDEX_DIR:-${DATA_DIR:-$REPO_DIR/data}/kb_human_standard}"
KB="${KB:-kb}"
THREADS="${THREADS:-8}"

# kallisto and bustools binaries passed to kb only when set (on macOS the kallisto bundled with kb-python 0.30.2 needs an
# HDF5 library that Homebrew no longer ships; the official kallisto 0.52.0 release binary is used instead)
KB_BINARIES=()
[ -n "${KALLISTO:-}" ] && KB_BINARIES+=(--kallisto "$KALLISTO")
[ -n "${BUSTOOLS:-}" ] && KB_BINARIES+=(--bustools "$BUSTOOLS")

md5_of() { if command -v md5 > /dev/null; then md5 -q "$1"; else md5sum "$1" | cut -d ' ' -f 1; fi; }
expected_md5() { awk -F'\t' -v f="$1" 'NR>1{n=split($5,u,";"); split($6,m,";"); for(i=1;i<=n;i++){k=split(u[i],p,"/"); if(p[k]==f) print m[i]}}' "$MANIFEST"; }

[ $# -ge 1 ] || { echo "Usage: $0 RUN [RUN ...]"; exit 1; }

for RUN in "$@"; do
  if [ -f "$OUT_DIR/$RUN/run_info.json" ]; then
    echo "[$RUN] already counted, skipping"
    continue
  fi

  ##### FASTQs of the run and their MD5 checksums
  R1="$FASTQ_DIR/${RUN}_1.fastq.gz"
  R2="$FASTQ_DIR/${RUN}_2.fastq.gz"
  for f in "$R1" "$R2"; do
    if [ "$(md5_of "$f")" != "$(expected_md5 "$(basename "$f")")" ]; then
      echo "[$RUN] ERROR: MD5 mismatch for $f"
      exit 1
    fi
  done

  ##### Pseudoalignment, barcode correction (10x v3 on-list) and counting, without cell filtering
  mkdir -p "$WORK_DIR/$RUN"
  "$KB" count -i "$INDEX_DIR/index.idx" -g "$INDEX_DIR/t2g.txt" -x 10xv3 -t "$THREADS" --cellranger \
    ${KB_BINARIES[@]+"${KB_BINARIES[@]}"} -o "$WORK_DIR/$RUN" "$R1" "$R2"

  ##### Count matrix in the Cell Ranger v3 layout expected by Read10X (genes without a symbol keep their Ensembl ID)
  mkdir -p "$OUT_DIR/$RUN/cellranger"
  cp "$WORK_DIR/$RUN/counts_unfiltered/cellranger/matrix.mtx.gz" "$WORK_DIR/$RUN/counts_unfiltered/cellranger/barcodes.tsv.gz" "$OUT_DIR/$RUN/cellranger/"
  gzip -dc "$WORK_DIR/$RUN/counts_unfiltered/cellranger/genes.tsv.gz" \
    | awk -F'\t' 'BEGIN{OFS="\t"} {name=($2=="" ? $1 : $2); print $1, name, "Gene Expression"}' \
    | gzip > "$OUT_DIR/$RUN/cellranger/features.tsv.gz"

  ##### Run statistics (reads processed, % pseudoaligned, barcodes, kb call); run_info.json last, as it marks the run as done
  cp "$WORK_DIR/$RUN/inspect.json" "$WORK_DIR/$RUN/kb_info.json" "$OUT_DIR/$RUN/"
  cp "$WORK_DIR/$RUN/run_info.json" "$OUT_DIR/$RUN/"
  echo "[$RUN] done"
done
