#!/bin/bash
########## This code downloads the FASTQs of the 8 analysed E-MTAB-11188 runs from ENA and verifies their MD5 checksums ##########
# E-MTAB-11188 (ENA PRJEB48995): 10x Genomics 3' v3 libraries of FACS-sorted neutrophils. Only the runs of healthy donors
# are downloaded (about 30-46 GB per run):
#   BM, steady state          : ERR7425213 (Sample 3), ERR7425214 (Sample 4)
#   PB, steady state          : ERR7425211 (Sample 1), ERR7425212 (Sample 2)
#   PB, G-CSF-treated donors  : ERR7425188 (Sample 5), ERR7425190 (Sample 6), ERR7425191 (Sample 7), ERR7425189 (Sample 8)
# The G-CSF-treated donors are processed with the others up to the QC (3.3.4) and removed in 3.3.5.
# Interrupted transfers are resumed; files that are already complete (MD5 match) are skipped. The ENA file report
# (ENA_manifest.tsv) is kept for the MD5 check in 3.3.2.
# Usage : bash 3.3.1.scRNAseq_EMTAB11188_download_FASTQ.sh   (RUNS="ERR... ERR..." downloads a subset; EMTAB11188_DIR or
#         DATA_DIR change the data folder, as in 0.config.R)

REPO_DIR="$(cd "$(dirname "$0")" && pwd)"
DATASET_DIR="${EMTAB11188_DIR:-${DATA_DIR:-$REPO_DIR/data}/E-MTAB-11188}"
FASTQ_DIR="${FASTQ_DIR:-$DATASET_DIR/fastq}"
MANIFEST="$DATASET_DIR/ENA_manifest.tsv"
RUNS="${RUNS:-ERR7425211 ERR7425212 ERR7425213 ERR7425214 ERR7425188 ERR7425190 ERR7425191 ERR7425189}"

md5_of() { if command -v md5 > /dev/null; then md5 -q "$1"; else md5sum "$1" | cut -d ' ' -f 1; fi; }

mkdir -p "$FASTQ_DIR"

##### ENA file report of the study: FASTQ URLs and MD5 checksums of every run
curl -L --fail \
  "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB48995&result=read_run&fields=run_accession,sample_accession,sample_title,library_layout,fastq_ftp,fastq_md5&format=tsv&download=true" \
  -o "$MANIFEST"

##### FASTQs of the selected runs, with MD5 check
tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r run sample title layout fastq_ftp fastq_md5; do
  [[ " $RUNS " == *" $run "* ]] || continue
  IFS=';' read -ra urls <<< "$fastq_ftp"
  IFS=';' read -ra md5s <<< "$fastq_md5"
  for i in "${!urls[@]}"; do
    url="${urls[$i]}"
    [[ "$url" == http* ]] || url="https://$url"
    outfile="$FASTQ_DIR/$(basename "$url")"
    if [ -f "$outfile" ] && [ "$(md5_of "$outfile")" = "${md5s[$i]}" ]; then
      echo "Already complete: $(basename "$outfile")"
      continue
    fi
    echo "Downloading $(basename "$outfile")"
    until curl -L --fail --connect-timeout 30 -C - -o "$outfile" "$url"; do
      echo "Transfer interrupted; retrying in 15 s"
      sleep 15
    done
    if [ "$(md5_of "$outfile")" != "${md5s[$i]}" ]; then
      echo "ERROR: MD5 check failed for $(basename "$outfile")"
      exit 1
    fi
    echo "MD5 OK: $(basename "$outfile")"
  done
done
