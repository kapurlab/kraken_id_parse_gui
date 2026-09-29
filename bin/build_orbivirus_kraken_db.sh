#!/usr/bin/env bash
#
# build_orbivirus_kraken_db.sh
#
# Build a custom, UNCAPPED Kraken2 database = RefSeq viral genomes
# + heavy enrichment with the entire NCBI *Orbivirus* genus + (optionally)
# the lab's own validated reference-guided consensus sequences.
#
# Why: the pipeline extracts the whole "Orbivirus" genus clade
# (extract_kraken_reads.py --include-children) before assembly/serotyping,
# so the number of reads that reach VP2 assembly is a direct function of how
# many reads Kraken classifies as Orbivirus -- i.e. of database content.
# k2_standard_08gb under-extracts because it is downsampled to an 8 GB cap and
# built only from RefSeq (sparse orbivirus diversity). This DB is uncapped and
# orbivirus-rich, so it recovers far more BTV/EHD reads/segments.
#
# Downloads go over HTTPS via NCBI `datasets` + a `curl` of the taxonomy dump,
# and `kraken2-build` is used ONLY for --add-to-library and --build. This
# deliberately avoids kraken2's bundled rsync/FTP downloaders, which fail on
# networks (and this Mac) where NCBI rsync (port 873) and FTP are blocked.
#
# Once verified locally, the database directory is self-contained (*.k2d +
# taxonomy/) and can be rsync'd to the HPC; point a preset's `kraken_db` at it.
#
# Run inside the `kraken_id_parse_env` conda env (needs: kraken2/kraken2-build,
# blast/dustmasker, seqkit, ncbi-datasets-cli [datasets/dataformat], curl).
#
# Usage:
#   bin/build_orbivirus_kraken_db.sh \
#       [--db DIR] [--panel-dir DIR] [--panel-taxid N] \
#       [--threads N] [--no-masking] \
#       [--test-r1 FASTQ --test-r2 FASTQ]
#
set -euo pipefail

# Resolve the repo root from this script's location so the default DB and the
# committed reference FASTAs are found regardless of the caller's CWD.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# ---------------------------------------------------------------------------
# Parameters (override via flags)
# ---------------------------------------------------------------------------
# Default: build INTO the repo so the DB is self-contained (binaries are
# git-ignored via database/kraken_orbivirus/.gitignore; only the source FASTAs
# and this script are committed). Override with --db to build elsewhere.
DB="$REPO_ROOT/database/kraken_orbivirus"
PANEL_DIR=""                 # dir of *_reference_guided.fasta lab consensuses (optional)
PANEL_TAXID="40051"          # Bluetongue virus (genus Orbivirus=10892, EHDV=40054)
# Committed lab reference panels (BTV/EHD genome segments) injected into the DB
# so orbivirus/EHDV recall is reproducible from the repo alone.
BTV_REF_FASTA="$REPO_ROOT/database/source_sequences/BTV_sequence.fasta"
EHD_REF_FASTA="$REPO_ROOT/database/source_sequences/EHD_sequence.fasta"
BTV_TAXID="40051"            # Bluetongue virus (species)
EHD_TAXID="40054"            # Epizootic hemorrhagic disease virus (species)
THREADS="$( (command -v sysctl >/dev/null 2>&1 && sysctl -n hw.ncpu) || nproc || echo 4)"
NO_MASKING=0                 # set with --no-masking only if dustmasker masking fails
TEST_R1=""
TEST_R2=""
VIRUSES_TAXID="10239"        # coarse "Viruses" taxid applied to the RefSeq viral background
ORBI_GENUS_TAXID="10892"     # genus Orbivirus -- fallback taxid for unmapped enrichment seqs
TAXDUMP_URL="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
VIRAL_REFSEQ_URL="https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --db)          DB="$2"; shift 2 ;;
    --panel-dir)   PANEL_DIR="$2"; shift 2 ;;
    --panel-taxid) PANEL_TAXID="$2"; shift 2 ;;
    --btv-ref)     BTV_REF_FASTA="$2"; shift 2 ;;
    --ehd-ref)     EHD_REF_FASTA="$2"; shift 2 ;;
    --no-lab-refs) BTV_REF_FASTA=""; EHD_REF_FASTA=""; shift ;;
    --threads)     THREADS="$2"; shift 2 ;;
    --no-masking)  NO_MASKING=1; shift ;;
    --test-r1)     TEST_R1="$2"; shift 2 ;;
    --test-r2)     TEST_R2="$2"; shift 2 ;;
    -h|--help)     grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; exit 1 ;;
  esac
done

WORK="$DB/_build_tmp"
MASK_FLAG=""
[[ "$NO_MASKING" -eq 1 ]] && MASK_FLAG="--no-masking"

log() { printf '\n[%s] %s\n' "$(date '+%H:%M:%S')" "$*"; }
die() { echo "ERROR: $*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "required tool '$1' not found in PATH (activate the kraken_id_parse_env env?)"; }

# Download all records of an NCBI Virus taxon (HTTPS) and emit a FASTA whose
# headers are >ACCESSION|kraken:taxid|<taxid>, taxid from datasets metadata
# (falling back to $3 so no sequence is dropped). $1=taxon $2=extra datasets
# flags $3=fallback-taxid $4=output-fasta $5=required(1)/optional(0).
# Returns nonzero (instead of dying) when optional and the download fails after
# retries -- NCBI HTTP/2 streams reset on large transfers, so we retry with
# HTTP/1.1 forced (GODEBUG=http2client=0).
fetch_and_tag() {
  local taxon="$1" extra="$2" fallback="$3" out="$4" required="${5:-1}"
  local d="$WORK/dl_${taxon}"
  local fna="$d/ncbi_dataset/data/genomic.fna"
  local rep="$d/ncbi_dataset/data/data_report.jsonl"
  local attempt
  if [[ ! -s "$fna" ]]; then
    for attempt in 1 2 3 4; do
      rm -rf "$d" "$WORK/${taxon}.zip"
      if GODEBUG=http2client=0 datasets download virus genome taxon "$taxon" \
           $extra --include genome --no-progressbar --filename "$WORK/${taxon}.zip" \
         && unzip -q -o "$WORK/${taxon}.zip" -d "$d" && [[ -s "$fna" ]]; then
        break
      fi
      echo "  download attempt $attempt for '$taxon' failed; retrying..." >&2
      sleep $((attempt * 15))
    done
  fi
  if [[ ! -s "$fna" ]]; then
    [[ "$required" -eq 1 ]] && die "no genomic.fna produced for taxon '$taxon' after retries"
    echo "WARNING: download for '$taxon' failed after retries; continuing without it." >&2
    return 1
  fi
  echo "  $taxon: $(grep -c '^>' "$fna") sequences downloaded"
  [[ -s "$rep" ]] || die "datasets data_report.jsonl missing for '$taxon'"
  dataformat tsv virus-genome --inputfile "$rep" \
    --fields accession,virus-tax-id > "$d/acc2taxid.tsv"
  awk -v fb="$fallback" '
    NR==FNR { if (FNR>1 && $1!="") map[$1]=$2; next }
    /^>/    { acc=substr($1,2); t=(acc in map)?map[acc]:fb;
              print ">"acc"|kraken:taxid|"t; next }
            { print }
  ' "$d/acc2taxid.tsv" "$fna" > "$out"
  return 0
}

# Fetch the RefSeq viral release as one HTTPS file (curl is far more robust than
# the datasets API for the whole-viral set, which resets mid-stream). Exclude any
# accession present in $2 (the orbivirus set) so orbivirus minimizers are never
# LCA-coarsened up to Viruses, then coarse-tag the remainder with taxid 10239.
# $1=output-fasta $2=orbivirus-accession-list. Returns nonzero on failure.
fetch_viral_refseq() {
  local out="$1" orbi_accs="$2"
  local gz="$WORK/viral.1.1.genomic.fna.gz" fna="$WORK/viral_refseq.fna"
  if [[ ! -s "$fna" ]]; then
    curl -fSL --retry 5 --retry-delay 10 -o "$gz" "$VIRAL_REFSEQ_URL" || return 1
    gunzip -c "$gz" > "$fna" || return 1
  fi
  [[ -s "$fna" ]] || return 1
  echo "  RefSeq viral: $(grep -c '^>' "$fna") sequences downloaded"
  local filtered="$WORK/viral_no_orbi.fna"
  if [[ -s "$orbi_accs" ]]; then
    seqkit grep -v -f "$orbi_accs" "$fna" > "$filtered" 2>/dev/null || cp "$fna" "$filtered"
  else
    cp "$fna" "$filtered"
  fi
  awk -v t="$VIRUSES_TAXID" '
    /^>/ { acc=substr($1,2); print ">"acc"|kraken:taxid|"t; next }
         { print }
  ' "$filtered" > "$out"
  return 0
}

# Tag a plain reference FASTA (headers carry no taxid) with an explicit
# |kraken:taxid|<taxid> and add it to the library. Used for the committed
# BTV/EHD lab reference panels. $1=fasta $2=taxid $3=label.
inject_ref_fasta() {
  local f="$1" taxid="$2" label="$3"
  [[ -n "$f" ]] || return 0
  if [[ ! -s "$f" ]]; then
    echo "NOTE: $label reference '$f' not found -- skipping." >&2
    return 0
  fi
  local tagged="$WORK/$(basename "$f").tagged.fna"
  # sub(/\r$/) strips CRLF carriage returns: otherwise the \r lands inside the
  # tagged header (>id\r|kraken:taxid|N) and kraken2's scan_fasta_file.pl can't
  # parse the taxid (it treats \r as whitespace ending the seqid).
  awk -v t="$taxid" '
    { sub(/\r$/, "") }
    /^>/ { acc=substr($1,2); print ">"acc"|kraken:taxid|"t; next }
         { print }
  ' "$f" > "$tagged"
  echo "  $label: $(grep -c '^>' "$tagged") sequences tagged taxid $taxid"
  kraken2-build --add-to-library "$tagged" $MASK_FLAG --db "$DB"
}

# ---------------------------------------------------------------------------
# Pre-flight
# ---------------------------------------------------------------------------
log "Pre-flight checks"
need kraken2-build; need kraken2; need seqkit; need dustmasker
need datasets;      need dataformat; need bracken-build
need curl;          need unzip
mkdir -p "$DB" "$WORK" "$DB/taxonomy"
echo "DB          = $DB"
echo "PANEL_DIR   = ${PANEL_DIR:-<none>}"
echo "BTV_REF     = ${BTV_REF_FASTA:-<none>} (taxid $BTV_TAXID)"
echo "EHD_REF     = ${EHD_REF_FASTA:-<none>} (taxid $EHD_TAXID)"
echo "THREADS     = $THREADS"
echo "MASKING     = $([[ $NO_MASKING -eq 1 ]] && echo off || echo on)"
echo "TRANSPORT   = HTTPS (datasets + curl); kraken2-build used only for add/build"

# ---------------------------------------------------------------------------
# Step 1 - taxonomy (HTTPS; nodes.dmp/names.dmp only -- no accession2taxid needed
#          because every added sequence carries an explicit |kraken:taxid|)
# ---------------------------------------------------------------------------
if [[ ! -s "$DB/taxonomy/nodes.dmp" || ! -s "$DB/taxonomy/names.dmp" ]]; then
  log "Step 1/6: downloading NCBI taxonomy dump over HTTPS"
  curl -fSL --retry 3 -o "$WORK/taxdump.tar.gz" "$TAXDUMP_URL"
  tar -xzf "$WORK/taxdump.tar.gz" -C "$DB/taxonomy" nodes.dmp names.dmp
  [[ -s "$DB/taxonomy/nodes.dmp" && -s "$DB/taxonomy/names.dmp" ]] || die "taxonomy extract failed"
else
  log "Step 1/6: taxonomy already present -- skipping"
fi

# ---------------------------------------------------------------------------
# Step 2 - Orbivirus genus enrichment (HTTPS via datasets; all GenBank + RefSeq).
#          Done first so its accession list can be excluded from the viral set.
# ---------------------------------------------------------------------------
log "Step 2/6: fetching entire NCBI Orbivirus genus (species-level taxids)"
fetch_and_tag orbivirus "" "$ORBI_GENUS_TAXID" "$WORK/orbi_tagged.fna" 1
ORBI_ACCS="$WORK/orbi_accs.txt"
tail -n +2 "$WORK/dl_orbivirus/acc2taxid.tsv" | cut -f1 | sort -u > "$ORBI_ACCS"

# ---------------------------------------------------------------------------
# Step 3 - RefSeq viral background (HTTPS via curl; coarse "Viruses" taxid).
#          Best-effort: orbivirus recall is preserved even if this is skipped.
# ---------------------------------------------------------------------------
log "Step 3/6: fetching RefSeq viral release for non-orbivirus context"
VIRAL_OK=1
if ! fetch_viral_refseq "$WORK/viral_tagged.fna" "$ORBI_ACCS"; then
  VIRAL_OK=0
fi

# Combine orbivirus (species-tagged) + viral background (orbivirus-excluded,
# coarse-tagged), then add to the library in one pass.
if [[ "$VIRAL_OK" -eq 1 ]]; then
  log "Adding orbivirus + viral background to library"
  cat "$WORK/orbi_tagged.fna" "$WORK/viral_tagged.fna" > "$WORK/viral_orbi_raw.fna"
else
  echo "NOTE: viral background unavailable; building Orbivirus-only base (primary goal preserved)." >&2
  cp "$WORK/orbi_tagged.fna" "$WORK/viral_orbi_raw.fna"
fi
seqkit rmdup -n "$WORK/viral_orbi_raw.fna" -o "$WORK/viral_orbi.fna" 2> "$WORK/rmdup.log" \
  || cp "$WORK/viral_orbi_raw.fna" "$WORK/viral_orbi.fna"
echo "  combined unique sequences: $(grep -c '^>' "$WORK/viral_orbi.fna")"
kraken2-build --add-to-library "$WORK/viral_orbi.fna" $MASK_FLAG --db "$DB"

# ---------------------------------------------------------------------------
# Step 4 - lab panel injection (optional)
# ---------------------------------------------------------------------------
if [[ -n "$PANEL_DIR" ]]; then
  shopt -s nullglob
  PANEL_FILES=("$PANEL_DIR"/*reference_guided.fasta)
  [[ ${#PANEL_FILES[@]} -eq 0 ]] && PANEL_FILES=("$PANEL_DIR"/*.fasta)
  shopt -u nullglob
  if [[ ${#PANEL_FILES[@]} -eq 0 ]]; then
    echo "WARNING: --panel-dir '$PANEL_DIR' has no FASTA files; skipping lab panel." >&2
  else
    log "Step 4/6: injecting ${#PANEL_FILES[@]} lab consensus file(s) (taxid $PANEL_TAXID)"
    LAB_TAGGED="$WORK/lab_panel_tagged.fna"
    : > "$LAB_TAGGED"
    for f in "${PANEL_FILES[@]}"; do
      base="$(basename "$f" | sed 's/\.[^.]*$//' | tr ' .' '__')"
      awk -v t="$PANEL_TAXID" -v s="$base" '
        /^>/ { n++; print ">"s"_"n"|kraken:taxid|"t; next }
             { print }
      ' "$f" >> "$LAB_TAGGED"
    done
    kraken2-build --add-to-library "$LAB_TAGGED" $MASK_FLAG --db "$DB"
  fi
else
  log "Step 4/6: no --panel-dir given -- skipping lab panel injection"
fi

# ---------------------------------------------------------------------------
# Step 4b - committed BTV/EHD reference panels (database/source_sequences/).
#           These ship in the repo so orbivirus/EHDV recall is reproducible
#           without any external folder. BTV -> taxid 40051, EHD -> 40054.
# ---------------------------------------------------------------------------
log "Step 4b/6: injecting committed BTV/EHD reference panels"
inject_ref_fasta "$BTV_REF_FASTA" "$BTV_TAXID" "BTV_sequence"
inject_ref_fasta "$EHD_REF_FASTA" "$EHD_TAXID" "EHD_sequence"

# ---------------------------------------------------------------------------
# Step 5 - build index, then Bracken distributions, then clean.
#          Bracken must run BEFORE --clean: it needs library/ and taxonomy/,
#          both of which --clean removes. The pipeline calls `bracken ... -r 250`
#          (bin/kraken.py), so database250mers.kmer_distrib is REQUIRED; the
#          other lengths are built for flexibility.
# ---------------------------------------------------------------------------
log "Step 5/6: building index (uncapped) with $THREADS threads"
kraken2-build --build --threads "$THREADS" --db "$DB"

log "Step 5/6: building Bracken read-length distributions (250 required; 100/150/200 extra)"
for L in 250 150 100 200; do
  bracken-build -d "$DB" -t "$THREADS" -k 35 -l "$L" \
    || echo "WARNING: bracken-build -l $L failed (pipeline needs at least -l 250)" >&2
done
[[ -s "$DB/database250mers.kmer_distrib" ]] \
  || die "database250mers.kmer_distrib missing -- the pipeline's 'bracken -r 250' step would fail"

log "Step 5/6: cleaning intermediate library files"
kraken2-build --clean --db "$DB"

# ---------------------------------------------------------------------------
# Step 6 - Krona taxonomy link + sanity check
# ---------------------------------------------------------------------------
log "Step 6/6: final checks"
[[ -s "$DB/hash.k2d" ]] || die "build finished but $DB/hash.k2d is missing"
echo "  DB size: $(du -sh "$DB" | cut -f1)"

# Repoint Krona's taxonomy at this DB (mirrors README setup), best-effort.
if command -v ktUpdateTaxonomy.sh >/dev/null 2>&1; then
  KRONA_TAX="$(dirname "$(command -v ktUpdateTaxonomy.sh)")/../opt/krona/taxonomy"
  if [[ -d "$(dirname "$KRONA_TAX")" ]]; then
    rm -rf "$KRONA_TAX" && ln -s "$DB" "$KRONA_TAX" || true
    echo "  Krona taxonomy linked -> $DB (run ktUpdateTaxonomy.sh if needed)"
  fi
fi

if [[ -n "$TEST_R1" && -n "$TEST_R2" ]]; then
  log "Sanity check: classifying test reads against the new DB"
  kraken2 --db "$DB" --threads "$THREADS" --paired "$TEST_R1" "$TEST_R2" \
    --output /dev/null --report "$WORK/sanity_report.txt"
  echo "  --- Orbivirus clade in sanity report ---"
  grep -i 'orbivirus' "$WORK/sanity_report.txt" || echo "  (no Orbivirus lines -- check input)"
fi

log "DONE. Custom DB ready at: $DB"
echo "Portable footprint: only *.k2d are needed at runtime ($(du -ch "$DB"/*.k2d | tail -1 | cut -f1))."
echo "The _build_tmp/ cache ($(du -sh "$DB/_build_tmp" 2>/dev/null | cut -f1)) is safe to delete before transfer:"
echo "  rsync -av --exclude _build_tmp '$DB' user@hpc:/dest/path/"
echo "Next: run the pipeline with --preset macos_dev_orbivirus_custom, then add the HPC preset."
