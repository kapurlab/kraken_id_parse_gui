#!/usr/bin/env bash
#
# add_orbivirus_to_blast.sh
#
# Make the committed BTV/EHD reference sequences available to the pipeline's
# nucleotide BLAST step (bin/kraken_id_parse.py -> Blast_Fasta), which BLASTs
# assembled scaffolds against the `blast_db` named in the preset (normally the
# large external `nt_viruses`).
#
# BLAST databases are NOT appendable in place, and `nt_viruses` is ~70 GB and
# external (it does not live in the repo). So this script, non-destructively:
#   1. checks whether the BTV/EHD GenBank accessions are already in nt_viruses,
#   2. builds a small supplement BLAST DB from the committed source FASTAs
#      (database/source_sequences/BTV_sequence.fasta + EHD_sequence.fasta), and
#   3. creates a BLAST *alias* DB that unions nt_viruses + the supplement.
# Point a preset's `blast_db` at the alias to BLAST against both at once.
#
# The supplement + alias are small; nt_viruses itself is never modified.
#
# Run inside the kraken_id_parse_env conda env (needs: makeblastdb, blastdbcmd,
# blastdb_aliastool, seqkit).
#
# Usage:
#   bin/add_orbivirus_to_blast.sh \
#       [--nt-viruses /path/to/nt_viruses] \
#       [--out-dir DIR] [--alias-name NAME] \
#       [--btv-ref FASTA] [--ehd-ref FASTA]
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# nt_viruses location: honor --nt-viruses, else $BLASTDB/nt_viruses, else the
# local dev default used elsewhere in this repo.
NT_VIRUSES="${BLASTDB:+$BLASTDB/nt_viruses}"
NT_VIRUSES="${NT_VIRUSES:-/Users/todstuber/databases/blast/nt_viruses}"
BTV_REF="$REPO_ROOT/database/source_sequences/BTV_sequence.fasta"
EHD_REF="$REPO_ROOT/database/source_sequences/EHD_sequence.fasta"
OUT_DIR=""                       # default decided after arg parse (next to nt_viruses)
SUPP_NAME="orbivirus_lab_refs"   # supplement DB name
ALIAS_NAME="nt_viruses_orbivirus"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --nt-viruses) NT_VIRUSES="$2"; shift 2 ;;
    --out-dir)    OUT_DIR="$2"; shift 2 ;;
    --alias-name) ALIAS_NAME="$2"; shift 2 ;;
    --btv-ref)    BTV_REF="$2"; shift 2 ;;
    --ehd-ref)    EHD_REF="$2"; shift 2 ;;
    -h|--help)    grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; exit 1 ;;
  esac
done

# Default the output next to nt_viruses so the alias members resolve via one
# directory ("add to the external db" in place, without touching nt_viruses).
OUT_DIR="${OUT_DIR:-$(dirname "$NT_VIRUSES")}"

log() { printf '\n==> %s\n' "$*"; }
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not on PATH (activate kraken_id_parse_env?)" >&2; exit 1; }; }
need makeblastdb; need blastdbcmd; need blastdb_aliastool; need seqkit

mkdir -p "$OUT_DIR"
OUT_DIR="$(cd "$OUT_DIR" && pwd)"   # make absolute so alias members resolve
COMBINED="$OUT_DIR/${SUPP_NAME}_input.fasta"

log "Inputs"
echo "  nt_viruses = $NT_VIRUSES"
echo "  BTV ref    = $BTV_REF"
echo "  EHD ref    = $EHD_REF"
echo "  out dir    = $OUT_DIR"
echo "  alias      = $OUT_DIR/$ALIAS_NAME"

[[ -s "$BTV_REF" || -s "$EHD_REF" ]] || { echo "ERROR: no source FASTAs found" >&2; exit 1; }

# ---------------------------------------------------------------------------
# 1. Report which GenBank accessions are already in nt_viruses (informational).
#    Lab consensuses (e.g. *_CONS) are not in GenBank and are expected to miss.
# ---------------------------------------------------------------------------
NT_PRESENT=0
if compgen -G "${NT_VIRUSES}.*" >/dev/null 2>&1 || [[ -e "${NT_VIRUSES}.nal" ]]; then
  NT_PRESENT=1
  log "Checking whether the BTV/EHD accessions are already in nt_viruses"
  # Pull GenBank-style accessions (e.g. MH845259.1) out of the headers.
  ACCS=$(grep -h '^>' "$BTV_REF" "$EHD_REF" 2>/dev/null \
         | grep -oE '[A-Z]{1,2}[0-9]{5,8}\.[0-9]+' | sort -u || true)
  if [[ -n "$ACCS" ]]; then
    n_total=$(echo "$ACCS" | wc -l | tr -d ' ')
    n_found=0
    while read -r acc; do
      [[ -z "$acc" ]] && continue
      if blastdbcmd -db "$NT_VIRUSES" -entry "$acc" -outfmt '%a' >/dev/null 2>&1; then
        n_found=$((n_found+1))
      fi
    done <<< "$ACCS"
    echo "  GenBank accessions in headers: $n_total; already in nt_viruses: $n_found"
  else
    echo "  (no GenBank-style accessions in headers; all look like lab consensuses)"
  fi
else
  echo "NOTE: nt_viruses not found at '$NT_VIRUSES' -- will build the supplement DB only." >&2
fi

# ---------------------------------------------------------------------------
# 2. Build the small supplement BLAST DB from the committed FASTAs.
#    seqkit rmdup avoids duplicate IDs that would break -parse_seqids.
# ---------------------------------------------------------------------------
log "Building supplement BLAST DB ($SUPP_NAME)"
cat ${BTV_REF:+"$BTV_REF"} ${EHD_REF:+"$EHD_REF"} 2>/dev/null \
  | seqkit rmdup -n -o "$COMBINED" 2>/dev/null
echo "  supplement sequences: $(grep -c '^>' "$COMBINED")"
makeblastdb -in "$COMBINED" -dbtype nucl -parse_seqids \
  -title "$SUPP_NAME" -out "$OUT_DIR/$SUPP_NAME"

# ---------------------------------------------------------------------------
# 3. Create the union alias (nt_viruses + supplement) when nt_viruses exists.
# ---------------------------------------------------------------------------
if [[ "$NT_PRESENT" -eq 1 ]]; then
  log "Creating alias DB unioning nt_viruses + $SUPP_NAME"
  blastdb_aliastool -dblist "$NT_VIRUSES $OUT_DIR/$SUPP_NAME" \
    -dbtype nucl -out "$OUT_DIR/$ALIAS_NAME" -title "$ALIAS_NAME"
  echo
  echo "DONE. Point a preset's blast_db at the alias to search both:"
  echo "    blast_db: \"$OUT_DIR/$ALIAS_NAME\""
else
  echo
  echo "DONE. Supplement DB built at: $OUT_DIR/$SUPP_NAME"
  echo "(no nt_viruses found, so no union alias was created)."
fi
