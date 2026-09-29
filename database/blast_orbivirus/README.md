# Orbivirus BLAST references

The pipeline's BLAST step (`bin/kraken_id_parse.py` → `Blast_Fasta`) BLASTs assembled
scaffolds against the `blast_db` named in the active preset — normally NCBI **`nt_viruses`**.
That database is ~70 GB and lives **outside the repo** (it cannot be committed). This folder
holds the small, committed source FASTAs and a script that makes the lab's BTV/EHD reference
sequences searchable alongside `nt_viruses`.

## How `nt_viruses` itself is obtained

`nt_viruses` is a standard NCBI BLAST database, downloaded with the BLAST+ helper:
```bash
update_blastdb.pl --decompress nt_viruses        # into $BLASTDB (e.g. /…/databases/blast)
export BLASTDB=/path/to/databases/blast
```
On the HPC it is at `/software/public/databases/BLAST/db/nt_viruses`.

## Adding the BTV/EHD reference sequences

BLAST databases are **not appendable in place**, and `nt_viruses` is large/external, so
`bin/add_orbivirus_to_blast.sh` adds the sequences **non-destructively**:

1. Reports how many of the BTV/EHD **GenBank accessions** are already in `nt_viruses`
   (most are — in testing, 66/67). Lab consensuses (`*_CONS`) are not in GenBank and are the
   sequences this supplement actually contributes.
2. Builds a small **supplement** BLAST DB (`orbivirus_lab_refs`) from
   `database/source_sequences/BTV_sequence.fasta` + `EHD_sequence.fasta`.
3. Creates a BLAST **alias** DB (`nt_viruses_orbivirus`) that unions `nt_viruses` + the
   supplement. `nt_viruses` is never modified.

```bash
conda activate kraken_id_parse_env
export BLASTDB=/path/to/databases/blast            # where nt_viruses lives
bin/add_orbivirus_to_blast.sh                      # writes supplement+alias next to nt_viruses
# or keep the small parts in-repo and reference the external nt_viruses by path:
bin/add_orbivirus_to_blast.sh --out-dir database/blast_orbivirus
```

Then point the preset's `blast_db` at the alias, e.g.:
```yaml
blast_db: "/path/to/databases/blast/nt_viruses_orbivirus"
```

Built DB files in this directory are **git-ignored** (see `.gitignore`); only the source
FASTAs and the build script are committed.

## How to add additional genomes later

Append sequences to the matching committed FASTA, then rebuild the supplement + alias:
```bash
cat new_seqs.fasta >> database/source_sequences/BTV_sequence.fasta   # (or EHD_sequence.fasta)
bin/add_orbivirus_to_blast.sh
```
`makeblastdb -parse_seqids` requires unique sequence IDs; the script runs `seqkit rmdup -n`
first to drop duplicate-named records.

## Note on serotyping vs. this BLAST DB

BTV **serotyping** (`bin/btv_serotyping.py`) does **not** use `nt_viruses` or this alias — it
BLASTs VP2/VP5 with `blastn -subject` directly against the curated panel in the repo's
top-level `reference_sequences/`. This folder is only about the scaffold-identification BLAST.
