# Build the Database

Stage 1 verifies the raw snapshot and creates a normalized, provenance-aware SQLite
database.

## Production build

```bash
python build_database.py \
  --raw-dir data/raw \
  --out data/build/genemapkit.db \
  --report data/build/build_report.json
```

The builder uses every parsable source recorded in the manifest unless `--source` or
`--exclude` is supplied.

## Production rules

For publication or release builds:

- Do not use `--max-records`.
- Do not use `--no-verify`.
- Do not use `--allow-missing`.
- Use a new output path while testing a changed schema.
- Archive the raw manifest and JSON build report.
- Confirm `complete_build`, `verified_snapshot`, and `missing_sources`.

## Full database

After downloading with `--all --validation`, build all available sources:

```bash
python build_database.py \
  --raw-dir data/raw \
  --out data/build/genemapkit-full.db \
  --report data/build/genemapkit-full-report.json
```

## Select or exclude sources

```bash
# Build from selected sources only.
python build_database.py \
  --raw-dir data/raw \
  --source hgnc \
  --source mane_summary \
  --source ensembl_gtf \
  --out data/build/core-example.db

# Build without functional and structural annotations.
python build_database.py \
  --raw-dir data/raw \
  --exclude ncbi_gene2go \
  --exclude sifts_pdb \
  --out data/build/no-go-or-pdb.db
```

## Held-out benchmark databases

Independent validation requires excluding the source used as truth.

```bash
# Hold out NCBI relationship truth.
python build_database.py \
  --raw-dir data/raw \
  --out data/build/genemapkit-ncbi-heldout.db \
  --report data/build/genemapkit-ncbi-heldout-report.json \
  --exclude ncbi_gene2ensembl \
  --exclude ncbi_gene2refseq

# Hold out MANE truth.
python build_database.py \
  --raw-dir data/raw \
  --out data/build/genemapkit-mane-heldout.db \
  --report data/build/genemapkit-mane-heldout-report.json \
  --exclude mane_summary
```

## Development build

Use record limits only for parser and schema development:

```bash
python build_database.py \
  --raw-dir data/raw \
  --out data/build/development.db \
  --report data/build/development-report.json \
  --max-records 1000 \
  --force
```

## Build report checks

```json
{
  "build_schema_version": "2.1.0",
  "missing_sources": [],
  "verified_snapshot": true,
  "complete_build": true
}
```

Schema 2.1 rejects Ensembl `EntrezGene_trans_name` and `MISC` records instead of
misclassifying transcript names as Entrez IDs. Rejected counts and examples appear in
`rejected_records`.

## Options

| Option | Meaning |
|---|---|
| `--raw-dir PATH` | Snapshot directory containing `manifest.json` |
| `--out PATH` | Output SQLite database |
| `--report PATH` | Output JSON build report |
| `--source NAME` | Include only this source; repeatable |
| `--exclude NAME` | Exclude this source; repeatable |
| `--force` | Replace an existing output |
| `--max-records N` | Development-only limit per source |
| `--no-verify` | Skip checksum verification; not for publication |
| `--allow-missing` | Continue with missing sources; not for publication |
| `--check VALUE INPUT OUTPUT` | Run one conversion after building |
