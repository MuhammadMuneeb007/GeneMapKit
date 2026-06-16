# Building the GeneMapKit Identifier Database

GeneMapKit uses a two-stage data pipeline:

1. **Stage 0: download** immutable upstream files into `data/raw`.
2. **Stage 1: build** a normalized SQLite identifier index in `data/build`.

The build stage verifies the snapshot manifest, parses the downloaded sources,
normalizes versioned identifiers, stores direct source-asserted relationships, and
records source/evidence provenance.

Schema version 2 separates broad identifier groups from direct conversion
relationships. Rebuild databases created with schema version 1 before using the
current converter.

Schema version 2.1 also validates Ensembl Entrez cross-reference namespaces. Only
numeric xrefs with `db_name=EntrezGene` are stored as Entrez IDs. Ensembl
`EntrezGene_trans_name`/`MISC` rows are rejected rather than misclassified as Entrez
IDs or used to create relationships.

## Recommended Production Workflow

From the repository directory:

```bash
# Stage 0: download the normal conversion snapshot
python download_databases.py --out data/raw

# Stage 1: verify and build the complete SQLite index
python build_database.py \
  --raw-dir data/raw \
  --out data/build/genemapkit.db \
  --report data/build/build_report.json
```

The package CLI provides the same build command:

```bash
python -m genemapkit.cli build \
  --raw-dir data/raw \
  --out data/build/genemapkit.db \
  --report data/build/build_report.json
```

By default, the builder:

- Builds from every downloaded source recorded in `data/raw/manifest.json`.
- Verifies file sizes and SHA-256 checksums.
- Fails if a selected source is missing, malformed, or empty.
- Refuses to overwrite an existing database.
- Writes a complete build report.

## Validation and Full Research Builds

Download validation sources before building when preparing benchmarks:

```bash
python download_databases.py --validation

python build_database.py \
  --raw-dir data/raw \
  --out data/build/genemapkit-validation.db \
  --report data/build/validation-build-report.json
```

Because the builder uses every source present in the manifest, downloaded validation
sources are included automatically.

For every registered source, including GO and PDB mappings:

```bash
python download_databases.py --all --validation

python build_database.py \
  --raw-dir data/raw \
  --out data/build/genemapkit-full.db \
  --report data/build/full-build-report.json
```

The full build can require substantial disk space, memory, and processing time.

## Development Build

Use `--max-records` only to test parsers and schemas:

```bash
python build_database.py \
  --raw-dir data/raw \
  --out data/build/development.db \
  --report data/build/development-report.json \
  --max-records 1000 \
  --force
```

This produces an intentionally incomplete database. The report and database metadata
record `complete_build: false`. Never use a limited build for benchmarks or the paper.

## Selecting Sources

Use `--source` repeatedly to build from only specific sources:

```bash
python build_database.py \
  --source hgnc \
  --source mane_summary \
  --source ensembl_gtf \
  --out data/build/core-example.db
```

Use `--exclude` repeatedly to omit downloaded sources:

```bash
python build_database.py \
  --exclude ncbi_gene2go \
  --exclude sifts_pdb \
  --out data/build/no-functional-or-structure.db
```

Useful source names:

```text
hgnc
ncbi_gene_info
ensembl_gtf
uniprot_human
mane_summary
ensembl_entrez
ensembl_refseq
ensembl_uniprot
ccds
ncbi_gene_history
ncbi_gene2ensembl
ncbi_gene2refseq
ncbi_gene2go
sifts_pdb
```

## Build Options

| Option | Meaning |
|---|---|
| `--raw-dir PATH` | Raw snapshot directory containing `manifest.json` |
| `--out PATH` | Output SQLite database |
| `--report PATH` | JSON build report |
| `--source NAME` | Include only a source; repeatable |
| `--exclude NAME` | Exclude a source; repeatable |
| `--force` | Replace an existing output database |
| `--allow-missing` | Continue when selected sources are missing or empty |
| `--no-verify` | Skip checksum verification |
| `--max-records N` | Limit records per source for development only |
| `--check VALUE INPUT_TYPE OUTPUT_TYPE` | Run one conversion check after building |

For publication builds, do not use:

```text
--no-verify
--allow-missing
--max-records
```

## Conversion Check

Run a simple query immediately after building:

```bash
python build_database.py \
  --raw-dir data/raw \
  --out data/build/genemapkit.db \
  --report data/build/build_report.json \
  --force \
  --check TP53 symbol ensembl_gene_id
```

Other examples:

```bash
--check 7157 entrez_id refseq_rna_id
--check NM_000546.6 refseq_rna_id symbol
--check P04637 uniprot_id ensembl_protein_id
```

## Supported Identifier Types

```text
symbol
alias_symbol
previous_symbol
hgnc_id
entrez_id
ensembl_gene_id
ensembl_transcript_id
ensembl_protein_id
refseq_rna_id
refseq_protein_id
uniprot_id
ccds_id
omim_id
ucsc_id
vega_id
ena_id
pdb_id
go_id
```

Versioned Ensembl, RefSeq, and CCDS inputs are normalized for lookup while their
original values are preserved in the database.

## SQLite Schema

### `mappings`

Stores normalized and original identifiers with provenance:

```text
group_id
id_type
id_value
original_value
source
evidence
```

Symbols, aliases, and previous symbols are queryable labels but are not used to merge
groups. This prevents false merges caused by ambiguous or overlapping aliases.

### `relationships`

Stores direct source-asserted identifier links:

```text
source_type
source_value
target_type
target_value
relationship_type
source
evidence
```

The converter queries this table instead of traversing every identifier in a broad
group. This prevents one transcript from mapping to all proteins associated with its
parent gene. Transcript-to-protein links are created only from sources that explicitly
assert the pair, including Ensembl GTF/CDS, MANE, Ensembl cross-references, and NCBI
`gene2ensembl`.

### `metadata`

Records snapshot release information, selected sources, build schema version, build
timestamp, and whether the build is complete.

### `source_summary`

Records parsed rows, rows without stable anchors, and inserted mapping rows per source.

## Build Report

The JSON report contains:

- Output database SHA-256
- Number of identifier groups
- Total mapping rows
- Distinct identifier counts by type
- Parsed records by source
- Inserted mapping rows by source
- Inserted direct relationship rows by source
- Missing sources
- Rejected source-row counts, reasons, and examples
- Snapshot verification status
- Whether the build is complete

Archive the report with every benchmark and publication release.

## Safety and Reproducibility Rules

1. Keep the raw snapshot and its `manifest.json` unchanged.
2. Verify checksums for every production build.
3. Use a separate output path for different snapshots or source selections.
4. Keep build reports with generated databases.
5. Do not use development-limited databases for scientific evaluation.
6. Record the exact GeneMapKit software release used to build publication artifacts.
