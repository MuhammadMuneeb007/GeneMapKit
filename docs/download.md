# Download Databases

Stage 0 downloads immutable upstream source files into `data/raw` and writes
`data/raw/manifest.json`. Compressed upstream files are stored unchanged.

## Recommended modes

| Mode | Command | Intended use |
|---|---|---|
| Default | `python download_databases.py --out data/raw` | Routine conversion snapshot |
| Core only | `python download_databases.py --out data/raw --core-only` | Small initial setup |
| Validation | `python download_databases.py --out data/raw --validation` | Default sources plus independent NCBI truth |
| Full research | `python download_databases.py --out data/raw --all --validation` | Every registered mapping, GO, PDB, and validation source |

!!! warning "Validation sources and leakage"
    Downloading a validation source does not make a benchmark independent. If that
    source is included when building the evaluated database, the result is diagnostic,
    not held-out. Exclude truth sources from benchmark database builds.

## Inspect available sources

```bash
python download_databases.py --list
python download_databases.py --coverage
```

## Download selected sources

```bash
python download_databases.py \
  --out data/raw \
  --sources hgnc mane_summary ensembl_gtf
```

## Download for a conversion

Download recommended sources that support a requested conversion:

```bash
python download_databases.py \
  --out data/raw \
  --for-conversion ensembl_transcript_id refseq_rna_id
```

Request the smallest direct source set:

```bash
python download_databases.py \
  --out data/raw \
  --for-conversion ensembl_transcript_id refseq_rna_id \
  --minimal
```

## Source tiers

| Tier | Sources | Purpose |
|---|---|---|
| Core | HGNC, NCBI Gene Info, Ensembl GTF, UniProt human, MANE | Essential conversion backbone |
| Cross-reference | Ensembl xrefs, CCDS, NCBI Gene History | Broader cross-database coverage |
| Extended | NCBI gene2ensembl, gene2refseq, gene2go, SIFTS PDB | Validation, functional, and structural mappings |

## Manifest

The manifest records:

- Source name and description
- Download URL
- Filename and byte size
- SHA-256 checksum
- Organism and taxon
- Pinned Ensembl release
- Download and verification timestamps
- Identifier namespaces provided by each source

Keep the manifest with every database build and benchmark result.

## Package CLI equivalents

```bash
python -m genemapkit.cli download --validation
python -m genemapkit.cli download --all --validation
python -m genemapkit.cli download-for \
  --input-type symbol \
  --output-type ensembl_gene_id \
  --minimal
```

## Important options

| Option | Meaning |
|---|---|
| `--ensembl-release N` | Pin the Ensembl release |
| `--force` | Download files again even if they exist |
| `--sources NAME...` | Download explicitly selected sources |
| `--validation` | Include NCBI validation sources |
| `--all` | Include every registered source |
| `--minimal` | Minimize sources for `--for-conversion` |
