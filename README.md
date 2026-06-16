# GeneMapKit

**Reproducible, provenance-aware conversion between human gene, transcript, protein, and annotation identifier namespaces.**

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg?style=flat-square&logo=python&logoColor=white)](https://python.org)
[![License](https://img.shields.io/badge/License-MIT-green.svg?style=flat-square)](LICENSE)

**Full documentation:** [muhammadmuneeb007.github.io/GeneMapKit](https://muhammadmuneeb007.github.io/GeneMapKit/)

---

## Overview

GeneMapKit downloads authoritative public source snapshots, builds a normalized SQLite database, converts identifiers using direct source-asserted relationships, and measures accuracy against held-out truth sources. It is designed for genomics research, transcriptomics analysis, and any bioinformatics workflow that needs reproducible gene ID standardization.

**Key properties:**

- **Reproducible snapshots** — raw files are stored unchanged with URLs, SHA-256 checksums, and timestamps in `manifest.json`.
- **Direct relationships** — conversion queries source-asserted links; no speculative transitive inference.
- **Explicit policies** — choose all supported mappings, MANE Select only, corroborated, or authoritative-source preference.
- **Provenance-rich output** — long-format results show the source, evidence, policy, and reason for every mapping.
- **Independent validation** — build held-out databases to quantify precision, recall, and exact-set agreement.

---

## Supported Identifier Namespaces (18 total)

```
symbol              alias_symbol        previous_symbol
hgnc_id             entrez_id
ensembl_gene_id     ensembl_transcript_id   ensembl_protein_id
refseq_rna_id       refseq_protein_id
uniprot_id          ccds_id             omim_id
ucsc_id             vega_id             ena_id
pdb_id              go_id
```

All conversions are any-to-any across these namespaces.

---

## Install

```bash
git clone https://github.com/MuhammadMuneeb007/GeneMapKit.git
cd GeneMapKit

python -m venv .venv
# Windows:  .\.venv\Scripts\Activate.ps1
# Linux/Mac: source .venv/bin/activate

python -m pip install .
genemapkit --help
```

Requires Python 3.10 or newer.

---

## Workflow

```
Public databases
      |
      v
Stage 0: genemapkit download      ->  data/raw/ + manifest.json
      |
      v
Stage 1: genemapkit build         ->  normalized SQLite + build report
      |
      +--> Stage 2: genemapkit convert    (identifier mappings with provenance)
      |
      +--> Stage 3: genemapkit validate   (held-out accuracy benchmarks)
```

---

## Stage 0 — Download Databases

```bash
# Recommended conversion snapshot (core + cross-reference sources)
genemapkit download --out data/raw

# Also download independent benchmark-validation sources
genemapkit download --out data/raw --validation

# Every registered source including GO and PDB mappings
genemapkit download --out data/raw --all --validation

# Smallest source set for one specific conversion
genemapkit download-for \
  --input-type ensembl_transcript_id \
  --output-type refseq_rna_id \
  --out data/raw
```

| Mode | Sources | Use |
|---|---|---|
| Default | 10 core and cross-reference sources | Routine ID conversion |
| `--validation` | Default + `gene2ensembl`, `gene2refseq` | Independent benchmark validation |
| `--all --validation` | Every registered source | Full research snapshot |
| `--core-only` | 5 essential backbone sources | Minimal initial setup |

Check what was downloaded:

```bash
genemapkit status --raw-dir data/raw
genemapkit info
```

The standalone script `python download_databases.py --help` exposes the same options plus `--list`, `--coverage`, and `--sources`.

---

## Stage 1 — Build the SQLite Database

```bash
genemapkit build \
  --raw-dir data/raw \
  --out data/build/genemapkit.db \
  --report data/build/build_report.json
```

The builder verifies SHA-256 checksums, groups identifiers, stores direct source-asserted relationships, and writes a JSON build report. Confirm the report contains:

```json
{
  "build_schema_version": "2.1.0",
  "missing_sources": [],
  "verified_snapshot": true,
  "complete_build": true
}
```

**Held-out databases for independent validation:**

```bash
# Exclude MANE from the build, then validate against it
genemapkit build \
  --raw-dir data/raw \
  --out data/build/genemapkit-mane-heldout.db \
  --exclude mane_summary

# Exclude NCBI relationship sources from the build, then validate against them
genemapkit build \
  --raw-dir data/raw \
  --out data/build/genemapkit-ncbi-heldout.db \
  --exclude ncbi_gene2ensembl \
  --exclude ncbi_gene2refseq
```

---

## Stage 2 — Convert Identifiers

### Convert one identifier

```bash
genemapkit convert-one TP53 \
  --from symbol \
  --to ensembl_gene_id \
  --to entrez_id \
  --policy authoritative \
  --db data/build/genemapkit.db
```

### Convert a CSV/TSV column

Important conversion parameters:

| Option | Meaning |
|---|---|
| `--column` | Input file column containing identifiers |
| `--input-type` / `-i` | Input namespace |
| `--output-type` / `-o` | Output namespace; repeat for multiple targets |
| `--format wide\|long` | Compact one-row-per-input output or provenance-preserving evidence rows |
| `--policy` | Mapping policy |
| `--db` | Built SQLite database |

Long format (one row per mapping evidence row, full provenance):

```bash
genemapkit convert input.csv \
  --column gene_symbol \
  --input-type symbol \
  --output-type ensembl_gene_id \
  --output-type entrez_id \
  --policy all-supported \
  --format long \
  --db data/build/genemapkit.db \
  --output mapped-long.csv
```

Wide format (one row per input, multiple outputs joined with `;`):

```bash
genemapkit convert input.csv \
  --column gene_symbol \
  --input-type symbol \
  --output-type ensembl_gene_id \
  --policy authoritative \
  --format wide \
  --output mapped-wide.csv
```

Use long format when you need source-level provenance. The same `output_id` may
appear more than once if multiple databases support the same mapping. Use wide
format when you want a compact table with one row per input gene.

Convert all approved human symbols to multiple common namespaces in one command:

```bash
genemapkit convert data/sample/all_human_genes.csv \
  --column symbol \
  --input-type symbol \
  --output-type hgnc_id \
  --output-type entrez_id \
  --output-type ensembl_gene_id \
  --output-type ensembl_transcript_id \
  --output-type ensembl_protein_id \
  --output-type refseq_rna_id \
  --output-type refseq_protein_id \
  --output-type uniprot_id \
  --db data/build/genemapkit_final_v21.db \
  --policy all-supported \
  --format wide \
  --output results/all_human_genes_all_ids_wide.csv
```

### Mapping policies

| Policy | Behavior |
|---|---|
| `all-supported` | Every direct source-asserted mapping (default) |
| `direct-only` | Direct relationships only; explicit no-inference guarantee |
| `corroborated` | Mappings supported by at least two source databases |
| `mane-select` | Only mappings with MANE Select evidence; no silent fallback |
| `authoritative` | All mappings from the highest-priority available authority |

### Symbol policies

| Policy | Matching |
|---|---|
| `approved` | Approved symbols only |
| `historical` | Approved and previous symbols |
| `all` | Approved, previous, and alias symbols (default) |

---

## Stage 3 — Validate Accuracy

### Independent MANE benchmark

```bash
genemapkit validate \
  --db data/build/genemapkit-mane-heldout.db \
  --raw-dir data/raw \
  --source mane_summary \
  --require-held-out \
  --policy all-supported \
  --out results/mane-heldout
```

### Independent NCBI benchmark

```bash
genemapkit validate \
  --db data/build/genemapkit-ncbi-heldout.db \
  --raw-dir data/raw \
  --source ncbi_gene2ensembl \
  --source ncbi_gene2refseq \
  --require-held-out \
  --policy all-supported \
  --out results/ncbi-heldout
```

The `--require-held-out` flag skips any truth source that was included in the build, preventing silent contamination. Validation reports include precision, recall, exact-set agreement, over-mapping, under-mapping, and Wilson 95% confidence intervals.

### External concordance

Compare final symbol conversions with six live external services (MyGene.info, Ensembl, NCBI, UniProt, g:Profiler, BridgeDb):

```bash
genemapkit benchmark-external data/sample/all_human_genes.csv \
  --column symbol \
  --db data/build/genemapkit.db \
  --policy all-supported \
  --refresh-cache \
  --out results/external_benchmark
```

Note: external agreement measures concordance, not independent accuracy, because the services share upstream sources.

### Publication benchmark suite

Runs local performance, memory scaling, coordinate concordance, external API concordance, held-out validation aggregation, and all-to-all path consistency together:

```bash
genemapkit benchmark-publication data/sample/all_human_genes.csv \
  --column symbol \
  --db data/build/genemapkit.db \
  --gtf data/raw/Homo_sapiens.GRCh38.116.gtf.gz \
  --held-out-report results/mane-heldout \
  --held-out-report results/ncbi-heldout \
  --policy all-supported \
  --path-pivot-type entrez_id \
  --path-sample 0 \
  --path-seed 20260616 \
  --repeats 1 \
  --coordinate-batch-size 1000 \
  --out results/publication_benchmark
```

### All-to-all path consistency

```bash
genemapkit benchmark-paths \
  --db data/build/genemapkit.db \
  --pivot-type entrez_id \
  --sample 1000 \
  --out results/path_benchmark
```

---

## Python API

```python
from genemapkit import GeneConverter

converter = GeneConverter(
    "data/build/genemapkit.db",
    mapping_policy="authoritative",
)

# Query one identifier with full provenance
results = converter.query("TP53", "symbol", "ensembl_gene_id")
for r in results:
    print(r.output_id, r.output_source, r.selection_reason)

# Get all distinct output values
result = converter.convert_single("TP53", "symbol", ["ensembl_gene_id", "entrez_id"])
print(result)  # {'ensembl_gene_id': ['ENSG00000141510'], 'entrez_id': ['7157']}

# Batch conversion to a DataFrame
frame = converter.convert_batch(
    ["TP53", "BRCA1", "EGFR"],
    "symbol",
    ["ensembl_gene_id", "entrez_id"],
    output_format="long",
)
```

---

## CLI Reference

```
genemapkit --help
genemapkit COMMAND --help
```

| Command | Purpose |
|---|---|
| `download` | Download snapshot sources |
| `download-for` | Download sources for a requested conversion |
| `build` | Build the normalized SQLite database |
| `convert-one` | Convert one identifier |
| `convert` | Convert a CSV/TSV column |
| `validate` | Benchmark against held-out relationship truth |
| `compare-mygene` | Compare with MyGene.info |
| `benchmark-external` | Compare with six live external APIs |
| `benchmark-publication` | Full publication benchmark suite |
| `benchmark-paths` | All-to-all coverage and path consistency |
| `status` | Show downloaded source status |
| `info` | Show supported identifier namespaces |

Full option reference: [muhammadmuneeb007.github.io/GeneMapKit/cli-reference/](https://muhammadmuneeb007.github.io/GeneMapKit/cli-reference/)

---

## Requirements

- Python 3.10 or newer
- Dependencies (installed automatically with `pip install .`):
  - `pandas >= 2.0.0`
  - `requests >= 2.32.0`
  - `click >= 8.0.0`
  - `tqdm >= 4.64.0`
  - `psutil >= 5.9.0`
  - `colorlog >= 6.6.0`

For documentation: `pip install -r requirements-docs.txt` then `mkdocs serve`.

---

## Source Databases

| Tier | Source | Provides |
|---|---|---|
| Core | HGNC complete set | symbol, HGNC ID, Entrez, Ensembl gene, RefSeq, UniProt, CCDS, OMIM |
| Core | NCBI Gene Info | symbol, alias, Entrez, Ensembl gene, HGNC, OMIM |
| Core | Ensembl GTF (GRCh38) | Ensembl gene/transcript/protein, symbol |
| Core | UniProt human idmapping | UniProt, Entrez, Ensembl gene/transcript/protein, RefSeq, PDB |
| Core | MANE summary | Ensembl/RefSeq canonical transcript crosswalk, symbol, HGNC |
| Cross-ref | Ensembl–Entrez xref | Ensembl gene/transcript/protein, Entrez |
| Cross-ref | Ensembl–RefSeq xref | Ensembl gene/transcript/protein, RefSeq RNA/protein |
| Cross-ref | Ensembl–UniProt xref | Ensembl gene/transcript/protein, UniProt |
| Cross-ref | CCDS | CCDS ID, Entrez, symbol |
| Cross-ref | NCBI Gene History | Retired Entrez IDs, previous symbols |
| Extended | NCBI gene2ensembl | Entrez, Ensembl gene/transcript/protein |
| Extended | NCBI gene2refseq | Entrez, RefSeq RNA/protein |
| Extended | NCBI gene2go | Entrez, GO IDs |
| Extended | SIFTS PDB | PDB IDs, UniProt |

---

## License

MIT — see [LICENSE](LICENSE).

---

## Author

- **Muhammad Muneeb**, The University of Queensland
- Email: [m.muneeb@uq.edu.au](mailto:m.muneeb@uq.edu.au)
- GitHub: [MuhammadMuneeb007](https://github.com/MuhammadMuneeb007/)
- Google Scholar: [scholar.google.com](https://scholar.google.com/citations?hl=en&user=X0xdltIAAAAJ&view_op=list_works&sortby=pubdate)
- Supervisor: [David Ascher](https://scmb.uq.edu.au/profile/8654/david-ascher), [BioSig Lab](https://biosig.lab.uq.edu.au/)

---

## Citation

```bibtex
@software{muneeb2025genemapkit,
  title={GeneMapKit: A Comprehensive Gene Identifier Mapping Toolkit},
  author={Muneeb, Muhammad and Ascher, David B.},
  year={2025},
  institution={The University of Queensland},
  url={https://github.com/MuhammadMuneeb007/GeneMapKit},
  doi={https://doi.org/10.5281/zenodo.16731821}
}
```
