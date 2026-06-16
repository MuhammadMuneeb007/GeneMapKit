# CLI Reference

GeneMapKit exposes standalone pipeline scripts and a unified package CLI.

## Unified CLI

```bash
python -m genemapkit.cli --help
```

Raw snapshot directory options are interchangeable in the package download command:
`--out-dir`, `--raw-dir`, and `--out`. The status command accepts `--raw-dir` and
`--out-dir`.

```bash
python -m genemapkit.cli status --raw-dir data/raw
```

If an existing environment reports `RequestsDependencyWarning`, synchronize its HTTP
dependencies before running live-API benchmarks:

```bash
python -m pip install -U "requests>=2.32,<3" "chardet<6"
```

| Command | Purpose |
|---|---|
| `download` | Download snapshot sources |
| `download-for` | Download sources for a requested conversion |
| `build` | Build the normalized SQLite database |
| `convert-one` | Convert one identifier |
| `convert` | Convert a CSV/TSV column |
| `validate` | Benchmark conversion relationships |
| `compare-mygene` | Compare final symbol conversions with MyGene.info |
| `benchmark-external` | Compare final symbol conversions with multiple live APIs |
| `benchmark-publication` | Run performance, coordinate, API, and held-out benchmarks |
| `benchmark-paths` | Measure all-to-all coverage, round trips, and path consistency |
| `status` | Show downloaded database status |
| `info` | Show supported identifier information |

## Standalone scripts

### Download

```bash
python download_databases.py --help
```

Common options:

```text
--out PATH
--ensembl-release N
--core-only
--all
--validation
--sources NAME...
--for-conversion INPUT OUTPUT...
--minimal
--list
--coverage
--force
```

### Build

```bash
python build_database.py --help
```

Common options:

```text
--raw-dir PATH
--out PATH
--report PATH
--source NAME
--exclude NAME
--force
--max-records N
--check VALUE INPUT_TYPE OUTPUT_TYPE
```

### Convert one

```bash
python -m genemapkit.cli convert-one --help
```

Common options:

```text
VALUE
--input-type TYPE / --from TYPE / -i TYPE
--output-type TYPE / --to TYPE / -o TYPE
--db PATH
--symbol-policy approved|historical|all
--policy all-supported|direct-only|corroborated|mane-select|authoritative
--include-retired
--provenance / --no-provenance
```

Repeat `--to` or `--output-type` to request multiple output namespaces in one
single-ID conversion.

### Convert a file

```bash
python -m genemapkit.cli convert --help
```

Common options:

```text
INPUT_FILE
--column NAME
--input-type TYPE / -i TYPE
--output-type TYPE / -o TYPE
--output PATH
--delimiter VALUE
--db PATH
--format wide|long
--separator VALUE
--symbol-policy approved|historical|all
--policy all-supported|direct-only|corroborated|mane-select|authoritative
--include-retired
--include-unmapped / --drop-unmapped
--progress / --no-progress
```

`--format long` preserves source-level provenance and may repeat an `output_id`
when multiple databases support the same mapping. `--format wide` collapses
outputs into one row per input identifier.

Repeat `--output-type` to convert a single input column to multiple target
namespaces in one run:

```bash
python -m genemapkit.cli convert data/sample/all_human_genes.csv \
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

Supported conversion namespaces are:

```text
symbol, alias_symbol, previous_symbol, hgnc_id, entrez_id,
ensembl_gene_id, ensembl_transcript_id, ensembl_protein_id,
refseq_rna_id, refseq_protein_id, uniprot_id, ccds_id, omim_id,
ucsc_id, vega_id, ena_id, pdb_id, go_id
```

### Validate

```bash
python validate_database.py --help
```

Common options:

```text
--db PATH
--raw-dir PATH
--out PATH
--source mane_summary|ncbi_gene2ensembl|ncbi_gene2refseq
--sample N
--max-examples N
--require-held-out
--policy all-supported|direct-only|corroborated|mane-select|authoritative
```

### Compare with MyGene.info

```bash
python compare_mygene.py --help
python -m genemapkit.cli compare-mygene --help
```

Common options:

```text
INPUT_FILE
--column NAME
--delimiter VALUE
--db PATH
--output-type TYPE
--policy all-supported|direct-only|corroborated|mane-select|authoritative
--out PATH
```

### Benchmark external APIs

```bash
python benchmark_external.py --help
python -m genemapkit.cli benchmark-external --help
```

Common options:

```text
INPUT_FILE
--column NAME
--delimiter VALUE
--db PATH
--provider mygene|ensembl|ncbi|uniprot|gprofiler|bridgedb
--output-type TYPE
--policy all-supported|direct-only|corroborated|mane-select|authoritative
--out PATH
--refresh-cache
```

The default run queries all six supported providers. Summary CSV, JSON, and Markdown
reports include Wilson 95% confidence intervals for concordance precision and recall.

## Progress and logs

Long-running package CLI commands show `tqdm` progress bars by default:

- `convert`: identifiers processed
- `validate`: validation sources and truth identifiers processed
- `compare-mygene`: API batches and genes compared
- `benchmark-external`: providers queried, API requests, and genes compared
- `benchmark-paths`: namespaces processed
- `benchmark-publication`: stage logs plus external and path progress

Use `--no-progress` when redirecting output to a machine-readable log or running in a
workflow system that supplies its own progress reporting. Use the top-level `--verbose`
option before the command to enable detailed package logs:

```bash
python -m genemapkit.cli --verbose benchmark-publication input.csv --no-progress
```

Standalone benchmark scripts show progress by default and accept `--no-progress`.

### Run the publication benchmark suite

```bash
python benchmark_publication.py --help
python -m genemapkit.cli benchmark-publication --help
```

This command produces local runtime and memory scaling, forced-offline operation,
coordinate concordance against live Ensembl, multi-API concordance and query timings,
independent held-out validation reports, and all-to-all path consistency.

### Benchmark conversion paths

```bash
python benchmark_paths.py --help
python -m genemapkit.cli benchmark-paths --help
```

The direct conversion coverage matrix always uses all identifiers in the database.
Round-trip and multi-step paths use a deterministic sample by default:

```bash
python benchmark_paths.py \
  --db data/build/genemapkit_final_v21.db \
  --pivot-type entrez_id \
  --sample 1000 \
  --out results/path_benchmark
```

## Exit behavior

`convert-one` exits with a non-zero status when no requested target mapping is found.
Pipeline scripts fail loudly for missing, malformed, or incompatible inputs unless an
explicit development-only option changes that behavior.
