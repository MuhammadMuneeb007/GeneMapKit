# Convert Identifiers

GeneMapKit can convert one identifier, a CSV/TSV column, or identifiers through the
Python API.

## Convert one identifier

```bash
python -m genemapkit.cli convert-one TP53 \
  --from symbol \
  --to ensembl_gene_id \
  --policy authoritative \
  --db data/build/genemapkit.db
```

Request multiple target namespaces by repeating `--to`:

```bash
python -m genemapkit.cli convert-one TP53 \
  --from symbol \
  --to entrez_id \
  --to ensembl_gene_id \
  --to uniprot_id \
  --policy all-supported \
  --db data/build/genemapkit.db
```

## Convert a CSV or TSV column

The key input options are:

| Option | Meaning |
|---|---|
| `INPUT_FILE` | CSV/TSV file to convert |
| `--column NAME` | Column containing the input identifiers |
| `--input-type TYPE` / `-i TYPE` | Namespace of the input identifiers, for example `symbol`, `entrez_id`, or `ensembl_gene_id` |
| `--output-type TYPE` / `-o TYPE` | Target namespace; repeat this option for multiple conversions |
| `--db PATH` | Built SQLite database |
| `--format wide\|long` | `wide` gives one row per input; `long` preserves evidence/provenance rows |
| `--policy POLICY` | Mapping policy such as `all-supported`, `authoritative`, `corroborated`, or `mane-select` |

Long format preserves source-level provenance. It can emit multiple rows for the
same input and output identifier when several databases support the same mapping:

```bash
python -m genemapkit.cli convert input.csv \
  --column gene_symbol \
  --input-type symbol \
  --output-type ensembl_gene_id \
  --output-type entrez_id \
  --policy all-supported \
  --format long \
  --db data/build/genemapkit.db \
  --output mapped-long.csv
```

Wide format preserves one row per input and joins multiple outputs:

```bash
python -m genemapkit.cli convert input.csv \
  --column gene_symbol \
  --input-type symbol \
  --output-type ensembl_gene_id \
  --policy authoritative \
  --format wide \
  --output mapped-wide.csv
```

Use long format for scientific analyses that require complete provenance. Use wide
format when you want a compact table with one row per input gene. In long format,
repeated `output_id` values are expected when multiple evidence sources support
the same mapping.

## Multiple target conversions

Repeat `--output-type` to convert one input column to several identifier
namespaces in a single command:

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

For a provenance-preserving table, use the same command with `--format long`:

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
  --format long \
  --output results/all_human_genes_all_ids_long.csv
```

The same pattern works for any supported input namespace. For example, convert
Ensembl gene IDs to symbols, Entrez IDs, and UniProt IDs:

```bash
python -m genemapkit.cli convert data/sample/ensembl_gene_id_sample.csv \
  --column ensembl_gene_id \
  --input-type ensembl_gene_id \
  --output-type symbol \
  --output-type entrez_id \
  --output-type uniprot_id \
  --db data/build/genemapkit_final_v21.db \
  --policy authoritative \
  --format wide \
  --output results/ensembl_gene_to_symbol_entrez_uniprot.csv
```

## Convert one identifier with multiple outputs

For single-ID conversion, `--from` is an alias for `--input-type`, and `--to` is
an alias for `--output-type`:

```bash
python -m genemapkit.cli convert-one TP53 \
  --from symbol \
  --to hgnc_id \
  --to entrez_id \
  --to ensembl_gene_id \
  --to uniprot_id \
  --db data/build/genemapkit_final_v21.db \
  --policy authoritative
```

## Symbol policies

Symbol lookup is token-exact.

| Policy | Matching behavior |
|---|---|
| `approved` | Approved symbols only |
| `historical` | Approved and previous symbols |
| `all` | Approved, previous, and alias symbols |

```bash
python -m genemapkit.cli convert-one P53 \
  --from symbol \
  --to entrez_id \
  --symbol-policy all
```

## Historical identifiers

Retired target identifiers are excluded by default:

```bash
python -m genemapkit.cli convert-one 7157 \
  --from entrez_id \
  --to entrez_id \
  --include-retired
```

## Python API

```python
from genemapkit import GeneConverter

converter = GeneConverter(
    "data/build/genemapkit.db",
    mapping_policy="authoritative",
)

results = converter.query(
    "ENST00000269305",
    "ensembl_transcript_id",
    "refseq_rna_id",
)

for result in results:
    print(result.output_id, result.output_source, result.selection_reason)
```

Batch conversion:

```python
frame = converter.convert_batch(
    ["TP53", "BRCA1"],
    "symbol",
    ["ensembl_gene_id", "entrez_id"],
    output_format="long",
)
```

## Supported namespaces

```text
symbol, alias_symbol, previous_symbol, hgnc_id, entrez_id,
ensembl_gene_id, ensembl_transcript_id, ensembl_protein_id,
refseq_rna_id, refseq_protein_id, uniprot_id, ccds_id, omim_id,
ucsc_id, vega_id, ena_id, pdb_id, go_id
```
