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

Long format preserves one row per mapping and its provenance:

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

Use long format for scientific analyses that require complete provenance.

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
