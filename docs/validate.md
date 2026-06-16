# Validate Accuracy

Validation reuses one immutable read-only SQLite connection per conversion benchmark.
On HPC systems, copying the built database to node-local scratch before a complete
held-out run can further reduce shared-filesystem latency.

Validation compares complete predicted mapping sets with relationship-level truth.
Coverage and round-trip tests alone cannot detect over-mapping.

## Diagnostic validation

```bash
python validate_database.py \
  --db data/build/genemapkit-full.db \
  --raw-dir data/raw \
  --out results/validation-diagnostic \
  --policy all-supported \
  --sample 1000
```

If truth sources were included in the database build, rows are marked
`independent: false`. These results diagnose behavior but must not be presented as
independent accuracy estimates.

## Independent NCBI benchmark

First build a database that excludes the NCBI truth sources:

```bash
python build_database.py \
  --raw-dir data/raw \
  --out data/build/genemapkit-ncbi-heldout.db \
  --report data/build/genemapkit-ncbi-heldout-report.json \
  --exclude ncbi_gene2ensembl \
  --exclude ncbi_gene2refseq
```

Then validate:

```bash
python validate_database.py \
  --db data/build/genemapkit-ncbi-heldout.db \
  --raw-dir data/raw \
  --source ncbi_gene2ensembl \
  --source ncbi_gene2refseq \
  --require-held-out \
  --policy authoritative \
  --sample 1000 \
  --out results/ncbi-heldout-authoritative
```

## Independent MANE benchmark

```bash
python build_database.py \
  --raw-dir data/raw \
  --out data/build/genemapkit-mane-heldout.db \
  --report data/build/genemapkit-mane-heldout-report.json \
  --exclude mane_summary

python validate_database.py \
  --db data/build/genemapkit-mane-heldout.db \
  --raw-dir data/raw \
  --source mane_summary \
  --require-held-out \
  --policy authoritative \
  --sample 1000 \
  --out results/mane-heldout-authoritative
```

Do not evaluate `mane-select` against a MANE-held-out database. The policy requires
MANE evidence, which the held-out database intentionally lacks.

## Compare policies

Run each policy against the same held-out database and write results to separate
directories:

```bash
python validate_database.py \
  --db data/build/genemapkit-ncbi-heldout.db \
  --raw-dir data/raw \
  --source ncbi_gene2ensembl \
  --source ncbi_gene2refseq \
  --require-held-out \
  --policy all-supported \
  --sample 1000 \
  --out results/ncbi-all-supported

python validate_database.py \
  --db data/build/genemapkit-ncbi-heldout.db \
  --raw-dir data/raw \
  --source ncbi_gene2ensembl \
  --source ncbi_gene2refseq \
  --require-held-out \
  --policy authoritative \
  --sample 1000 \
  --out results/ncbi-authoritative
```

## Metrics

| Metric | Meaning |
|---|---|
| Precision | Fraction of returned mappings present in truth |
| Recall | Fraction of truth mappings returned |
| Exact-set match | Predicted set exactly equals expected set |
| Over-mapped | At least one returned mapping is absent from truth |
| Under-mapped | At least one truth mapping is missing |
| Unmapped | No mapping was returned |

## Generated reports

- `validation_summary.md`: readable benchmark table and diagnostics
- `validation_metrics.json`: complete machine-readable report
- `validation_benchmarks.csv`: one row per conversion benchmark
- `validation_disagreements.csv`: expected, predicted, extra, and missing mappings

Use `--sample 0` only after the logic stabilizes; it evaluates complete truth sources.

## External concordance

After held-out validation, compare symbol conversions from the final full database with
MyGene.info. This is useful as an external integration check, but it is not independent
ground truth because MyGene.info integrates several overlapping public databases.

See [Compare with MyGene.info](external-concordance.md).
