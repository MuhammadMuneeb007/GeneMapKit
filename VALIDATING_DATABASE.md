# Validating the GeneMapKit Database

Validation must compare complete predicted mapping sets with relationship-level truth.
Coverage and round-trip tests alone are insufficient because an over-broad database can
return the correct mapping alongside many incorrect mappings.

## Quick Diagnostic Benchmark

```powershell
python validate_database.py `
  --db data/build/genemapkit.db `
  --raw-dir data/raw `
  --out results/validation `
  --policy all-supported `
  --sample 1000
```

Equivalent package command:

```powershell
python -m genemapkit.cli validate `
  --db data/build/genemapkit.db `
  --raw-dir data/raw `
  --out results/validation `
  --sample 1000
```

Generated files:

- `validation_summary.md`: readable metric table and structural diagnostics.
- `validation_metrics.json`: complete machine-readable results.
- `validation_benchmarks.csv`: one row per conversion benchmark.
- `validation_disagreements.csv`: expected, predicted, extra, and missing mappings.

## Independent Validation

A truth source is independent only when it was excluded from the database build.
Use `--require-held-out` to skip sources present in the build metadata:

```powershell
python -m genemapkit.cli validate `
  --require-held-out `
  --source ncbi_gene2refseq `
  --sample 1000
```

For a stronger held-out benchmark, build a separate database:

```powershell
python build_database.py `
  --raw-dir data/raw `
  --out data/build/genemapkit_heldout.db `
  --exclude ncbi_gene2ensembl `
  --exclude ncbi_gene2refseq

python validate_database.py `
  --db data/build/genemapkit_heldout.db `
  --raw-dir data/raw `
  --source ncbi_gene2ensembl `
  --source ncbi_gene2refseq `
  --require-held-out `
  --out results/validation_heldout `
  --sample 1000
```

## Metrics

- **Precision:** fraction of returned mappings present in the truth source.
- **Recall:** fraction of truth mappings returned by GeneMapKit.
- **Exact-set match:** prediction exactly equals the expected mapping set.
- **Over-mapped input:** at least one returned mapping is absent from truth.
- **Under-mapped input:** at least one expected mapping was not returned.
- **Unmapped input:** no mapping was returned.

Use sampled validation while developing. Use complete validation with `--sample 0`
only after the schema and matching logic stabilize because large NCBI sources take
substantial time to parse.

## Compare Mapping Policies

Validation accepts the same mapping policies as conversion. Write each result to a
separate directory:

```powershell
python validate_database.py `
  --db data/build/genemapkit_mane_heldout_v21.db `
  --raw-dir data/raw `
  --source mane_summary `
  --require-held-out `
  --policy all-supported `
  --out results/policy_all_supported `
  --sample 1000

python validate_database.py `
  --db data/build/genemapkit_mane_heldout_v21.db `
  --raw-dir data/raw `
  --source mane_summary `
  --require-held-out `
  --policy authoritative `
  --out results/policy_authoritative `
  --sample 1000
```

Do not evaluate `mane-select` against a MANE-held-out database: by definition that
database contains no MANE assertions, so the policy should return no mappings.
Use the full database for a diagnostic MANE-policy check, and clearly label it as
non-independent.

## External MyGene.info Concordance

After independent held-out validation, compare final symbol conversions with
MyGene.info:

```powershell
python compare_mygene.py data/sample/symbol_sample.csv `
  --column symbol `
  --db data/build/genemapkit_final_v21.db `
  --policy all-supported `
  --out results/mygene_all_supported
```

This compares complete identifier sets and writes detailed differences. MyGene.info
integrates overlapping upstream sources, so report this as external concordance rather
than independent accuracy.
