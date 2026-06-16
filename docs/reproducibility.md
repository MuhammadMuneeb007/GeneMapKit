# Reproducibility

GeneMapKit is designed so a database and benchmark can be reconstructed from source
snapshots, software code, and reports.

## Archive these artifacts

For every release or paper benchmark, archive:

1. Raw `manifest.json`
2. Raw source files or a durable snapshot location
3. GeneMapKit software version or Git commit
4. Build command
5. JSON build report
6. Built database SHA-256
7. Validation commands
8. Validation reports and disagreement examples
9. Mapping policy used for each result

## Production build checklist

- [ ] The Ensembl release is pinned.
- [ ] Raw source checksums are verified.
- [ ] `missing_sources` is empty.
- [ ] `verified_snapshot` is `true`.
- [ ] `complete_build` is `true`.
- [ ] No development-only flags were used.
- [ ] Rejected source records were reviewed.
- [ ] The database SHA-256 was recorded.

## Benchmark checklist

- [ ] Truth sources are excluded from the evaluated build.
- [ ] Validation rows report `independent: true`.
- [ ] Precision, recall, exact-set agreement, over-mapping, and under-mapping are reported.
- [ ] Policy comparisons use the same truth inputs.
- [ ] Disagreement examples are manually reviewed before claims are written.
- [ ] Diagnostic and independent results are clearly separated.

## Publication commands

Use complete validation only after sampled development benchmarks stabilize:

```bash
python validate_database.py \
  --db data/build/genemapkit-ncbi-heldout.db \
  --raw-dir data/raw \
  --source ncbi_gene2ensembl \
  --source ncbi_gene2refseq \
  --require-held-out \
  --policy authoritative \
  --sample 0 \
  --out results/publication-ncbi-authoritative
```

## Avoid leakage

If a truth source is present in the database, matching it does not demonstrate
independent accuracy. The validator records selected build sources and labels each
benchmark as independent or non-independent.

## Database redistribution

Do not assume the combined SQLite database can be redistributed without restrictions.
Review licenses, attribution requirements, and redistribution terms for every upstream
source. Source manifests, build scripts, and checksums can still support reproducible
local construction.

## Documentation site

Build the documentation locally:

```bash
python -m pip install -r requirements-docs.txt
mkdocs serve
```

The GitHub Pages workflow builds documentation with strict link and configuration
checks whenever documentation changes are pushed to `main`.

## Publish GitHub Pages

The repository includes `.github/workflows/pages.yml`. To publish the site:

1. Push the documentation files to the `main` branch.
2. Open the GitHub repository.
3. Go to **Settings > Pages**.
4. Under **Build and deployment**, select **GitHub Actions** as the source.
5. Run the **Deploy documentation** workflow, or push another documentation change.

After deployment, the site is available at:

```text
https://muhammadmuneeb007.github.io/GeneMapKit/
```
