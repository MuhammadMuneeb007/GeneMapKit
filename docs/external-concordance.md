# Compare with MyGene.info

The final full GeneMapKit database can be compared with live MyGene.info annotations
for human gene symbols.

!!! important "Concordance is not accuracy"
    MyGene.info integrates data from sources including NCBI and Ensembl that overlap
    with GeneMapKit. It is not an independent gold standard. Report this analysis as
    **external concordance**, not validation accuracy.

## Run the comparison

```bash
python compare_mygene.py data/sample/symbol_sample.csv \
  --column symbol \
  --db data/build/genemapkit_final_v21.db \
  --policy all-supported \
  --out results/mygene_all_supported
```

Compare mappings supported by at least two GeneMapKit sources:

```bash
python compare_mygene.py data/sample/symbol_sample.csv \
  --column symbol \
  --db data/build/genemapkit_final_v21.db \
  --policy corroborated \
  --out results/mygene_corroborated
```

Equivalent package command:

```bash
python -m genemapkit.cli compare-mygene data/sample/symbol_sample.csv \
  --column symbol \
  --db data/build/genemapkit_final_v21.db \
  --policy authoritative \
  --out results/mygene_authoritative
```

## Compared namespaces

By default, the command compares:

- Ensembl gene IDs
- Entrez Gene IDs
- RefSeq RNA accessions
- RefSeq protein accessions
- UniProt accessions

Select specific outputs by repeating `--output-type`:

```bash
python compare_mygene.py data/sample/symbol_sample.csv \
  --output-type ensembl_gene_id \
  --output-type entrez_id
```

## Set-based comparison

The command compares complete output sets instead of taking the first identifier.
Each symbol and output namespace receives one status:

| Status | Meaning |
|---|---|
| `exact_set` | GeneMapKit and MyGene return the same complete set |
| `overlap` | Both return results and at least one mapping agrees |
| `disjoint` | Both return results, but none agree |
| `genemapkit_only` | Only GeneMapKit returns mappings |
| `mygene_only` | Only MyGene returns mappings |
| `both_unmapped` | Neither returns a mapping |

## Generated reports

- `mygene_concordance_summary.md`: readable summary table
- `mygene_concordance_summary.csv`: one summary row per namespace
- `mygene_concordance_details.csv`: complete per-symbol set differences
- `mygene_concordance.json`: complete metadata and results

The report records the GeneMapKit database metadata, mapping policy, live MyGene
endpoint, comparison timestamp, and a warning against interpreting concordance as
independent accuracy.

## Recommended publication evaluation

Use three complementary analyses:

1. **Held-out MANE/NCBI validation:** independent relationship-level accuracy.
2. **Multi-service external concordance:** comparison with live MyGene.info, Ensembl,
   NCBI, UniProt, g:Profiler, and BridgeDb services.
3. **Source ablation:** measure what coverage is uniquely contributed by each database.

## Multi-service external benchmark

Run the final database against all supported external services:

```bash
python benchmark_external.py data/sample/symbol_sample.csv \
  --column symbol \
  --db data/build/genemapkit_final_v21.db \
  --policy all-supported \
  --out results/external_benchmark
```

The command prints combined and per-provider tables. It caches normalized API responses
under `results/external_benchmark/cache`. Use `--refresh-cache` to deliberately query
every service again.

Select particular services or namespaces by repeating options:

```bash
python benchmark_external.py data/sample/symbol_sample.csv \
  --provider ensembl \
  --provider uniprot \
  --output-type ensembl_gene_id \
  --output-type uniprot_id
```

| Service | Supported comparisons |
|---|---|
| MyGene.info | Ensembl gene, Entrez, RefSeq RNA/protein, UniProt |
| Ensembl REST | Ensembl gene, transcript, and protein |
| NCBI E-utilities | Entrez Gene |
| UniProt REST | UniProt, linked Ensembl gene/transcript/protein, linked RefSeq RNA/protein |
| g:Profiler conversion API | Ensembl gene/transcript/protein, Entrez, RefSeq RNA/protein, UniProt |
| BridgeDb web service | Ensembl gene/transcript/protein, Entrez, RefSeq RNA/protein, UniProt |

The Ensembl comparator uses the batch symbol-lookup endpoint. It requests up to 1,000
symbols per HTTP request and disables expensive transcript expansion unless Ensembl
transcript or protein outputs are explicitly requested.

The other external adapters also use bulk requests: MyGene batches symbols, NCBI uses
attributed ESearch/ESummary batches, UniProt uses paginated exact-symbol OR queries,
BridgeDb uses `xrefsBatch`, and g:Profiler splits large conversion requests into
bounded batches. Transient HTTP failures and rate limits are retried with backoff.

Generated multi-service reports:

- `external_benchmark_summary.md`: combined and per-provider readable tables
- `external_benchmark_summary.csv`: combined external-union metrics
- `external_provider_summary.csv`: metrics for each provider and namespace
- `external_benchmark_details.csv`: complete per-symbol differences
- `external_mapping_support.csv`: each GeneMapKit mapping and supporting services
- `external_benchmark.json`: complete metadata and results

A mapping receives `consensus` support when at least two external services return it,
`single_external_support` when one returns it, and `unsupported_externally` otherwise.
Unsupported service/namespace combinations are not counted as disagreement.

Combined and per-provider precision and recall include Wilson 95% confidence
intervals. The JSON and CSV reports retain the interval method, confidence level,
numerator, and denominator. These mapping-level intervals treat mappings as
observations; because mappings from one input gene can be correlated, a
gene-clustered bootstrap is a useful additional sensitivity analysis for a paper.

External services integrate overlapping upstream databases. Agreement with them is
external concordance and corroboration, not independent accuracy.

## Complete publication benchmark

The multi-service command measures external concordance only. For the complete
publication-oriented evaluation, run:

```bash
python benchmark_publication.py data/sample/all_human_genes.csv \
  --column symbol \
  --db data/build/genemapkit_final_v21.db \
  --gtf data/raw/Homo_sapiens.GRCh38.116.gtf.gz \
  --held-out-report results/mane_heldout_complete \
  --held-out-report results/ncbi_heldout_complete \
  --policy all-supported \
  --path-pivot-type entrez_id \
  --path-sample 0 \
  --path-seed 20260616 \
  --repeats 1 \
  --coordinate-batch-size 1000 \
  --out results/publication_benchmark_all_genes
```

The suite writes:

- `publication_benchmark_summary.md`: combined readable report
- `system_comparison.csv`: GeneMapKit local versus live API timing and offline comparison
- `performance_scaling.csv`: batch runtime, throughput, process RSS, and Python memory
- `coordinate_concordance.csv`: pinned local GTF versus live Ensembl coordinates
- `external_concordance_summary.csv`: multi-service identifier concordance
- `external_provider_timings.csv`: live or cached API timing metadata
- `held_out_validation.csv`: independent accuracy rows from held-out databases
- `path_conversion_coverage.csv`: complete direct all-to-all conversion success
- `path_round_trip_consistency.csv`: original-ID recovery and ambiguity expansion
- `path_pivot_consistency.csv`: direct versus multi-step conversion agreement
- `publication_benchmark.json`: complete results and environment metadata

Coordinate queries are not currently available from the final SQLite database.
Schema 2.1.0 does not store coordinates, so the suite transparently measures coordinate
concordance from the pinned local GTF. Adding an indexed coordinate table is a separate
software feature and database-schema change.

## All-to-all path consistency

The path benchmark starts from every selected namespace and measures:

1. Direct conversion success to every other namespace.
2. Round-trip recovery for `A -> B -> A`.
3. Exact round trips versus legitimate one-to-many ambiguity expansion.
4. Agreement between direct `A -> C` conversion and `A -> Entrez -> C`.

Do not require every round trip to return exactly one identifier. Transcript, protein,
and gene relationships are often one-to-many. The stronger correctness requirement is
that the original identifier remains recoverable and that expansions are explicitly
reported.
