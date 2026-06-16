<div class="hero" markdown>

# GeneMapKit

**Reproducible, provenance-aware conversion between human gene, transcript, protein,
annotation, and structural identifier namespaces.**

GeneMapKit downloads authoritative public source snapshots, builds a normalized SQLite
database, converts identifiers using direct source-asserted relationships, and measures
accuracy against held-out truth sources.

[Get started](getting-started.md){ .md-button .md-button--primary }
[View on GitHub](https://github.com/muhammadmuneeb007/GeneMapKit){ .md-button }

</div>

<div class="metric-grid" markdown>
<div class="metric"><strong>18</strong>identifier namespaces</div>
<div class="metric"><strong>14</strong>supported data sources</div>
<div class="metric"><strong>8.5M</strong>mapping rows in the full v2.1 build</div>
<div class="metric"><strong>29M</strong>direct relationship rows</div>
</div>

## Why GeneMapKit?

Identifier conversion is not always one-to-one. A gene can have multiple transcripts,
a transcript may have multiple cross-database representations, and different sources
can disagree. GeneMapKit makes these cases visible rather than silently selecting the
first match.

- **Reproducible snapshots:** raw files remain unchanged and are recorded with URLs,
  versions, sizes, and SHA-256 checksums.
- **Direct relationships:** conversion queries source-asserted links instead of
  inferring every possible relationship from a broad gene group.
- **Explicit policies:** choose all supported mappings, MANE Select mappings, or an
  authoritative-source preference.
- **Provenance-rich output:** long-format results explain the source, evidence, policy,
  and reason for each returned mapping.
- **Independent validation:** build held-out databases and quantify precision, recall,
  exact-set agreement, over-mapping, and under-mapping.

## Workflow

```text
Public databases
      |
      v
Stage 0: download_databases.py  ->  data/raw + manifest.json
      |
      v
Stage 1: build_database.py      ->  normalized SQLite + build report
      |
      +--> Stage 2: convert identifiers with explicit policies
      |
      +--> Stage 3: validate against held-out relationship truth
```

## Supported identifiers

GeneMapKit supports approved, alias, and previous symbols; HGNC and Entrez IDs;
Ensembl gene, transcript, and protein IDs; RefSeq RNA and protein accessions; UniProt,
CCDS, OMIM, UCSC, Vega, ENA, PDB, and GO identifiers.

## Current scope

GeneMapKit currently targets **Homo sapiens**, GRCh38, with a pinned Ensembl release.
The database is built locally because the complete research snapshot is large and
because reproducibility requires preserving the exact source manifest.

!!! note "Database artifacts"
    Raw snapshots and built SQLite databases are generated artifacts and are not stored
    in Git. Archive publication snapshots separately with their manifests, build
    reports, checksums, software release, and DOI.
