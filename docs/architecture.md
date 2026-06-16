# Architecture

GeneMapKit separates data acquisition, database construction, conversion, and
validation so each scientific claim can be reproduced.

## Pipeline stages

| Stage | Entry point | Output |
|---|---|---|
| Stage 0 | `download_databases.py` | Immutable raw files and `manifest.json` |
| Stage 1 | `build_database.py` | Normalized SQLite database and build report |
| Stage 2 | `python -m genemapkit.cli convert...` | Identifier mappings with provenance |
| Stage 3 | `validate_database.py` | Accuracy metrics and disagreement reports |

## Direct relationships

Broad gene groups are useful for organizing identifiers, but they are unsafe for
conversion when treated as fully transitive. If two transcripts belong to the same
gene, that does not mean each transcript encodes every protein associated with that
gene.

Schema 2.1 therefore stores direct source assertions in `relationships`. Conversion
queries this table instead of expanding every identifier in a group.

## SQLite tables

### `mappings`

Stores normalized and original identifiers with source provenance:

```text
group_id, id_type, id_value, original_value, source, evidence
```

### `relationships`

Stores direct source-asserted links:

```text
source_type, source_value, target_type, target_value,
relationship_type, source, evidence
```

### `metadata`

Records schema version, snapshot details, selected sources, and build completeness.

### `source_summary`

Records parsed records, records without anchors, mapping rows, and relationship rows
for each source.

## Identifier normalization

Versioned Ensembl, RefSeq, and CCDS accessions are normalized for lookup while their
original values remain in the database. Symbols and aliases are matched exactly as
tokens; substring alias matching is not used.

## Provenance

Each long-format conversion result includes:

- Requested and matched input namespace
- Normalized input
- Output identifier and namespace
- Mapping count and status
- Mapping policy and selection reason
- Source and source-level evidence

## Current limitations

- Human identifiers only
- Reviewed versus unreviewed UniProt status is not yet represented
- Authority precedence is explicit but still requires independent evaluation
- Conflicting authoritative assertions are not yet represented as a separate status
- Built database redistribution depends on upstream source license review
