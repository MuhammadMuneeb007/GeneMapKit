# Mapping Policies

Identifier conversion can legitimately return multiple results. Policies make result
selection explicit and reproducible.

## Available policies

| Policy | Behavior | Recommended use |
|---|---|---|
| `all-supported` | Returns every direct source-asserted mapping | Complete source-level discovery and provenance |
| `direct-only` | Guarantees only direct relationships | Explicitly prevent future inferred traversal |
| `corroborated` | Returns mappings supported by at least two source databases | Higher-confidence cross-source mappings |
| `mane-select` | Returns only mappings with `MANE Select` evidence | Canonical Ensembl/RefSeq transcript and protein crosswalks |
| `authoritative` | Uses the highest-priority available authority for the target namespace | General-purpose selected output |

`all-supported` is the backward-compatible default.

## Example

```bash
python -m genemapkit.cli convert-one ENST00000269305 \
  --from ensembl_transcript_id \
  --to refseq_rna_id \
  --policy all-supported \
  --db data/build/genemapkit.db
```

For the TP53 transcript, `all-supported` can return several directly supported RefSeq
RNAs. `mane-select` returns the canonical MANE Select mapping when present.

```bash
python -m genemapkit.cli convert-one ENST00000269305 \
  --from ensembl_transcript_id \
  --to refseq_rna_id \
  --policy mane-select \
  --db data/build/genemapkit.db
```

## MANE Select behavior

`mane-select`:

- Requires the source to be `mane_summary`.
- Requires evidence status exactly equal to `MANE Select`.
- Excludes `MANE Plus Clinical`.
- Does not silently fall back to another source.

An empty result means no direct MANE Select assertion was present for that conversion.

## Authoritative behavior

`authoritative` selects the highest-priority source that has a direct mapping for the
requested target namespace. It returns every mapping asserted by that selected source;
it never picks an arbitrary first identifier.

Examples of preferred authorities include:

- HGNC for approved symbols and HGNC identifiers
- NCBI for Entrez and RefSeq identifiers
- Ensembl for Ensembl entities
- UniProt for UniProt identifiers
- SIFTS for PDB mappings
- MANE for canonical Ensembl/RefSeq transcript crosswalks when available

The exact precedence table is defined in `genemapkit/core/converter.py`.

## Policy provenance

Long-format results include:

```text
mapping_policy
selection_reason
source_support_count
supporting_sources
agreement_status
output_source
output_evidence
```

This allows downstream workflows and publications to explain why a mapping was
returned.

## Choosing a policy

- Use `all-supported` when completeness and disagreements matter.
- Use `corroborated` when at least two source databases must support a mapping.
- Use `mane-select` for canonical MANE crosswalks.
- Use `authoritative` for a selected general-purpose result.
- Use `direct-only` when explicitly documenting that inferred traversal is prohibited.

Benchmark policies independently before making accuracy claims.
