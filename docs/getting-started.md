# Getting Started

This guide creates a complete local GeneMapKit database and runs a first conversion.

## Requirements

- Python 3.10 or newer
- Git
- Sufficient disk space for downloaded sources, temporary SQLite journals, and builds
- A stable internet connection for downloading public databases

The full research snapshot contains multi-gigabyte NCBI files. Keep substantially more
free disk space than the final database size because SQLite uses temporary files while
building.

## Install

=== "Windows PowerShell"

    ```powershell
    git clone https://github.com/muhammadmuneeb007/GeneMapKit.git
    cd GeneMapKit

    python -m venv .venv
    .\.venv\Scripts\Activate.ps1
    python -m pip install .
    genemapkit --help
    ```

=== "Linux or macOS"

    ```bash
    git clone https://github.com/muhammadmuneeb007/GeneMapKit.git
    cd GeneMapKit

    python -m venv .venv
    source .venv/bin/activate
    python -m pip install .
    genemapkit --help
    ```

## Quick production workflow

```bash
# Download the recommended conversion snapshot.
genemapkit download --out data/raw

# Build and verify the SQLite database.
genemapkit build \
  --raw-dir data/raw \
  --out data/build/genemapkit.db \
  --report data/build/build_report.json

# Convert one identifier.
python -m genemapkit.cli convert-one TP53 \
  --from symbol \
  --to ensembl_gene_id \
  --policy authoritative \
  --db data/build/genemapkit.db
```

On Windows Command Prompt, write each command on one line. In PowerShell, use a
backtick instead of `\` for line continuation.

## Full research workflow

For GO, PDB, and validation mappings:

```bash
python download_databases.py --out data/raw --all --validation

python build_database.py \
  --raw-dir data/raw \
  --out data/build/genemapkit-full.db \
  --report data/build/genemapkit-full-report.json
```

## Check the build

The production build report should contain:

```json
{
  "build_schema_version": "2.1.0",
  "missing_sources": [],
  "verified_snapshot": true,
  "complete_build": true
}
```

Never use a database with `complete_build: false` for scientific benchmarking.

## Next

1. Learn the [download modes](download.md).
2. Understand [production and held-out builds](build.md).
3. Choose a [mapping policy](policies.md).
4. Run [independent validation](validate.md).
