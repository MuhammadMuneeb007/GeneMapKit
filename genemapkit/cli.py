#!/usr/bin/env python3
"""
GeneMapKit - Comprehensive Gene ID Mapping Toolkit

A powerful command-line tool for converting between different gene identifier formats
using data from multiple authoritative sources including HGNC, NCBI, Ensembl, and UniProt.

Author: GeneMapKit Development Team
License: MIT
"""

import click
import sys
import logging
from pathlib import Path
import pandas as pd

from genemapkit.core.downloader import (
    DEFAULT_ENSEMBL_RELEASE,
    TIERS,
    VALIDATION_SOURCES,
    DatabaseDownloader,
)
from genemapkit.core.builder import PARSERS, build_database, convert_identifier
from genemapkit.core.converter import (
    MAPPING_POLICIES,
    SUPPORTED_ID_TYPES,
    SYMBOL_POLICIES,
    GeneConverter,
)
from genemapkit.core.validation import VALIDATION_CONVERSIONS, validate_database
from genemapkit.core.external_validation import MYGENE_OUTPUT_FIELDS, compare_mygene
from genemapkit.core.external_benchmark import (
    DEFAULT_OUTPUTS,
    DEFAULT_PROVIDERS,
    PROVIDER_OUTPUTS,
    benchmark_external_services,
)
from genemapkit.core.publication_benchmark import run_publication_benchmark
from genemapkit.core.path_benchmark import DEFAULT_PATH_ID_TYPES, benchmark_conversion_paths
from genemapkit.utils.logger import setup_logging
from genemapkit.utils.validators import validate_input_file


# Supported ID types
@click.group()
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose logging')
@click.option('--quiet', '-q', is_flag=True, help='Suppress all output except errors')
@click.pass_context
def cli(ctx, verbose, quiet):
    """GeneMapKit - Gene ID Mapping Toolkit"""
    ctx.ensure_object(dict)
    
    # Setup logging
    log_level = 'DEBUG' if verbose else 'WARNING' if quiet else 'INFO'
    ctx.obj['logger'] = setup_logging(log_level)
    

@cli.command()
@click.option('--force', '-f', is_flag=True, help='Force re-download even if files exist')
@click.option('--sources', '-s', multiple=True,
              type=click.Choice(list(DatabaseDownloader.SOURCES) + ['core', 'default', 'all']),
              help='Explicit sources or source tier to download')
@click.option('--core-only', is_flag=True, help='Download only essential core sources')
@click.option('--all', 'all_sources', is_flag=True,
              help='Download every registered source, including large/niche sources')
@click.option('--validation', is_flag=True,
              help='Also download independent benchmark-validation sources')
@click.option('--out-dir', '--raw-dir', '--out', default='data/raw', show_default=True,
              help='Raw snapshot directory')
@click.option('--ensembl-release', default=DEFAULT_ENSEMBL_RELEASE, show_default=True,
              type=int, help='Pinned Ensembl release')
@click.pass_context
def download(ctx, force, sources, core_only, all_sources, validation,
             out_dir, ensembl_release):
    """Download immutable raw database snapshots with a provenance manifest."""
    logger = ctx.obj['logger']
    logger.info("Starting database download process...")

    downloader = DatabaseDownloader(
        data_dir=out_dir, logger=logger, ensembl_release=ensembl_release
    )
    
    try:
        if core_only and (all_sources or validation or sources):
            raise click.UsageError(
                "--core-only cannot be combined with --all, --validation, or --sources"
            )

        # Handle special cases
        if all_sources or 'all' in sources:
            results = downloader.download_all(force_download=force, include_optional=True)
        elif core_only or 'core' in sources:
            results = downloader.download_specific_sources(
                downloader.sources_for_tiers(('core',)), force_download=force
            )
        elif validation:
            if sources:
                source_list = []
                if 'core' in sources:
                    source_list.extend(downloader.sources_for_tiers(('core',)))
                if 'default' in sources:
                    source_list.extend(
                        downloader.sources_for_tiers(('core', 'crossref'))
                    )
                source_list.extend(
                    s for s in sources if s not in ['core', 'default', 'all']
                )
                source_list.extend(downloader.validation_sources())
                results = downloader.download_specific_sources(
                    source_list, force_download=force
                )
            else:
                results = downloader.download_validation(
                    force_download=force, include_default=True
                )
        elif 'default' in sources:
            results = downloader.download_default(force_download=force)
        elif sources:
            source_list = [s for s in sources if s not in ['core', 'default', 'all']]
            results = downloader.download_specific_sources(source_list, force_download=force)
        else:
            results = downloader.download_default(force_download=force)
        
        # Display download summary
        successful = sum(1 for success in results.values() if success)
        total = len(results)
        
        if successful == total:
            logger.info("🎉 All database downloads completed successfully!")
        else:
            logger.warning(f"⚠️  {successful}/{total} downloads completed successfully")
            failed_sources = [source for source, success in results.items() if not success]
            logger.warning(f"Failed sources: {failed_sources}")
        
    except click.ClickException:
        raise
    except Exception as error:
        raise click.ClickException(str(error)) from error


@cli.command(name='download-for')
@click.option('--input-type', '-i', required=True,
              type=click.Choice(SUPPORTED_ID_TYPES),
              help='Type of input gene identifiers')
@click.option('--output-type', '-o', multiple=True, required=True,
              type=click.Choice(SUPPORTED_ID_TYPES),
              help='Type(s) of output gene identifiers')
@click.option('--minimal', is_flag=True,
              help='Download only minimal required databases')
@click.option('--force', '-f', is_flag=True,
              help='Force re-download even if files exist')
@click.option('--out-dir', '--raw-dir', '--out', default='data/raw', show_default=True,
              help='Raw snapshot directory')
@click.option('--ensembl-release', default=DEFAULT_ENSEMBL_RELEASE, show_default=True,
              type=int, help='Pinned Ensembl release')
@click.pass_context
def download_for_conversion(ctx, input_type, output_type, minimal, force,
                            out_dir, ensembl_release):
    """Download databases specifically for a gene ID conversion."""
    logger = ctx.obj['logger']
    logger.info(f"Downloading databases for conversion: {input_type} -> {output_type}")
    
    try:
        downloader = DatabaseDownloader(
            data_dir=out_dir, logger=logger, ensembl_release=ensembl_release
        )
        results = downloader.download_for_conversion(
            input_type, list(output_type), force_download=force, minimal=minimal
        )
        
        # Display download summary
        successful = sum(1 for success in results.values() if success)
        total = len(results)
        
        if successful == total:
            logger.info("🎉 All required databases downloaded successfully!")
            logger.info("You can now perform your gene ID conversion.")
        else:
            logger.warning(f"⚠️  {successful}/{total} databases downloaded successfully")
            failed_sources = [source for source, success in results.items() if not success]
            logger.warning(f"Failed sources: {failed_sources}")
            
    except click.ClickException:
        raise
    except Exception as error:
        raise click.ClickException(str(error)) from error


@cli.command(name='build')
@click.option('--raw-dir', default='data/raw', show_default=True,
              type=click.Path(path_type=Path),
              help='Snapshot directory containing raw files and manifest.json')
@click.option('--out', default='data/build/genemapkit.db', show_default=True,
              type=click.Path(path_type=Path), help='Output SQLite database')
@click.option('--report', default='data/build/build_report.json', show_default=True,
              type=click.Path(path_type=Path), help='Output JSON build report')
@click.option('--source', 'sources', multiple=True, type=click.Choice(sorted(PARSERS)),
              help='Build from this source; repeat to select multiple sources')
@click.option('--exclude', multiple=True, type=click.Choice(sorted(PARSERS)),
              help='Exclude this source; repeat to exclude multiple sources')
@click.option('--no-verify', is_flag=True,
              help='Skip manifest checksum verification; not for publication builds')
@click.option('--allow-missing', is_flag=True,
              help='Continue when selected sources are missing or empty')
@click.option('--force', is_flag=True, help='Replace an existing output database')
@click.option('--max-records', type=click.IntRange(min=1),
              help='Development-only record limit per source; creates an incomplete index')
@click.option('--check', nargs=3, metavar='VALUE INPUT_TYPE OUTPUT_TYPE',
              help='Run one conversion check after building')
@click.option('--progress/--no-progress', default=True, show_default=True,
              help='Show per-source progress bars during the build')
@click.pass_context
def build_index(ctx, raw_dir, out, report, sources, exclude, no_verify,
                allow_missing, force, max_records, check, progress):
    """Build the normalized SQLite identifier index from a raw snapshot."""
    logger = ctx.obj['logger']
    try:
        result = build_database(
            raw_dir,
            out,
            sources=list(sources) or None,
            exclude=list(exclude) or None,
            verify=not no_verify,
            strict=not allow_missing,
            force=force,
            max_records=max_records,
            report_path=report,
            progress=progress,
        )
        logger.info(
            "Build complete: %s groups, %s mapping rows, and %s direct relationships",
            result['groups'],
            result['mapping_rows'],
            result.get('relationship_rows', 0),
        )
        logger.info("Database: %s", out)
        logger.info("Report: %s", report)
        if check:
            value, input_type, output_type = check
            results = convert_identifier(out, value, input_type, output_type)
            click.echo(
                f"{value} ({input_type}) -> {output_type}: "
                + (", ".join(results) if results else "(no mappings)")
            )
    except Exception as error:
        raise click.ClickException(str(error)) from error


@cli.command(name='validate')
@click.option('--db', default='data/build/genemapkit.db', show_default=True,
              type=click.Path(path_type=Path), help='Built Stage 1 SQLite database')
@click.option('--raw-dir', default='data/raw', show_default=True,
              type=click.Path(path_type=Path), help='Raw snapshot with truth sources')
@click.option('--out', default='results/validation', show_default=True,
              type=click.Path(path_type=Path), help='Validation report directory')
@click.option('--source', 'sources', multiple=True,
              type=click.Choice(sorted(VALIDATION_CONVERSIONS)),
              help='Truth source; repeat to select multiple')
@click.option('--sample', default=1000, show_default=True,
              type=click.IntRange(min=0),
              help='Unique truth inputs per conversion; 0 evaluates complete sources')
@click.option('--max-examples', default=25, show_default=True,
              type=click.IntRange(min=0), help='Disagreement examples per conversion')
@click.option('--require-held-out', is_flag=True,
              help='Skip truth sources included in the database build')
@click.option('--policy', 'mapping_policy', default='all-supported', show_default=True,
              type=click.Choice(MAPPING_POLICIES),
              help='Mapping-selection policy to benchmark')
@click.option('--progress/--no-progress', default=True, show_default=True,
              help='Show progress bars for long-running validation loops')
@click.pass_context
def validate_index(ctx, db, raw_dir, out, sources, sample, max_examples,
                   require_held_out, mapping_policy, progress):
    """Benchmark conversion sets against MANE and held-out NCBI relationships."""
    logger = ctx.obj['logger']
    try:
        report = validate_database(
            db,
            raw_dir,
            sources=list(sources) or None,
            sample=None if sample == 0 else sample,
            max_examples=max_examples,
            require_held_out=require_held_out,
            mapping_policy=mapping_policy,
            output_dir=out,
            progress=progress,
        )
        logger.info("Validation reports: %s", out)
        for row in report['benchmarks']:
            logger.info(
                "%s %s -> %s: precision=%s%% recall=%s%% overmapped=%s/%s independent=%s",
                row['truth_source'],
                row['input_type'],
                row['output_type'],
                row['precision_pct'],
                row['recall_pct'],
                row['overmapped_inputs'],
                row['tested_inputs'],
                row['independent'],
            )
        for skipped in report['skipped']:
            logger.warning("Skipped %s: %s", skipped['source'], skipped['reason'])
    except Exception as error:
        raise click.ClickException(str(error)) from error


@cli.command(name='compare-mygene')
@click.argument('input_file', type=click.Path(exists=True, path_type=Path))
@click.option('--column', default='symbol', show_default=True,
              help='Column containing human gene symbols')
@click.option('--delimiter', default=',', show_default=True,
              help=r"Input delimiter; use \t for TSV")
@click.option('--db', default='data/build/genemapkit_final_v21.db', show_default=True,
              type=click.Path(path_type=Path), help='Final built SQLite database')
@click.option('--output-type', 'output_types', multiple=True,
              type=click.Choice(sorted(MYGENE_OUTPUT_FIELDS)),
              help='Output namespace to compare; repeatable')
@click.option('--policy', 'mapping_policy', default='all-supported', show_default=True,
              type=click.Choice(MAPPING_POLICIES), help='GeneMapKit mapping policy')
@click.option('--out', default='results/mygene_concordance', show_default=True,
              type=click.Path(path_type=Path), help='Output report directory')
@click.option('--progress/--no-progress', default=True, show_default=True)
def compare_mygene_command(input_file, column, delimiter, db, output_types,
                           mapping_policy, out, progress):
    """Compare symbol conversions with MyGene.info external annotations."""
    try:
        delimiter = '\t' if delimiter == r'\t' else delimiter
        frame = pd.read_csv(input_file, delimiter=delimiter, dtype=str)
        if column not in frame:
            raise click.UsageError(
                f"Column '{column}' not found. Available columns: "
                f"{', '.join(frame.columns)}"
            )
        outputs = list(output_types) or [
            'ensembl_gene_id',
            'entrez_id',
            'refseq_rna_id',
            'refseq_protein_id',
            'uniprot_id',
        ]
        report = compare_mygene(
            db,
            frame[column].dropna(),
            outputs,
            mapping_policy=mapping_policy,
            output_dir=out,
            progress=progress,
        )
        click.echo(
            "MyGene.info is an external integrated annotation service; "
            "these results measure concordance, not independent accuracy."
        )
        for row in report['summaries']:
            click.echo(
                f"{row['output_type']}: exact={row['exact_set_inputs']}/"
                f"{row['tested_symbols']} precision={row['concordance_precision_pct']}% "
                f"recall={row['concordance_recall_pct']}%"
            )
        click.echo(f"Reports: {out}")
    except click.ClickException:
        raise
    except Exception as error:
        raise click.ClickException(str(error)) from error


@cli.command(name='benchmark-external')
@click.argument('input_file', type=click.Path(exists=True, path_type=Path))
@click.option('--column', default='symbol', show_default=True,
              help='Column containing human gene symbols')
@click.option('--delimiter', default=',', show_default=True,
              help=r"Input delimiter; use \t for TSV")
@click.option('--db', default='data/build/genemapkit_final_v21.db', show_default=True,
              type=click.Path(path_type=Path), help='Final built SQLite database')
@click.option('--provider', 'providers', multiple=True,
              type=click.Choice(sorted(PROVIDER_OUTPUTS)),
              help='External service; repeatable (default: all)')
@click.option('--output-type', 'output_types', multiple=True,
              type=click.Choice(sorted(set().union(*PROVIDER_OUTPUTS.values()))),
              help='Output namespace; repeatable')
@click.option('--policy', 'mapping_policy', default='all-supported', show_default=True,
              type=click.Choice(MAPPING_POLICIES), help='GeneMapKit mapping policy')
@click.option('--out', default='results/external_benchmark', show_default=True,
              type=click.Path(path_type=Path), help='Output report and cache directory')
@click.option('--refresh-cache', is_flag=True, help='Ignore cached API responses')
@click.option('--progress/--no-progress', default=True, show_default=True)
def benchmark_external_command(input_file, column, delimiter, db, providers,
                               output_types, mapping_policy, out, refresh_cache, progress):
    """Compare symbol mappings with six supported external services."""
    try:
        delimiter = '\t' if delimiter == r'\t' else delimiter
        frame = pd.read_csv(input_file, delimiter=delimiter, dtype=str)
        if column not in frame:
            raise click.UsageError(
                f"Column '{column}' not found. Available columns: "
                f"{', '.join(frame.columns)}"
            )
        report = benchmark_external_services(
            db,
            frame[column].dropna(),
            list(output_types) or DEFAULT_OUTPUTS,
            providers=list(providers) or DEFAULT_PROVIDERS,
            mapping_policy=mapping_policy,
            output_dir=out,
            refresh_cache=refresh_cache,
            progress=progress,
        )
        for row in report['summaries']:
            click.echo(
                f"{row['output_type']}: precision={row['concordance_precision_pct']}% "
                f"[{row['concordance_precision_ci_low_pct']}, "
                f"{row['concordance_precision_ci_high_pct']}] "
                f"recall={row['concordance_recall_pct']}% "
                f"[{row['concordance_recall_ci_low_pct']}, "
                f"{row['concordance_recall_ci_high_pct']}] "
                f"conflicts={row['conflict']} GMK-only={row['genemapkit_only']}"
            )
        for provider, error in report['provider_errors'].items():
            click.echo(f"Provider error - {provider}: {error}", err=True)
        click.echo(
            "External API agreement measures concordance, not independent accuracy."
        )
        click.echo(f"Reports and API cache: {out}")
    except click.ClickException:
        raise
    except Exception as error:
        raise click.ClickException(str(error)) from error


@cli.command(name='benchmark-publication')
@click.argument('input_file', type=click.Path(exists=True, path_type=Path))
@click.option('--column', default='symbol', show_default=True)
@click.option('--delimiter', default=',', show_default=True)
@click.option('--db', default='data/build/genemapkit_final_v21.db', show_default=True,
              type=click.Path(exists=True, path_type=Path))
@click.option('--gtf', default='data/raw/Homo_sapiens.GRCh38.116.gtf.gz',
              show_default=True, type=click.Path(exists=True, path_type=Path))
@click.option('--held-out-report', 'held_out_reports', multiple=True,
              type=click.Path(path_type=Path),
              help='Directory containing validation_metrics.json; repeatable')
@click.option('--batch-size', 'batch_sizes', multiple=True, type=click.IntRange(min=1),
              help='Local scalability batch size; repeatable')
@click.option('--repeats', default=3, show_default=True, type=click.IntRange(min=1))
@click.option('--policy', 'mapping_policy', default='all-supported', show_default=True,
              type=click.Choice(MAPPING_POLICIES))
@click.option('--out', default='results/publication_benchmark', show_default=True,
              type=click.Path(path_type=Path))
@click.option('--refresh-external-cache', is_flag=True)
@click.option('--coordinate-batch-size', default=1000, show_default=True,
              type=click.IntRange(min=1),
              help='Symbols per Ensembl coordinate API request')
@click.option('--path-pivot-type', default='entrez_id', show_default=True,
              type=click.Choice(DEFAULT_PATH_ID_TYPES))
@click.option('--path-sample', default=1000, show_default=True,
              type=click.IntRange(min=0),
              help='Identifiers per namespace for path tests; 0 means complete')
@click.option('--path-seed', default=20260616, show_default=True, type=int)
@click.option('--progress/--no-progress', default=True, show_default=True)
def benchmark_publication_command(input_file, column, delimiter, db, gtf,
                                  held_out_reports, batch_sizes, repeats,
                                  mapping_policy, out, refresh_external_cache,
                                  coordinate_batch_size, path_pivot_type,
                                  path_sample, path_seed, progress):
    """Run performance, coordinates, APIs, and held-out validation together."""
    try:
        delimiter = '\t' if delimiter == r'\t' else delimiter
        frame = pd.read_csv(input_file, delimiter=delimiter, dtype=str)
        if column not in frame:
            raise click.UsageError(
                f"Column '{column}' not found. Available columns: "
                f"{', '.join(frame.columns)}"
            )
        reports = list(held_out_reports) or [
            Path('results/mane_heldout_all_supported'),
            Path('results/ncbi_heldout_all_supported'),
        ]
        report = run_publication_benchmark(
            db,
            gtf,
            frame[column].dropna(),
            output_types=DEFAULT_OUTPUTS,
            providers=DEFAULT_PROVIDERS,
            held_out_report_dirs=reports,
            batch_sizes=list(batch_sizes) or [1, 10, 100, 1000],
            repeats=repeats,
            mapping_policy=mapping_policy,
            output_dir=out,
            refresh_external_cache=refresh_external_cache,
            coordinate_batch_size=coordinate_batch_size,
            path_pivot_type=path_pivot_type,
            path_sample=None if path_sample == 0 else path_sample,
            path_random_seed=path_seed,
            progress=progress,
        )
        click.echo(f"Publication benchmark: {out / Path('publication_benchmark_summary.md')}")
        click.echo(
            f"Offline local conversion: "
            f"{report['local_performance']['offline_operation_passed']}"
        )
        click.echo(
            f"Independent held-out rows: "
            f"{len(report['held_out_validation']['independent_benchmarks'])}"
        )
    except click.ClickException:
        raise
    except Exception as error:
        raise click.ClickException(str(error)) from error


@cli.command(name='benchmark-paths')
@click.option('--db', default='data/build/genemapkit_final_v21.db', show_default=True,
              type=click.Path(exists=True, path_type=Path))
@click.option('--id-type', 'id_types', multiple=True,
              type=click.Choice(DEFAULT_PATH_ID_TYPES),
              help='Namespace to include; repeatable')
@click.option('--pivot-type', default='entrez_id', show_default=True,
              type=click.Choice(DEFAULT_PATH_ID_TYPES))
@click.option('--sample', default=1000, show_default=True,
              type=click.IntRange(min=0), help='Path sample per namespace; 0 means complete')
@click.option('--seed', default=20260616, show_default=True, type=int)
@click.option('--max-examples', default=100, show_default=True,
              type=click.IntRange(min=0))
@click.option('--out', default='results/path_benchmark', show_default=True,
              type=click.Path(path_type=Path))
@click.option('--progress/--no-progress', default=True, show_default=True)
def benchmark_paths_command(db, id_types, pivot_type, sample, seed, max_examples, out,
                            progress):
    """Measure all-to-all coverage, round trips, and pivot path consistency."""
    try:
        report = benchmark_conversion_paths(
            db,
            id_types=list(id_types) or DEFAULT_PATH_ID_TYPES,
            pivot_type=pivot_type,
            path_sample=None if sample == 0 else sample,
            random_seed=seed,
            max_examples=max_examples,
            output_dir=out,
            progress=progress,
        )
        click.echo(f"Runtime: {report['elapsed_seconds']} seconds")
        click.echo(f"Pivot namespace: {report['pivot_type']}")
        click.echo(f"Reports: {out}")
    except Exception as error:
        raise click.ClickException(str(error)) from error


@cli.command(name='convert-one')
@click.argument('value')
@click.option('--input-type', '--from', '-i', required=True,
              type=click.Choice(SUPPORTED_ID_TYPES),
              help='Input identifier namespace')
@click.option('--output-type', '--to', '-o', multiple=True, required=True,
              type=click.Choice(SUPPORTED_ID_TYPES),
              help='Target namespace; repeat for multiple targets')
@click.option('--db', default='data/build/genemapkit.db', show_default=True,
              type=click.Path(path_type=Path), help='Built Stage 1 SQLite database')
@click.option('--symbol-policy', default='all', show_default=True,
              type=click.Choice(sorted(SYMBOL_POLICIES)),
              help='Which symbol categories are searched for symbol inputs')
@click.option('--policy', 'mapping_policy', default='all-supported', show_default=True,
              type=click.Choice(MAPPING_POLICIES),
              help='Mapping-selection policy')
@click.option('--include-retired', is_flag=True,
              help='Include retired target identifiers in results')
@click.option('--provenance/--no-provenance', default=True, show_default=True,
              help='Show source and evidence supporting each mapping')
def convert_one(value, input_type, output_type, db, symbol_policy, mapping_policy,
                include_retired, provenance):
    """Convert one identifier using the built SQLite database."""
    try:
        converter = GeneConverter(db, symbol_policy=symbol_policy)
        found = False
        for target in output_type:
            results = converter.query(
                value,
                input_type,
                target,
                mapping_policy=mapping_policy,
                include_retired=include_retired,
            )
            if not results:
                click.echo(f"{target}: (unmapped)")
                continue
            found = True
            seen = set()
            for result in results:
                key = (
                    (
                        result.output_id,
                        result.matched_input_type,
                        result.output_source,
                        result.output_evidence,
                    )
                    if provenance else result.output_id
                )
                if key in seen:
                    continue
                seen.add(key)
                line = (
                    f"{target}: {result.output_id} "
                    f"[{result.mapping_status}; policy={result.mapping_policy}; "
                    f"matched={result.matched_input_type}]"
                )
                if provenance:
                    line += (
                        f" source={result.output_source} "
                        f"evidence={result.output_evidence} "
                        f"reason={result.selection_reason}"
                    )
                click.echo(line)
        if not found:
            raise click.exceptions.Exit(2)
    except click.exceptions.Exit:
        raise
    except Exception as error:
        raise click.ClickException(str(error)) from error


@cli.command()
@click.argument('input_file', type=click.Path(exists=True, path_type=Path))
@click.option('--input-type', '-i', required=True,
              type=click.Choice(SUPPORTED_ID_TYPES),
              help='Input identifier namespace')
@click.option('--output-type', '-o', multiple=True, required=True,
              type=click.Choice(SUPPORTED_ID_TYPES),
              help='Target namespace; repeat for multiple targets')
@click.option('--output', '-out', default='output_mapped.csv', show_default=True,
              type=click.Path(path_type=Path), help='Output CSV path')
@click.option('--column', '-c', required=True,
              help='Column containing input identifiers')
@click.option('--delimiter', '-d', default=',', show_default=True,
              help=r"Input delimiter; use \t for TSV")
@click.option('--db', default='data/build/genemapkit.db', show_default=True,
              type=click.Path(path_type=Path), help='Built Stage 1 SQLite database')
@click.option('--format', 'output_format', default='wide', show_default=True,
              type=click.Choice(['wide', 'long']),
              help='Wide joins multiple IDs; long preserves mapping provenance')
@click.option('--separator', default=';', show_default=True,
              help='Separator for multiple IDs in wide output')
@click.option('--symbol-policy', default='all', show_default=True,
              type=click.Choice(sorted(SYMBOL_POLICIES)),
              help='Which symbol categories are searched for symbol inputs')
@click.option('--policy', 'mapping_policy', default='all-supported', show_default=True,
              type=click.Choice(MAPPING_POLICIES),
              help='Mapping-selection policy')
@click.option('--include-retired', is_flag=True,
              help='Include retired target identifiers in results')
@click.option('--include-unmapped/--drop-unmapped', default=True, show_default=True,
              help='Keep or remove unmapped inputs')
@click.option('--progress/--no-progress', default=True, show_default=True,
              help='Show identifier conversion progress')
@click.pass_context
def convert(ctx, input_file, input_type, output_type, output, column, delimiter, db,
            output_format, separator, symbol_policy, mapping_policy, include_retired,
            include_unmapped, progress):
    """Convert a CSV/TSV column using the built Stage 1 SQLite database."""
    logger = ctx.obj['logger']
    try:
        validate_input_file(str(input_file), logger)
        delimiter = '\t' if delimiter == r'\t' else delimiter
        dataframe = pd.read_csv(input_file, delimiter=delimiter, dtype=str)
        if column not in dataframe.columns:
            raise click.UsageError(
                f"Column '{column}' not found. Available columns: "
                f"{', '.join(dataframe.columns)}"
            )

        converter = GeneConverter(db, symbol_policy=symbol_policy)
        converted = converter.convert_batch(
            dataframe[column],
            input_type,
            output_type,
            include_unmapped=include_unmapped,
            output_format=output_format,
            separator=separator,
            mapping_policy=mapping_policy,
            include_retired=include_retired,
            progress=progress,
        )
        original = dataframe.reset_index().rename(columns={'index': 'input_row'})
        result = original.merge(converted, on='input_row', how='inner')
        output.parent.mkdir(parents=True, exist_ok=True)
        result.drop(columns=['input_row']).to_csv(output, index=False)

        if output_format == 'long':
            mapped = (
                converted.loc[
                    converted['mapping_status'] != 'unmapped', 'input_row'
                ].nunique()
                if 'mapping_status' in converted else 0
            )
        else:
            status_columns = [
                column_name for column_name in converted
                if column_name.endswith('_mapping_status')
            ]
            mapped = (
                converted[status_columns].ne('unmapped').any(axis=1).sum()
                if status_columns else 0
            )
        logger.info("Converted %d input rows; %d had at least one mapping", len(dataframe), mapped)
        logger.info("Output saved to: %s", output)
        if output_format == 'wide':
            logger.warning(
                "Wide output joins multiple mappings with '%s'; use --format long "
                "for full provenance.",
                separator,
            )
    except click.ClickException:
        raise
    except Exception as error:
        raise click.ClickException(str(error)) from error


def _print_portable_info() -> None:
    """Print current identifier and source information using portable ASCII."""
    click.echo("\nGeneMapKit - Gene ID Mapping Toolkit")
    click.echo("=" * 60)
    click.echo("\nSupported identifier types:")
    for id_type in SUPPORTED_ID_TYPES:
        click.echo(f"  - {id_type}")
    click.echo("\nRegistered database sources:")
    for tier in TIERS:
        click.echo(f"  {tier}:")
        sources = DatabaseDownloader(data_dir="data/raw").sources_for_tiers((tier,))
        for source in sources:
            suffix = " [validation]" if source in VALIDATION_SOURCES else ""
            click.echo(f"    - {source}{suffix}")
    click.echo("\nQuick start:")
    click.echo("  genemapkit download --all --validation")
    click.echo("  genemapkit status --raw-dir data/raw")
    click.echo("  genemapkit build --raw-dir data/raw --force")


def _print_portable_status(raw_dir: Path) -> None:
    """Print all registered source statuses using portable ASCII."""
    downloader = DatabaseDownloader(data_dir=str(raw_dir))
    all_sources = []
    click.echo("\nDatabase Status")
    click.echo("=" * 60)
    for tier in TIERS:
        sources = downloader.sources_for_tiers((tier,))
        all_sources.extend(sources)
        click.echo(f"\n{tier.upper()} sources:")
        for source in sources:
            file_info = downloader.get_file_info(source)
            validation = " [validation]" if source in VALIDATION_SOURCES else ""
            if file_info:
                click.echo(
                    f"  available  {source}{validation}: "
                    f"{file_info['filename']} ({file_info['size_mb']:.1f} MB)"
                )
            else:
                click.echo(f"  missing    {source}{validation}")
    available = sum(1 for source in all_sources if downloader.get_file_info(source))
    click.echo(f"\nSnapshot directory: {raw_dir}")
    click.echo(f"Summary: {available}/{len(all_sources)} sources available")


@cli.command()
@click.pass_context
def info(ctx):
    """Show information about available databases and supported ID types."""
    _print_portable_info()

@cli.command()
@click.option('--raw-dir', '--out-dir', default='data/raw', show_default=True,
              type=click.Path(path_type=Path),
              help='Raw snapshot directory to inspect')
@click.pass_context
def status(ctx, raw_dir):
    """Check status of downloaded databases."""
    _print_portable_status(raw_dir)


if __name__ == '__main__':
    cli()
