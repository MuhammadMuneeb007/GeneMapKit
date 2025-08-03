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
import os
import logging
from pathlib import Path
import pandas as pd

# Add the project root to Python path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from genemapkit.core.downloader import DatabaseDownloader
from genemapkit.core.converter import GeneConverter
from genemapkit.utils.logger import setup_logging
from genemapkit.utils.validators import validate_input_file, validate_id_type


# Supported ID types
SUPPORTED_ID_TYPES = [
    'symbol', 'ensembl_gene_id', 'ensembl_transcript_id', 
    'entrez_id', 'hgnc_id', 'refseq_id', 'uniprot_id'
]

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
              type=click.Choice(['hgnc', 'ncbi_gene_info', 'ensembl_gtf', 
                               'uniprot_idmapping', 'uniprot_human', 'refseq_human', 
                               'ensembl_biomart_sample', 'core', 'all']),
              default=['core'], help='Specify which databases to download')
@click.option('--include-optional', is_flag=True, help='Include optional/supplementary databases')
@click.pass_context
def download(ctx, force, sources, include_optional):
    """Download gene annotation databases from multiple sources."""
    logger = ctx.obj['logger']
    logger.info("Starting database download process...")
    
    downloader = DatabaseDownloader(logger=logger)
    
    try:
        # Handle special cases
        if 'all' in sources:
            results = downloader.download_all(force_download=force, include_optional=True)
        elif 'core' in sources:
            core_sources = ['hgnc', 'ncbi_gene_info', 'ensembl_gtf', 'uniprot_human']
            if include_optional:
                core_sources.extend(['uniprot_idmapping', 'refseq_human', 'ensembl_biomart_sample'])
            results = downloader.download_specific_sources(core_sources, force_download=force)
        else:
            # Download specific sources
            source_list = [s for s in sources if s not in ['core', 'all']]
            results = downloader.download_specific_sources(source_list, force_download=force)
        
        # Display download summary
        successful = sum(1 for success in results.values() if success)
        total = len(results)
        
        if successful == total:
            logger.info("🎉 All database downloads completed successfully!")
        else:
            logger.warning(f"⚠️  {successful}/{total} downloads completed successfully")
            failed_sources = [source for source, success in results.items() if not success]
            logger.warning(f"Failed sources: {failed_sources}")
        
    except Exception as e:
        logger.error(f"Error during download: {str(e)}")
        sys.exit(1)


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
@click.pass_context
def download_for_conversion(ctx, input_type, output_type, minimal, force):
    """Download databases specifically for a gene ID conversion."""
    logger = ctx.obj['logger']
    logger.info(f"Downloading databases for conversion: {input_type} -> {output_type}")
    
    try:
        downloader = DatabaseDownloader(logger=logger)
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
            
    except Exception as e:
        logger.error(f"Error during conversion-specific download: {str(e)}")
        sys.exit(1)


@cli.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.option('--input-type', '-i', required=True, 
              type=click.Choice(SUPPORTED_ID_TYPES),
              help='Type of input gene identifiers')
@click.option('--output-type', '-o', multiple=True, required=True,
              type=click.Choice(SUPPORTED_ID_TYPES),
              help='Type(s) of output gene identifiers (can specify multiple)')
@click.option('--output', '-out', default='output_mapped.csv',
              help='Output file path (default: output_mapped.csv)')
@click.option('--column', '-c', default=None,
              help='Specific column name containing gene IDs (auto-detect if not specified)')
@click.option('--delimiter', '-d', default=',',
              help='Input file delimiter (default: comma)')
@click.option('--include-unmapped', is_flag=True,
              help='Include rows with unmapped identifiers in output')
@click.pass_context
def convert(ctx, input_file, input_type, output_type, output, column, delimiter, include_unmapped):
    """Convert gene identifiers between different formats."""
    logger = ctx.obj['logger']
    logger.info(f"Starting gene ID conversion: {input_type} -> {output_type}")
    
    try:
        # Validate input file
        validate_input_file(input_file, logger)
        
        # Initialize pandas-based converter (no database needed)
        converter = GeneConverter(data_dir="data/databases")
        
        # Check if we have adequate databases for this conversion
        downloader = DatabaseDownloader(logger=logger)
        recommended_sources = downloader.get_recommended_sources_for_conversion(input_type, list(output_type))
        
        # Check which recommended sources are available
        missing_sources = []
        for source in recommended_sources[:3]:  # Check top 3 recommended sources
            source_info = downloader.SOURCES[source]
            filepath = downloader.data_dir / source_info['filename']
            if not filepath.exists():
                # Also check for decompressed version
                decompressed = downloader.data_dir / source_info['filename'].replace('.gz', '')
                if not decompressed.exists():
                    missing_sources.append(source)
        
        if missing_sources:
            logger.warning(f"Missing recommended databases for optimal conversion: {missing_sources}")
            logger.info("To download these databases, run:")
            logger.info(f"  python main.py download --sources {' '.join(missing_sources)}")
            logger.info("Proceeding with available databases...")
        
        # Load input data
        logger.info(f"Loading input file: {input_file}")
        if delimiter == '\\t':
            delimiter = '\t'
        
        df = pd.read_csv(input_file, delimiter=delimiter)
        logger.info(f"Loaded {len(df)} rows from input file")
        
        # Auto-detect column if not specified
        if column is None:
            possible_columns = [col for col in df.columns if any(
                keyword in col.lower() for keyword in ['gene', 'id', 'symbol', 'name']
            )]
            if len(possible_columns) == 1:
                column = possible_columns[0]
                logger.info(f"Auto-detected gene ID column: {column}")
            elif len(possible_columns) > 1:
                logger.warning(f"Multiple possible gene ID columns found: {possible_columns}")
                logger.warning(f"Using first column: {possible_columns[0]}")
                column = possible_columns[0]
            else:
                column = df.columns[0]
                logger.warning(f"Could not auto-detect gene ID column, using first column: {column}")
        
        if column not in df.columns:
            logger.error(f"Column '{column}' not found in input file")
            logger.info(f"Available columns: {list(df.columns)}")
            sys.exit(1)
        
        # Perform conversion
        logger.info(f"Converting {len(df)} gene IDs...")
        result_df = converter.convert_batch(
            gene_ids=df[column].tolist(),
            input_type=input_type,
            output_types=list(output_type),
            include_unmapped=include_unmapped
        )
        
        # Merge with original data
        df_result = df.copy()
        for out_type in output_type:
            df_result[out_type] = result_df[out_type]
        
        # Save results
        df_result.to_csv(output, index=False)
        
        # Log statistics
        total_input = len(df)
        mapped_count = result_df.dropna(subset=list(output_type), how='all').shape[0]
        unmapped_count = total_input - mapped_count
        
        logger.info(f"Conversion completed!")
        logger.info(f"Total input IDs: {total_input}")
        logger.info(f"Successfully mapped: {mapped_count} ({mapped_count/total_input*100:.1f}%)")
        logger.info(f"Unmapped: {unmapped_count} ({unmapped_count/total_input*100:.1f}%)")
        logger.info(f"Output saved to: {output}")
        
    except Exception as e:
        logger.error(f"Error during conversion: {str(e)}")
        sys.exit(1)


@cli.command()
@click.pass_context
def info(ctx):
    """Show information about available databases and supported ID types."""
    logger = ctx.obj['logger']
    
    click.echo("\n🧬 GeneMapKit - Gene ID Mapping Toolkit")
    click.echo("=" * 60)
    
    click.echo("\n📋 Supported Gene ID Types:")
    id_info = [
        ("symbol", "Gene Symbol (HGNC)", "TP53"),
        ("ensembl_gene_id", "Ensembl Gene ID", "ENSG00000141510"),
        ("ensembl_transcript_id", "Ensembl Transcript ID", "ENST00000269305"),
        ("entrez_id", "NCBI Entrez Gene ID", "7157"),
        ("hgnc_id", "HGNC ID", "HGNC:11998"),
        ("refseq_id", "RefSeq ID", "NM_000546"),
        ("uniprot_id", "UniProt ID", "P04637")
    ]
    
    for i, (id_type, description, example) in enumerate(id_info, 1):
        click.echo(f"  {i}. {id_type:20} {description:25} Example: {example}")
    
    click.echo("\n🗄️  Database Sources:")
    sources_info = [
        ("HGNC", "Hugo Gene Nomenclature Committee", "Official gene symbols and cross-references"),
        ("NCBI Gene Info", "National Center for Biotechnology Information", "Entrez IDs, symbols, synonyms"),
        ("Ensembl", "Ensembl Genome Browser", "Gene annotations and genomic coordinates"),
        ("UniProt", "Universal Protein Resource", "Protein-gene mappings")
    ]
    
    for name, full_name, description in sources_info:
        click.echo(f"  • {name:15} {full_name}")
        click.echo(f"    {' '*17} {description}")
    
    click.echo("\n� Supported Conversions:")
    click.echo("  Any ID type ↔ Any other ID type")
    click.echo("  Examples:")
    click.echo("  • Gene Symbol → Ensembl Gene ID")
    click.echo("  • Ensembl Gene ID → Entrez ID")
    click.echo("  • Entrez ID → UniProt ID")
    click.echo("  • HGNC ID → RefSeq ID")
    click.echo("  • Multiple outputs in single run")
    
    click.echo("\n📖 Example Usage:")
    examples = [
        "# Download all databases",
        "genemapkit download",
        "",
        "# Convert gene symbols to Ensembl IDs", 
        "genemapkit convert input.csv -i symbol -o ensembl_gene_id",
        "",
        "# Convert with multiple outputs",
        "genemapkit convert input.csv -i symbol -o ensembl_gene_id -o entrez_id",
        "",
        "# Use tab-delimited input",
        "genemapkit convert input.tsv -i entrez_id -o symbol --delimiter '\\t'",
        "",
        "# Include unmapped genes in output",
        "genemapkit convert input.csv -i symbol -o ensembl_gene_id --include-unmapped"
    ]
    
    for example in examples:
        click.echo(f"  {example}")
    
    click.echo("\n💾 Database Download Options:")
    click.echo("  Core databases (4 sources - optimized for human genes):")
    click.echo("    • HGNC - Gene symbols and cross-references")
    click.echo("    • NCBI - Entrez IDs and gene info")
    click.echo("    • Ensembl GTF - Gene annotations")
    click.echo("    • UniProt Human - Human protein mappings (preferred over general UniProt)")
    click.echo("")
    click.echo("  All databases (7 sources - excludes large general UniProt mapping):")
    click.echo("    • Core sources plus supplementary databases")
    click.echo("    • HGNC JSON, RefSeq Human, BioMart")
    click.echo("    • Note: Large general UniProt mapping excluded by default")
    click.echo("")
    click.echo("  Smart download for specific conversions:")
    click.echo("    genemapkit download-for -i symbol -o ensembl_gene_id")
    click.echo("")
    
    click.echo("\n🎯 Quick Start:")
    examples = [
        "# Download databases for your conversion type",
        "genemapkit download-for -i symbol -o ensembl_gene_id",
        "",
        "# Or download core databases",
        "genemapkit download --sources core",
        "",
        "# Convert your genes",
        "genemapkit convert input.csv -i symbol -o ensembl_gene_id",
        "",
        "# Include unmapped genes in output",
        "genemapkit convert input.csv -i symbol -o ensembl_gene_id --include-unmapped"
    ]
    
    for example in examples:
        click.echo(f"  {example}")
    
    click.echo("\n📊 Sample Data:")
    click.echo("  Sample files available in data/sample/ directory:")
    click.echo("  • gene_symbols.csv - Gene symbols")
    click.echo("  • ensembl_ids.csv - Ensembl gene IDs")
    click.echo("  • entrez_with_metadata.tsv - Entrez IDs with metadata")
    click.echo("  • clinical_data.csv - Clinical data example")
    

@cli.command()
@click.pass_context
def status(ctx):
    """Check status of downloaded databases."""
    logger = ctx.obj['logger']
    
    db_dir = Path("data/databases")
    downloader = DatabaseDownloader()
    
    click.echo("\n📊 Database Status:")
    click.echo("=" * 60)
    
    # Core databases
    click.echo("\n🔥 Core Databases:")
    core_sources = ['hgnc', 'ncbi_gene_info', 'ensembl_gtf', 'uniprot_human']
    
    for source in core_sources:
        info = downloader.get_file_info(source)
        if info:
            click.echo(f"✅ {source.upper()}: Available")
            click.echo(f"   📁 {info['filename']} ({info['size_mb']:.1f} MB)")
            click.echo(f"   📋 {downloader.SOURCES[source]['description']}")
        else:
            click.echo(f"❌ {source.upper()}: Not downloaded")
    
    # Optional databases
    click.echo("\n📚 Optional/Supplementary Databases:")
    optional_sources = ['uniprot_idmapping', 'refseq_human', 'ensembl_biomart_sample']
    
    for source in optional_sources:
        info = downloader.get_file_info(source)
        if info:
            click.echo(f"✅ {source.upper()}: Available")
            click.echo(f"   📁 {info['filename']} ({info['size_mb']:.1f} MB)")
            click.echo(f"   📋 {downloader.SOURCES[source]['description']}")
        else:
            click.echo(f"❌ {source.upper()}: Not downloaded")
    
    # Summary
    all_sources = core_sources + optional_sources
    available_count = sum(1 for source in all_sources if downloader.get_file_info(source))
    
    click.echo(f"\n📈 Summary: {available_count}/{len(all_sources)} databases available")
    
    if available_count == 0:
        click.echo("\n💡 No databases found. Run 'genemapkit download' to download databases")
    elif available_count < len(core_sources):
        click.echo("\n⚠️  Some core databases missing. Run 'genemapkit download --sources core' for basic functionality")
    else:
        click.echo("\n🎉 Core databases available. Ready for gene ID conversion!")
    
    # Database directory info
    if db_dir.exists():
        total_size = sum(f.stat().st_size for f in db_dir.glob('*') if f.is_file())
        click.echo(f"\n💾 Total database size: {total_size / (1024 * 1024):.1f} MB")
        click.echo(f"📂 Database directory: {db_dir.absolute()}")


if __name__ == '__main__':
    cli()
