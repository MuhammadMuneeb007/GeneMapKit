#!/usr/bin/env python3
"""
Comprehensive Database Download Script for GeneMapKit

This script downloads all gene annotation databases from multiple sources
and provides detailed logging of the process.

Usage:
    python download_databases.py [options]
"""

import sys
import logging
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from genemapkit.core.downloader import DatabaseDownloader
from genemapkit.utils.logger import setup_logging


def main():
    """Main download function with comprehensive logging."""
    
    print("🧬 GeneMapKit Comprehensive Database Downloader")
    print("=" * 60)
    
    # Setup logging
    logger = setup_logging('INFO')
    
    # Create downloader
    downloader = DatabaseDownloader(logger=logger)
    
    # Display available databases
    logger.info("📊 Available Database Sources:")
    logger.info("-" * 40)
    
    for source, info in downloader.SOURCES.items():
        logger.info(f"🗄️  {source.upper()}")
        logger.info(f"   Description: {info['description']}")
        logger.info(f"   Contains: {', '.join(info['contains'])}")
        logger.info(f"   URL: {info['url']}")
        logger.info(f"   File: {info['filename']}")
        logger.info("")
    
    # Ask user for download preference
    print("\nDownload Options:")
    print("1. Core databases only (HGNC, NCBI, Ensembl GTF, UniProt Human) - 4 databases")
    print("2. All databases (includes supplementary sources, excludes large general UniProt) - 7 databases")
    print("3. Custom selection")
    print("4. Download for specific conversion")
    print("5. Exit")
    
    try:
        choice = input("\nEnter your choice (1-5): ").strip()
        
        if choice == "1":
            logger.info("Downloading core databases...")
            core_sources = ['hgnc', 'ncbi_gene_info', 'ensembl_gtf', 'uniprot_human']
            results = downloader.download_specific_sources(core_sources, force_download=False)
            
        elif choice == "2":
            logger.info("Downloading all databases...")
            results = downloader.download_all(force_download=False, include_optional=True)
            
        elif choice == "3":
            logger.info("Available sources:")
            sources = list(downloader.SOURCES.keys())
            for i, source in enumerate(sources, 1):
                source_info = downloader.SOURCES[source]
                print(f"{i}. {source} - {source_info['description']}")
            
            selected_indices = input("Enter source numbers (comma-separated): ").strip()
            try:
                indices = [int(x.strip()) - 1 for x in selected_indices.split(',')]
                selected_sources = [sources[i] for i in indices if 0 <= i < len(sources)]
                
                if selected_sources:
                    logger.info(f"Downloading selected sources: {selected_sources}")
                    results = downloader.download_specific_sources(selected_sources, force_download=False)
                else:
                    logger.error("No valid sources selected")
                    return
            except (ValueError, IndexError):
                logger.error("Invalid selection")
                return
                
        elif choice == "4":
            print("\nAvailable gene ID types:")
            id_types = ['symbol', 'ensembl_gene_id', 'ensembl_transcript_id', 
                       'entrez_id', 'hgnc_id', 'refseq_id', 'uniprot_id']
            for i, id_type in enumerate(id_types, 1):
                print(f"{i}. {id_type}")
            
            try:
                input_idx = int(input("Enter input type number: ").strip()) - 1
                output_indices = input("Enter output type numbers (comma-separated): ").strip()
                output_idx_list = [int(x.strip()) - 1 for x in output_indices.split(',')]
                
                input_type = id_types[input_idx]
                output_types = [id_types[i] for i in output_idx_list if 0 <= i < len(id_types)]
                
                if output_types:
                    logger.info(f"Downloading for conversion: {input_type} -> {output_types}")
                    results = downloader.download_for_conversion(input_type, output_types, 
                                                               force_download=False, minimal=False)
                else:
                    logger.error("No valid output types selected")
                    return
            except (ValueError, IndexError):
                logger.error("Invalid selection")
                return
                
        elif choice == "5":
            logger.info("Download cancelled by user")
            return
            
        else:
            logger.error("Invalid choice")
            return
        
        # Display results
        successful = sum(1 for success in results.values() if success)
        total = len(results)
        
        logger.info("\n" + "=" * 60)
        logger.info("DOWNLOAD SUMMARY")
        logger.info("=" * 60)
        
        for source, success in results.items():
            status = "✅ SUCCESS" if success else "❌ FAILED"
            logger.info(f"{source.upper()}: {status}")
        
        logger.info(f"\nOverall: {successful}/{total} downloads completed successfully")
        
        if successful == total:
            logger.info("🎉 All downloads completed successfully!")
            logger.info("You can now run gene ID conversions with GeneMapKit.")
        else:
            logger.warning("⚠️  Some downloads failed. You may have limited functionality.")
        
    except KeyboardInterrupt:
        logger.info("\nDownload interrupted by user")
    except Exception as e:
        logger.error(f"Unexpected error: {str(e)}")


if __name__ == '__main__':
    main()
