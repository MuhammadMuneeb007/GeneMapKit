"""
Database Downloader Module

Handles downloading gene annotation databases from multiple sources:
- HGNC (Hugo Gene Nomenclature Committee)
- NCBI Gene Info
- Ensembl
- UniProt
"""

import os
import gzip
import shutil
import requests
from pathlib import Path
from typing import Dict, Optional, List
from urllib.parse import urlparse
import logging
from tqdm import tqdm
import time


class DatabaseDownloader:
    """Downloads and manages gene annotation databases from multiple sources."""
    
    # Database source URLs - Comprehensive collection for gene ID mapping
    SOURCES = {
        'hgnc': {
            'url': 'https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt',
            'filename': 'hgnc_complete_set.txt',
            'description': 'HGNC complete gene set with symbols, IDs, and cross-references',
            'contains': ['symbol', 'hgnc_id', 'ensembl_gene_id', 'entrez_id', 'refseq_id', 'uniprot_id', 'aliases'],
            'priority': 1,
            'best_for': ['symbol', 'hgnc_id', 'general_mapping']
        } ,
        'ncbi_gene_info': {
            'url': 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz',
            'filename': 'Homo_sapiens.gene_info.gz',
            'description': 'NCBI Gene Info for Homo sapiens with Entrez IDs, symbols, and synonyms',
            'contains': ['entrez_id', 'symbol', 'synonyms', 'chromosome', 'description'],
            'priority': 1,
            'best_for': ['entrez_id', 'general_mapping']
        },
        'ensembl_gtf': {
            'url': 'https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz',
            'filename': 'Homo_sapiens.GRCh38.111.gtf.gz',
            'description': 'Ensembl GTF annotation with gene IDs, transcript IDs, and genomic coordinates',
            'contains': ['ensembl_gene_id', 'ensembl_transcript_id', 'symbol', 'chromosome', 'coordinates'],
            'priority': 1,
            'best_for': ['ensembl_gene_id', 'ensembl_transcript_id']
        },
        # 'uniprot_idmapping': {
        #     'url': 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz',
        #     'filename': 'uniprot_idmapping_selected.tab.gz',
        #     'description': 'UniProt ID mapping to gene symbols, Ensembl IDs, and Entrez IDs (LARGE FILE - SKIPPED BY DEFAULT)',
        #     'contains': ['uniprot_id', 'symbol', 'ensembl_gene_id', 'entrez_id'],
        #     'priority': 4,
        #     'best_for': [],
        #     'skip_by_default': True
        # },
        'uniprot_human': {
            'url': 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz',
            'filename': 'HUMAN_9606_idmapping_selected.tab.gz',
            'description': 'UniProt ID mapping specifically for human proteins',
            'contains': ['uniprot_id', 'symbol', 'ensembl_gene_id', 'entrez_id'],
            'priority': 2,
            'best_for': ['uniprot_id']
        },
        'refseq_human': {
            'url': 'https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.1.rna.fna.gz',
            'filename': 'human_refseq_rna.fna.gz',
            'description': 'RefSeq human mRNA sequences with gene annotations',
            'contains': ['refseq_id', 'symbol', 'entrez_id'],
            'priority': 2,
            'best_for': ['refseq_id']
        },
        'ensembl_biomart_sample': {
            'url': 'https://ftp.ensembl.org/pub/release-111/tsv/homo_sapiens/Homo_sapiens.GRCh38.111.entrez.tsv.gz',
            'filename': 'ensembl_entrez_mapping.tsv.gz',
            'description': 'Ensembl to Entrez ID mapping from BioMart export',
            'contains': ['ensembl_gene_id', 'entrez_id', 'symbol'],
            'priority': 3,
            'best_for': ['ensembl_gene_id', 'entrez_id']
        }
    }
    
    def __init__(self, data_dir: str = "data/databases", logger: Optional[logging.Logger] = None):
        """
        Initialize the database downloader.
        
        Args:
            data_dir: Directory to store downloaded databases
            logger: Logger instance for logging messages
        """
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.logger = logger or logging.getLogger(__name__)
    
    def get_recommended_sources_for_conversion(self, input_type: str, output_types: List[str]) -> List[str]:
        """
        Get recommended database sources for a specific conversion.
        
        Args:
            input_type: Input gene ID type
            output_types: List of output gene ID types
            
        Returns:
            List of recommended source names, ordered by priority
        """
        all_types = [input_type] + output_types
        recommended = {}
        
        # Find sources that support the required ID types
        for source_name, source_info in self.SOURCES.items():
            # Skip sources marked as skip_by_default (like large uniprot_idmapping)
            if source_info.get('skip_by_default', False):
                continue
                
            score = 0
            supported_types = 0
            
            # Check how many required types this source supports
            for id_type in all_types:
                if id_type in source_info['contains']:
                    supported_types += 1
                    score += 10  # Base score for supporting the type
                
                # Bonus score if this source is specifically good for this type
                if id_type in source_info.get('best_for', []):
                    score += 20
            
            # Only recommend sources that support at least one required type
            if supported_types > 0:
                # Adjust score by priority (lower priority number = higher score)
                priority_bonus = (4 - source_info.get('priority', 3)) * 5
                final_score = score + priority_bonus
                recommended[source_name] = final_score
        
        # Sort by score (highest first) and return source names
        sorted_sources = sorted(recommended.items(), key=lambda x: x[1], reverse=True)
        return [source for source, score in sorted_sources]
    
    def get_minimal_sources_for_conversion(self, input_type: str, output_types: List[str]) -> List[str]:
        """
        Get the minimal set of sources needed for a conversion.
        
        Args:
            input_type: Input gene ID type
            output_types: List of output gene ID types
            
        Returns:
            Minimal list of source names needed for the conversion
        """
        all_types = set([input_type] + output_types)
        covered_types = set()
        selected_sources = []
        
        # Get sources ranked by their suitability
        recommended = self.get_recommended_sources_for_conversion(input_type, output_types)
        
        # Greedily select sources until all types are covered
        for source in recommended:
            source_info = self.SOURCES[source]
            source_types = set(source_info['contains'])
            
            # Check if this source adds new coverage
            new_coverage = source_types & all_types - covered_types
            
            if new_coverage:
                selected_sources.append(source)
                covered_types.update(new_coverage)
                
                # Stop if we've covered all required types
                if all_types.issubset(covered_types):
                    break
        
        return selected_sources
        
    def download_for_conversion(self, input_type: str, output_types: List[str], 
                               force_download: bool = False, minimal: bool = False) -> Dict[str, bool]:
        """
        Download databases specifically needed for a conversion.
        
        Args:
            input_type: Input gene ID type
            output_types: List of output gene ID types
            force_download: Force re-download even if files exist
            minimal: If True, download only minimal required sources
            
        Returns:
            Dictionary with source names as keys and success status as values
        """
        if minimal:
            sources_to_download = self.get_minimal_sources_for_conversion(input_type, output_types)
        else:
            sources_to_download = self.get_recommended_sources_for_conversion(input_type, output_types)
        
        self.logger.info(f"Downloading databases for conversion: {input_type} → {output_types}")
        self.logger.info(f"Recommended sources: {sources_to_download}")
        
        return self.download_specific_sources(sources_to_download, force_download)
    
    def download_source(self, source: str, force_download: bool = False) -> bool:
        """
        Download a specific database source.
        
        Args:
            source: Source name ('hgnc', 'ncbi', 'ensembl', 'uniprot')
            force_download: Force re-download even if file exists
            
        Returns:
            True if download successful, False otherwise
        """
        if source not in self.SOURCES:
            self.logger.error(f"Unknown source: {source}")
            return False
            
        source_info = self.SOURCES[source]
        url = source_info['url']
        filename = source_info['filename']
        filepath = self.data_dir / filename
        
        # Check if file already exists
        if filepath.exists() and not force_download:
            self.logger.info(f"File {filename} already exists. Use --force to re-download.")
            return True
            
        self.logger.info(f"Downloading {source_info['description']}...")
        self.logger.info(f"URL: {url}")
        
        try:
            # Download with progress bar
            response = requests.get(url, stream=True)
            response.raise_for_status()
            
            total_size = int(response.headers.get('content-length', 0))
            
            with open(filepath, 'wb') as file, tqdm(
                desc=f"Downloading {source}",
                total=total_size,
                unit='iB',
                unit_scale=True,
                unit_divisor=1024,
            ) as progress_bar:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        file.write(chunk)
                        progress_bar.update(len(chunk))
            
            # Decompress if gzipped
            if filename.endswith('.gz'):
                self._decompress_file(filepath)
            
            # Parse JSON files if needed
            if filename.endswith('.json'):
                self._validate_json_file(filepath)
            
            self.logger.info(f"Successfully downloaded {filename}")
            return True
            
        except Exception as e:
            self.logger.error(f"Error downloading {source}: {str(e)}")
            if filepath.exists():
                filepath.unlink()  # Remove partial download
            return False
    
    def _validate_json_file(self, filepath: Path) -> None:
        """Validate JSON file structure."""
        import json
        
        self.logger.info(f"Validating JSON file: {filepath.name}")
        
        try:
            with open(filepath, 'r') as f:
                data = json.load(f)
            
            if isinstance(data, dict) and 'response' in data:
                docs = data.get('response', {}).get('docs', [])
                self.logger.info(f"JSON file contains {len(docs)} gene records")
            else:
                self.logger.warning("JSON file structure may not be as expected")
                
        except Exception as e:
            self.logger.error(f"Error validating JSON file: {str(e)}")
    
    def download_specific_sources(self, source_list: List[str], force_download: bool = False) -> Dict[str, bool]:
        """
        Download specific database sources.
        
        Args:
            source_list: List of source names to download
            force_download: Force re-download even if files exist
            
        Returns:
            Dictionary with source names as keys and success status as values
        """
        results = {}
        failed_details = []
        
        for source in source_list:
            if source not in self.SOURCES:
                self.logger.error(f"Unknown source: {source}")
                self.logger.info(f"Available sources: {list(self.SOURCES.keys())}")
                results[source] = False
                failed_details.append(f"{source}: Unknown source")
                continue
            
            self.logger.info(f"Downloading {source}...")
            try:
                success = self.download_source(source, force_download)
                results[source] = success
                
                if not success:
                    failed_details.append(f"{source}: Download failed")
                else:
                    self.logger.info(f"✅ {source} downloaded successfully")
                    
            except Exception as e:
                self.logger.error(f"❌ Error downloading {source}: {str(e)}")
                results[source] = False
                failed_details.append(f"{source}: {str(e)}")
        
        # Summary
        successful = sum(1 for success in results.values() if success)
        total = len(results)
        
        if successful < total:
            self.logger.warning(f"Failed downloads ({total - successful}/{total}):")
            for detail in failed_details:
                self.logger.warning(f"  • {detail}")
        
        return results
    
    def _decompress_file(self, filepath: Path) -> None:
        """Decompress a gzipped file."""
        decompressed_path = filepath.with_suffix('')
        
        self.logger.info(f"Decompressing {filepath.name}...")
        
        try:
            with gzip.open(filepath, 'rb') as f_in:
                with open(decompressed_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            
            # Remove the compressed file
            filepath.unlink()
            self.logger.info(f"Decompressed to {decompressed_path.name}")
            
        except Exception as e:
            self.logger.error(f"Error decompressing {filepath.name}: {str(e)}")
            if decompressed_path.exists():
                decompressed_path.unlink()
    
    def download_all(self, force_download: bool = False, include_optional: bool = True) -> Dict[str, bool]:
        """
        Download all database sources.
        
        Args:
            force_download: Force re-download even if files exist
            include_optional: Include optional/supplementary databases
            
        Returns:
            Dictionary with source names as keys and success status as values
        """
        results = {}
        
        # Core databases (essential for basic functionality)
        core_sources = ['hgnc', 'ncbi_gene_info', 'ensembl_gtf', 'uniprot_human']
        
        # Optional databases (for enhanced functionality) - excluding large uniprot_idmapping
        optional_sources = ['refseq_human', 'ensembl_biomart_sample']
        
        sources_to_download = core_sources
        if include_optional:
            sources_to_download.extend(optional_sources)
        
        self.logger.info(f"Starting download of {len(sources_to_download)} database sources...")
        
        for i, source in enumerate(sources_to_download, 1):
            self.logger.info(f"\n{'='*60}")
            self.logger.info(f"Processing {source.upper()} database ({i}/{len(sources_to_download)})")
            self.logger.info(f"{'='*60}")
            
            source_info = self.SOURCES[source]
            self.logger.info(f"Description: {source_info['description']}")
            self.logger.info(f"Contains: {', '.join(source_info['contains'])}")
            
            results[source] = self.download_source(source, force_download)
            
            if results[source]:
                self.logger.info(f"✅ {source.upper()} download completed")
            else:
                self.logger.error(f"❌ {source.upper()} download failed")
            
            # Small delay between downloads to be respectful to servers
            time.sleep(2)
        
        # Summary
        successful = sum(1 for success in results.values() if success)
        total = len(results)
        
        self.logger.info(f"\n{'='*60}")
        self.logger.info(f"DOWNLOAD SUMMARY: {successful}/{total} sources completed successfully")
        self.logger.info(f"{'='*60}")
        
        return results
    
    def get_file_info(self, source: str) -> Optional[Dict]:
        """
        Get information about a downloaded file.
        
        Args:
            source: Source name
            
        Returns:
            Dictionary with file information or None if not found
        """
        if source not in self.SOURCES:
            return None
            
        source_info = self.SOURCES[source]
        filename = source_info['filename']
        
        # Check for both compressed and decompressed versions
        filepath = self.data_dir / filename
        if not filepath.exists():
            filepath = self.data_dir / filename.replace('.gz', '')
            
        if not filepath.exists():
            return None
            
        stat = filepath.stat()
        return {
            'filename': filepath.name,
            'size_mb': stat.st_size / (1024 * 1024),
            'modified': stat.st_mtime,
            'description': source_info['description']
        }
    
    def list_downloaded_files(self) -> Dict[str, Dict]:
        """List all downloaded database files with their information."""
        files_info = {}
        
        for source in self.SOURCES:
            info = self.get_file_info(source)
            if info:
                files_info[source] = info
                
        return files_info
    
    def cleanup_old_files(self, keep_days: int = 30) -> None:
        """
        Remove database files older than specified days.
        
        Args:
            keep_days: Number of days to keep files
        """
        import time
        
        cutoff_time = time.time() - (keep_days * 24 * 60 * 60)
        
        for file_path in self.data_dir.glob('*'):
            if file_path.is_file() and file_path.stat().st_mtime < cutoff_time:
                self.logger.info(f"Removing old file: {file_path.name}")
                file_path.unlink()
