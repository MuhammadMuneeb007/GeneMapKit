"""
Gene Database Module

Handles loading and managing gene annotation data from multiple sources.
Provides unified interface for querying gene mappings.
"""

import sqlite3
import pandas as pd
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set
import logging
from collections import defaultdict


class GeneDatabase:
    """Manages gene annotation databases and provides mapping functionality."""
    
    def __init__(self, data_dir: str = "data/databases", logger: Optional[logging.Logger] = None):
        """
        Initialize the gene database.
        
        Args:
            data_dir: Directory containing downloaded database files
            logger: Logger instance for logging messages
        """
        self.data_dir = Path(data_dir)
        self.logger = logger or logging.getLogger(__name__)
        self.db_path = self.data_dir / "genemapkit.db"
        
        # Mapping tables
        self.mappings = {}
        self.loaded_sources = set()
        
        # Initialize database
        self._init_database()
        self._load_all_sources()
    
    def _init_database(self) -> None:
        """Initialize SQLite database for fast lookups."""
        self.data_dir.mkdir(parents=True, exist_ok=True)
        
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        # Create main mapping table
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS gene_mappings (
                symbol TEXT,
                ensembl_gene_id TEXT,
                ensembl_transcript_id TEXT,
                entrez_id TEXT,
                hgnc_id TEXT,
                refseq_id TEXT,
                uniprot_id TEXT,
                aliases TEXT,
                description TEXT,
                chromosome TEXT,
                source TEXT,
                UNIQUE(symbol, ensembl_gene_id, entrez_id) ON CONFLICT IGNORE
            )
        ''')
        
        # Create indexes for fast lookups
        indexes = [
            'CREATE INDEX IF NOT EXISTS idx_symbol ON gene_mappings(symbol)',
            'CREATE INDEX IF NOT EXISTS idx_ensembl ON gene_mappings(ensembl_gene_id)',
            'CREATE INDEX IF NOT EXISTS idx_entrez ON gene_mappings(entrez_id)',
            'CREATE INDEX IF NOT EXISTS idx_hgnc ON gene_mappings(hgnc_id)',
            'CREATE INDEX IF NOT EXISTS idx_refseq ON gene_mappings(refseq_id)',
            'CREATE INDEX IF NOT EXISTS idx_uniprot ON gene_mappings(uniprot_id)'
        ]
        
        for index_sql in indexes:
            cursor.execute(index_sql)
        
        conn.commit()
        conn.close()
        
        self.logger.info(f"Database initialized at: {self.db_path}")
    
    def clear_database(self) -> None:
        """Clear all data from the database."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        cursor.execute('DELETE FROM gene_mappings')
        conn.commit()
        conn.close()
        self.loaded_sources.clear()
        self.logger.info("Database cleared")
    
    def _load_all_sources(self) -> None:
        """Load data from all available database sources."""
        # Check which files are available - updated for comprehensive sources
        available_files = {
            'hgnc': self.data_dir / 'hgnc_complete_set.txt',
            'ncbi_gene_info': self.data_dir / 'Homo_sapiens.gene_info',
            'ensembl_gtf': self.data_dir / 'Homo_sapiens.GRCh38.111.gtf',
            'uniprot_human': self.data_dir / 'HUMAN_9606_idmapping_selected.tab',
            'uniprot_idmapping': self.data_dir / 'uniprot_idmapping_selected.tab',
            'refseq_human': self.data_dir / 'human_refseq_rna.fna',
            'ensembl_biomart': self.data_dir / 'ensembl_entrez_mapping.tsv'
        }
        
        for source, filepath in available_files.items():
            if filepath.exists():
                try:
                    self.logger.info(f"Loading {source.upper()} database...")
                    self._load_source(source, filepath)
                    self.loaded_sources.add(source)
                    self.logger.info(f"✅ Successfully loaded {source}")
                except Exception as e:
                    self.logger.error(f"❌ Error loading {source}: {str(e)}")
        
        if not self.loaded_sources:
            self.logger.warning("No database sources loaded. Run 'genemapkit download' first.")
        else:
            self.logger.info(f"Loaded {len(self.loaded_sources)} database sources: {', '.join(self.loaded_sources)}")
    
    def _load_source(self, source: str, filepath: Path) -> None:
        """Load data from a specific source file."""
        if source == 'hgnc':
            self._load_hgnc(filepath)
        elif source == 'ncbi_gene_info':
            self._load_ncbi(filepath)
        elif source == 'ensembl_gtf':
            self._load_ensembl(filepath)
        elif source in ['uniprot_human', 'uniprot_idmapping']:
            self._load_uniprot(filepath)
        elif source == 'refseq_human':
            self._load_refseq(filepath)
        elif source == 'ensembl_biomart':
            self._load_ensembl_biomart(filepath)
    
    def _load_hgnc(self, filepath: Path) -> None:
        """Load HGNC data."""
        try:
            df = pd.read_csv(filepath, sep='\t', dtype=str, low_memory=False)
            
            # Map HGNC columns to our standard format
            column_mapping = {
                'symbol': 'symbol',
                'hgnc_id': 'hgnc_id',
                'ensembl_gene_id': 'ensembl_gene_id',
                'entrez_id': 'entrez_id',
                'refseq_accession': 'refseq_id',
                'uniprot_ids': 'uniprot_id',
                'alias_symbol': 'aliases',
                'name': 'description',
                'location': 'chromosome'
            }
            
            # Select and rename columns
            available_cols = {k: v for k, v in column_mapping.items() if k in df.columns}
            if not available_cols:
                self.logger.warning("No recognized columns found in HGNC data")
                return
                
            df_mapped = df[list(available_cols.keys())].rename(columns=available_cols)
            df_mapped['source'] = 'hgnc'
            
            # Remove completely empty rows
            df_mapped = df_mapped.dropna(how='all')
            
            # Insert into database
            self._insert_mappings(df_mapped)
            self.logger.info(f"Loaded {len(df_mapped)} entries from HGNC")
            
        except Exception as e:
            self.logger.error(f"Error loading HGNC data: {str(e)}")
            import traceback
            self.logger.debug(f"Full traceback: {traceback.format_exc()}")
    
    def _load_ncbi(self, filepath: Path) -> None:
        """Load NCBI Gene Info data."""
        try:
            df = pd.read_csv(filepath, sep='\t', dtype=str, comment='#', low_memory=False)
            
            # Map NCBI columns
            column_mapping = {
                'GeneID': 'entrez_id',
                'Symbol': 'symbol',
                'Synonyms': 'aliases',
                'description': 'description',
                'chromosome': 'chromosome'
            }
            
            available_cols = {k: v for k, v in column_mapping.items() if k in df.columns}
            df_mapped = df[list(available_cols.keys())].rename(columns=available_cols)
            
            # Extract Ensembl Gene IDs from dbXrefs column if available
            if 'dbXrefs' in df.columns:
                ensembl_ids = []
                hgnc_ids = []
                for dbxrefs in df['dbXrefs']:
                    ensembl_id = None
                    hgnc_id = None
                    if pd.notna(dbxrefs) and dbxrefs != '-':
                        # Parse dbXrefs like "MIM:191170|HGNC:HGNC:11998|Ensembl:ENSG00000141510|AllianceGenome:HGNC:11998"
                        for ref in dbxrefs.split('|'):
                            if ref.startswith('Ensembl:ENSG'):
                                ensembl_id = ref.replace('Ensembl:', '')
                            elif ref.startswith('HGNC:HGNC:'):
                                hgnc_id = ref.replace('HGNC:', '')
                    ensembl_ids.append(ensembl_id)
                    hgnc_ids.append(hgnc_id)
                
                df_mapped['ensembl_gene_id'] = ensembl_ids
                df_mapped['hgnc_id'] = hgnc_ids
            
            df_mapped['source'] = 'ncbi_gene_info'
            
            self._insert_mappings(df_mapped)
            self.logger.info(f"Loaded {len(df_mapped)} entries from NCBI")
            
        except Exception as e:
            self.logger.error(f"Error loading NCBI data: {str(e)}")
            import traceback
            self.logger.debug(f"Full traceback: {traceback.format_exc()}")
    
    def _load_ensembl(self, filepath: Path) -> None:
        """Load Ensembl GTF data."""
        try:
            gene_records = {}  # Store by gene_id to merge gene and transcript info
            transcript_mappings = {}  # Map transcript_id to gene_id
            
            with open(filepath, 'r') as f:
                for line_num, line in enumerate(f):
                    if line.startswith('#'):
                        continue
                    
                    if line_num > 500000:  # Increased limit to capture more transcript data
                        break
                        
                    fields = line.strip().split('\t')
                    if len(fields) < 9:
                        continue
                    
                    feature_type = fields[2]
                    attr_str = fields[8]
                    
                    # Process gene entries
                    if feature_type == 'gene':
                        gene_id = self._extract_gtf_attribute(attr_str, 'gene_id')
                        gene_name = self._extract_gtf_attribute(attr_str, 'gene_name')
                        gene_biotype = self._extract_gtf_attribute(attr_str, 'gene_biotype')
                        
                        if gene_id and gene_name:
                            gene_records[gene_id] = {
                                'ensembl_gene_id': gene_id,
                                'symbol': gene_name,
                                'description': gene_biotype,
                                'chromosome': fields[0],
                                'source': 'ensembl'
                            }
                    
                    # Process transcript entries to get transcript-to-gene mappings
                    elif feature_type == 'transcript':
                        gene_id = self._extract_gtf_attribute(attr_str, 'gene_id')
                        transcript_id = self._extract_gtf_attribute(attr_str, 'transcript_id')
                        gene_name = self._extract_gtf_attribute(attr_str, 'gene_name')
                        
                        if gene_id and transcript_id:
                            transcript_mappings[transcript_id] = gene_id
                            
                            # Also create/update gene record if not exists
                            if gene_id not in gene_records and gene_name:
                                gene_records[gene_id] = {
                                    'ensembl_gene_id': gene_id,
                                    'symbol': gene_name,
                                    'chromosome': fields[0],
                                    'source': 'ensembl'
                                }
            
            # Create final records with transcript mappings
            final_records = []
            
            # Add gene records
            for gene_id, gene_record in gene_records.items():
                final_records.append(gene_record)
            
            # Add transcript-to-gene mappings
            for transcript_id, gene_id in transcript_mappings.items():
                if gene_id in gene_records:
                    # Create a record that maps transcript to gene info
                    transcript_record = gene_records[gene_id].copy()
                    transcript_record['ensembl_transcript_id'] = transcript_id
                    final_records.append(transcript_record)
            
            if final_records:
                df_mapped = pd.DataFrame(final_records)
                self._insert_mappings(df_mapped)
                self.logger.info(f"Loaded {len(df_mapped)} entries from Ensembl ({len(gene_records)} genes, {len(transcript_mappings)} transcripts)")
            
        except Exception as e:
            self.logger.error(f"Error loading Ensembl data: {str(e)}")
            import traceback
            self.logger.debug(f"Full traceback: {traceback.format_exc()}")
    
    def _load_uniprot(self, filepath: Path) -> None:
        """Load UniProt ID mapping data."""
        try:
            self.logger.info(f"Loading UniProt data from {filepath.name}...")
            
            # Check if it's the large general UniProt file (skip if too large)
            file_size_mb = filepath.stat().st_size / (1024 * 1024)
            if file_size_mb > 1000:  # Skip files larger than 1GB
                self.logger.warning(f"Skipping large UniProt file {filepath.name} ({file_size_mb:.1f}MB)")
                return
            
            # UniProt ID mapping format: UniProtKB-AC, ID_type, ID
            # Use python engine for better error handling with malformed data
            df = pd.read_csv(filepath, sep='\t', header=None, dtype=str,
                           names=['uniprot_id', 'id_type', 'id_value'],
                           engine='python', on_bad_lines='skip')
            
            self.logger.info(f"Read {len(df)} UniProt mapping entries")
            
            # Group by UniProt ID and collect mappings
            uniprot_mappings = defaultdict(dict)
            
            processed = 0
            for _, row in df.iterrows():
                uniprot_id = row['uniprot_id']
                id_type = row['id_type']
                id_value = row['id_value']
                
                if pd.isna(uniprot_id) or pd.isna(id_type) or pd.isna(id_value):
                    continue
                
                if id_type == 'Gene_Name':
                    uniprot_mappings[uniprot_id]['symbol'] = id_value
                elif id_type == 'Ensembl':
                    uniprot_mappings[uniprot_id]['ensembl_gene_id'] = id_value
                elif id_type == 'GeneID':
                    uniprot_mappings[uniprot_id]['entrez_id'] = id_value
                
                processed += 1
                if processed % 100000 == 0:
                    self.logger.info(f"Processed {processed} UniProt entries...")
            
            # Convert to DataFrame
            records = []
            for uniprot_id, mappings in uniprot_mappings.items():
                record = {'uniprot_id': uniprot_id, 'source': 'uniprot'}
                record.update(mappings)
                records.append(record)
            
            if records:
                df_mapped = pd.DataFrame(records)
                self._insert_mappings(df_mapped)
                self.logger.info(f"Loaded {len(df_mapped)} entries from UniProt")
            
        except Exception as e:
            self.logger.error(f"Error loading UniProt data: {str(e)}")
    
    def _load_refseq(self, filepath: Path) -> None:
        """Load RefSeq FASTA data."""
        try:
            records = []
            
            with open(filepath, 'r') as f:
                for line_num, line in enumerate(f):
                    if line.startswith('>'):
                        # Parse FASTA header for RefSeq ID and gene info
                        header = line.strip()[1:]  # Remove '>'
                        parts = header.split('|')
                        
                        if len(parts) >= 2:
                            refseq_id = parts[1]  # RefSeq accession
                            
                            # Extract gene symbol if present in description
                            description = header.split(' ', 1)[1] if ' ' in header else ''
                            
                            record = {
                                'refseq_id': refseq_id,
                                'description': description,
                                'source': 'refseq'
                            }
                            records.append(record)
                    
                    # Limit for demo - remove in production
                    if line_num > 50000:
                        break
            
            if records:
                df_mapped = pd.DataFrame(records)
                self._insert_mappings(df_mapped)
                self.logger.info(f"Loaded {len(df_mapped)} entries from RefSeq")
            
        except Exception as e:
            self.logger.error(f"Error loading RefSeq data: {str(e)}")
    
    def _load_ensembl_biomart(self, filepath: Path) -> None:
        """Load Ensembl BioMart export data."""
        try:
            df = pd.read_csv(filepath, sep='\t', dtype=str, low_memory=False)
            
            # Map BioMart columns (adjust based on actual export format)
            column_mapping = {
                'Gene stable ID': 'ensembl_gene_id',
                'NCBI gene (formerly Entrezgene) ID': 'entrez_id',
                'Gene name': 'symbol',
                'Gene description': 'description'
            }
            
            available_cols = {k: v for k, v in column_mapping.items() if k in df.columns}
            df_mapped = df[list(available_cols.keys())].rename(columns=available_cols)
            df_mapped['source'] = 'ensembl_biomart'
            
            self._insert_mappings(df_mapped)
            self.logger.info(f"Loaded {len(df_mapped)} entries from Ensembl BioMart")
            
        except Exception as e:
            self.logger.error(f"Error loading Ensembl BioMart data: {str(e)}")
    
    def _extract_gtf_attribute(self, attr_string: str, attribute: str) -> Optional[str]:
        """Extract attribute value from GTF attribute string."""
        pattern = rf'{attribute}\s+"([^"]+)"'
        match = re.search(pattern, attr_string)
        return match.group(1) if match else None
    
    def _insert_mappings(self, df: pd.DataFrame) -> None:
        """Insert gene mappings into the database."""
        conn = sqlite3.connect(self.db_path)
        
        # Fill missing columns with None
        required_columns = [
            'symbol', 'ensembl_gene_id', 'ensembl_transcript_id', 'entrez_id',
            'hgnc_id', 'refseq_id', 'uniprot_id', 'aliases', 'description',
            'chromosome', 'source'
        ]
        
        for col in required_columns:
            if col not in df.columns:
                df[col] = None
        
        # Clean data and handle duplicates
        df_to_insert = df[required_columns].copy()
        
        # Remove rows where all key fields are None (would violate primary key)
        key_fields = ['symbol', 'ensembl_gene_id', 'entrez_id']
        # Only drop if ALL key fields are None/empty - be more permissive
        def has_valid_key(row):
            return any(pd.notna(row[field]) and str(row[field]).strip() not in ['', 'nan', 'None', '-'] 
                      for field in key_fields)
        
        df_to_insert = df_to_insert[df_to_insert.apply(has_valid_key, axis=1)]
        
        # Handle duplicates by keeping the first occurrence
        # Use a combination of key fields to identify duplicates
        df_to_insert['_key'] = (
            df_to_insert['symbol'].fillna('') + '|' +
            df_to_insert['ensembl_gene_id'].fillna('') + '|' +
            df_to_insert['entrez_id'].fillna('')
        )
        df_to_insert = df_to_insert.drop_duplicates(subset=['_key'], keep='first')
        df_to_insert = df_to_insert.drop(columns=['_key'])
        
        if len(df_to_insert) == 0:
            self.logger.warning("No valid records to insert after cleaning")
            conn.close()
            return
        
        # Insert in smaller chunks to avoid "too many SQL variables" error
        chunk_size = 500  # SQLite limit is ~999 variables, so 500 rows * 11 columns = ~5500 variables
        
        try:
            for i in range(0, len(df_to_insert), chunk_size):
                chunk = df_to_insert.iloc[i:i + chunk_size]
                # Use INSERT OR IGNORE to handle any remaining duplicates gracefully
                chunk.to_sql('gene_mappings', conn, if_exists='append', 
                            index=False, method='multi')
        except sqlite3.IntegrityError as e:
            self.logger.warning(f"Some duplicate records were skipped: {str(e)}")
        except Exception as e:
            self.logger.error(f"Error inserting data: {str(e)}")
            raise
        
        conn.close()
    
    def lookup_gene(self, gene_id: str, id_type: str) -> List[Dict]:
        """
        Look up gene information by ID.
        
        Args:
            gene_id: Gene identifier
            id_type: Type of identifier (symbol, ensembl_gene_id, etc.)
            
        Returns:
            List of matching gene records
        """
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        
        # Handle special cases for ID lookup
        if id_type == 'symbol':
            # Also check aliases for symbols
            query = '''
                SELECT * FROM gene_mappings 
                WHERE symbol = ? OR aliases LIKE ?
            '''
            cursor = conn.execute(query, (gene_id, f'%{gene_id}%'))
        elif id_type == 'ensembl_transcript_id':
            # Special handling for transcript IDs - try multiple approaches
            # Method 1: Direct match
            query = f'SELECT * FROM gene_mappings WHERE {id_type} = ?'
            cursor = conn.execute(query, (gene_id,))
            results = [dict(row) for row in cursor.fetchall()]
            
            if not results:
                # Method 2: Try finding in aliases or description
                query = '''
                    SELECT * FROM gene_mappings 
                    WHERE aliases LIKE ? OR description LIKE ?
                '''
                cursor = conn.execute(query, (f'%{gene_id}%', f'%{gene_id}%'))
                results = [dict(row) for row in cursor.fetchall()]
            
            conn.close()
            return results
        else:
            # Standard lookup for other ID types
            column = id_type
            # Add NULL/empty checks for better matching
            query = f'SELECT * FROM gene_mappings WHERE {column} = ? AND {column} IS NOT NULL AND {column} != ""'
            cursor = conn.execute(query, (gene_id,))
        
        results = [dict(row) for row in cursor.fetchall()]
        conn.close()
        
        return results
    
    def batch_lookup(self, gene_ids: List[str], id_type: str) -> Dict[str, List[Dict]]:
        """
        Perform batch lookup for multiple gene IDs.
        
        Args:
            gene_ids: List of gene identifiers
            id_type: Type of identifier
            
        Returns:
            Dictionary mapping input IDs to their records
        """
        results = {}
        
        for gene_id in gene_ids:
            results[gene_id] = self.lookup_gene(gene_id, id_type)
        
        return results
    
    def get_statistics(self) -> Dict[str, int]:
        """Get database statistics."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        stats = {}
        
        # Total records
        cursor.execute('SELECT COUNT(*) FROM gene_mappings')
        stats['total_records'] = cursor.fetchone()[0]
        
        # Records per source
        cursor.execute('SELECT source, COUNT(*) FROM gene_mappings GROUP BY source')
        stats['by_source'] = dict(cursor.fetchall())
        
        # Non-null counts for each ID type
        id_types = ['symbol', 'ensembl_gene_id', 'entrez_id', 'hgnc_id', 
                   'refseq_id', 'uniprot_id']
        
        for id_type in id_types:
            cursor.execute(f'SELECT COUNT(*) FROM gene_mappings WHERE {id_type} IS NOT NULL')
            stats[f'{id_type}_count'] = cursor.fetchone()[0]
        
        conn.close()
        return stats
