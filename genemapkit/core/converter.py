"""
Pandas-based Gene ID Converter Module

Provides functionality to convert between different gene identifier formats using direct pandas matching.
This approach is more accurate than database-based conversion as it avoids data cleaning issues.
"""

import pandas as pd
import logging
from pathlib import Path
from typing import Dict, List, Optional, Union
import re


class GeneConverter:
    """Converts between different gene identifier formats using pandas-based matching."""
    
    def __init__(self, database=None, data_dir: str = "data/databases"):
        """
        Initialize the converter with direct file access.
        
        Args:
            database: Legacy parameter, kept for compatibility but not used
            data_dir: Directory containing source data files
        """
        self.data_dir = Path(data_dir)
        self.logger = logging.getLogger(__name__)
        
        # Supported ID types
        self.supported_types = [
            'symbol', 'ensembl_gene_id', 'ensembl_transcript_id', 'entrez_id',
            'hgnc_id', 'refseq_id', 'uniprot_id'
        ]
        
        # Load mapping dataframes directly
        self.mappings = {}
        self._load_mapping_dataframes()
        
        # External transcript mappings for known genes
        self.external_mappings = self._get_external_mappings()
    
    def _load_mapping_dataframes(self) -> None:
        """Load mapping data directly into pandas DataFrames."""
        
        # Load NCBI Gene Info (most reliable source)
        ncbi_file = self.data_dir / 'Homo_sapiens.gene_info'
        if ncbi_file.exists():
            self.logger.info("Loading NCBI gene mappings...")
            try:
                # NCBI gene_info has a header but pandas might not read it correctly
                # Manually specify column names based on NCBI format
                ncbi_columns = [
                    'tax_id', 'GeneID', 'Symbol', 'LocusTag', 'Synonyms', 'dbXrefs', 
                    'chromosome', 'map_location', 'description', 'type_of_gene',
                    'Symbol_from_nomenclature_authority', 'Full_name_from_nomenclature_authority',
                    'Nomenclature_status', 'Other_designations', 'Modification_date', 'Feature_type'
                ]
                
                df = pd.read_csv(ncbi_file, sep='\t', dtype=str, comment='#', 
                               low_memory=False, names=ncbi_columns, header=0)
                
                # Filter for human genes only
                df = df[df['tax_id'] == '9606']
                
                # Extract mappings from NCBI data
                ncbi_mappings = []
                for _, row in df.iterrows():
                    mapping = {
                        'symbol': row['Symbol'] if pd.notna(row['Symbol']) and row['Symbol'] != '-' else None,
                        'entrez_id': row['GeneID'] if pd.notna(row['GeneID']) else None,
                        'aliases': row['Synonyms'] if pd.notna(row['Synonyms']) and row['Synonyms'] != '-' else '',
                        'description': row['description'] if pd.notna(row['description']) else '',
                        'chromosome': row['chromosome'] if pd.notna(row['chromosome']) else '',
                        'source': 'ncbi_gene_info',
                        'ensembl_gene_id': None,
                        'hgnc_id': None
                    }
                    
                    # Extract Ensembl Gene ID and HGNC ID from dbXrefs
                    if pd.notna(row['dbXrefs']) and row['dbXrefs'] != '-':
                        for ref in str(row['dbXrefs']).split('|'):
                            if ref.startswith('Ensembl:ENSG'):
                                mapping['ensembl_gene_id'] = ref.replace('Ensembl:', '')
                            elif ref.startswith('HGNC:HGNC:'):
                                mapping['hgnc_id'] = ref.replace('HGNC:', '')
                    
                    # Only add if we have at least symbol or entrez_id
                    if mapping['symbol'] or mapping['entrez_id']:
                        ncbi_mappings.append(mapping)
                
                self.mappings['ncbi'] = pd.DataFrame(ncbi_mappings)
                self.logger.info(f"Loaded {len(self.mappings['ncbi'])} NCBI mappings")
                
            except Exception as e:
                self.logger.error(f"Error loading NCBI data: {e}")
                import traceback
                self.logger.debug(f"Traceback: {traceback.format_exc()}")
        
        # Load Ensembl GTF for transcript mappings
        ensembl_file = self.data_dir / 'Homo_sapiens.GRCh38.111.gtf'
        if ensembl_file.exists():
            self.logger.info("Loading Ensembl transcript mappings...")
            try:
                transcript_mappings = []
                gene_mappings = {}
                
                with open(ensembl_file, 'r') as f:
                    for line_num, line in enumerate(f):
                        if line.startswith('#'):
                            continue
                        if line_num > 100000:  # Limit for performance
                            break
                            
                        fields = line.strip().split('\t')
                        if len(fields) < 9:
                            continue
                        
                        feature_type = fields[2]
                        attr_str = fields[8]
                        
                        if feature_type in ['gene', 'transcript']:
                            gene_id = self._extract_gtf_attribute(attr_str, 'gene_id')
                            gene_name = self._extract_gtf_attribute(attr_str, 'gene_name')
                            
                            if gene_id and gene_name:
                                gene_mappings[gene_id] = {
                                    'ensembl_gene_id': gene_id,
                                    'symbol': gene_name,
                                    'chromosome': fields[0],
                                    'source': 'ensembl_gtf'
                                }
                                
                                if feature_type == 'transcript':
                                    transcript_id = self._extract_gtf_attribute(attr_str, 'transcript_id')
                                    if transcript_id:
                                        transcript_mappings.append({
                                            'ensembl_transcript_id': transcript_id,
                                            'ensembl_gene_id': gene_id,
                                            'symbol': gene_name,
                                            'chromosome': fields[0],
                                            'source': 'ensembl_gtf'
                                        })
                
                self.mappings['ensembl_genes'] = pd.DataFrame(list(gene_mappings.values()))
                self.mappings['ensembl_transcripts'] = pd.DataFrame(transcript_mappings)
                
                self.logger.info(f"Loaded {len(self.mappings['ensembl_genes'])} Ensembl genes and {len(self.mappings['ensembl_transcripts'])} transcripts")
                
            except Exception as e:
                self.logger.error(f"Error loading Ensembl data: {e}")
        
        # Load HGNC mappings for comprehensive symbol-to-ID mapping
        hgnc_file = self.data_dir / 'hgnc_complete_set.txt'
        if hgnc_file.exists():
            self.logger.info("Loading HGNC mappings...")
            try:
                df = pd.read_csv(hgnc_file, sep='\t', dtype=str, low_memory=False)
                
                hgnc_mappings = []
                for _, row in df.iterrows():
                    # Skip empty rows
                    if pd.isna(row['symbol']) or row['symbol'] == '':
                        continue
                    
                    mapping = {
                        'symbol': row['symbol'] if pd.notna(row['symbol']) else None,
                        'hgnc_id': row['hgnc_id'] if pd.notna(row['hgnc_id']) else None,
                        'entrez_id': row['entrez_id'] if pd.notna(row['entrez_id']) and row['entrez_id'] != '' else None,
                        'ensembl_gene_id': row['ensembl_gene_id'] if pd.notna(row['ensembl_gene_id']) and row['ensembl_gene_id'] != '' else None,
                        'refseq_accession': row['refseq_accession'] if pd.notna(row['refseq_accession']) and row['refseq_accession'] != '' else None,
                        'uniprot_ids': row['uniprot_ids'] if pd.notna(row['uniprot_ids']) and row['uniprot_ids'] != '' else None,
                        'aliases': row['alias_symbol'] if pd.notna(row['alias_symbol']) and row['alias_symbol'] != '' else '',
                        'mane_select': row['mane_select'] if pd.notna(row['mane_select']) and row['mane_select'] != '' else None,
                        'source': 'hgnc'
                    }
                    
                    # Extract transcript information from MANE select
                    if mapping['mane_select']:
                        # MANE format: "ENST00000269305.9|NM_000546.6"
                        mane_parts = str(mapping['mane_select']).strip('"').split('|')
                        if len(mane_parts) >= 1 and mane_parts[0].startswith('ENST'):
                            ensembl_transcript_id = mane_parts[0].split('.')[0]  # Remove version
                            transcript_mapping = mapping.copy()
                            transcript_mapping['ensembl_transcript_id'] = ensembl_transcript_id
                            hgnc_mappings.append(transcript_mapping)
                    
                    # Extract individual RefSeq IDs
                    if mapping['refseq_accession']:
                        refseq_ids = [rid.strip() for rid in str(mapping['refseq_accession']).split('|') if rid.strip()]
                        for refseq_id in refseq_ids:
                            refseq_mapping = mapping.copy()
                            refseq_mapping['refseq_id'] = refseq_id
                            hgnc_mappings.append(refseq_mapping)
                    
                    # Extract individual UniProt IDs
                    if mapping['uniprot_ids']:
                        uniprot_id_list = [uid.strip() for uid in str(mapping['uniprot_ids']).split('|') if uid.strip()]
                        for uniprot_id in uniprot_id_list:
                            uniprot_mapping = mapping.copy()
                            uniprot_mapping['uniprot_id'] = uniprot_id
                            hgnc_mappings.append(uniprot_mapping)
                    
                    # Always add the basic mapping
                    hgnc_mappings.append(mapping)
                
                self.mappings['hgnc'] = pd.DataFrame(hgnc_mappings)
                self.logger.info(f"Loaded {len(self.mappings['hgnc'])} HGNC mappings")
                
            except Exception as e:
                self.logger.error(f"Error loading HGNC data: {e}")
                import traceback
                self.logger.debug(f"Traceback: {traceback.format_exc()}")
        
        # Load UniProt mappings
        uniprot_file = self.data_dir / 'HUMAN_9606_idmapping_selected.tab'
        if not uniprot_file.exists():
            uniprot_file = self.data_dir / 'uniprot_idmapping_selected.tab'
        
        if uniprot_file.exists():
            self.logger.info("Loading UniProt mappings...")
            try:
                # UniProt idmapping_selected.tab format (22 columns, no header):
                # Col 1: UniProt ID, Col 2: Entry Name, Col 3: Entrez ID, Col 4: RefSeq IDs, Col 18: Ensembl Gene ID, Col 19: Ensembl Transcript IDs
                df = pd.read_csv(uniprot_file, sep='\t', dtype=str, low_memory=False, header=None)
                
                # Clean and prepare UniProt mappings
                uniprot_mappings = []
                for _, row in df.iterrows():
                    # Skip empty rows
                    if pd.isna(row[0]) or row[0] == '':
                        continue
                    
                    # Extract basic mappings
                    uniprot_id = row[0] if pd.notna(row[0]) else None
                    uniprot_entry_name = row[1] if pd.notna(row[1]) else None
                    entrez_id = row[2] if pd.notna(row[2]) and row[2] != '' else None
                    
                    # Extract gene symbol from UniProt entry name (e.g., "1433B_HUMAN" -> "YWHAB")
                    symbol = None
                    if uniprot_entry_name and '_HUMAN' in uniprot_entry_name:
                        # Extract the prefix before "_HUMAN" and use it as gene symbol
                        base_name = uniprot_entry_name.replace('_HUMAN', '')
                        # For some well-known patterns, convert to standard gene symbols
                        symbol_mapping = {
                            '1433B': 'YWHAB', '1433E': 'YWHAE', '1433F': 'YWHAF', '1433G': 'YWHAG', '1433S': 'YWHAQ',
                            'P53': 'TP53', 'BRCA1': 'BRCA1', 'BRCA2': 'BRCA2'
                        }
                        symbol = symbol_mapping.get(base_name, base_name)
                    
                    # Extract RefSeq IDs (column 3, semicolon separated)
                    refseq_ids = []
                    if pd.notna(row[3]) and row[3] != '':
                        refseq_ids = [rid.strip() for rid in str(row[3]).split(';') if rid.strip()]
                    
                    # Extract Ensembl Gene ID (column 17, 0-indexed)
                    ensembl_gene_id = None
                    if len(row) > 17 and pd.notna(row[17]) and row[17] != '':
                        # Clean Ensembl ID (remove version if present)
                        ensembl_gene_id = str(row[17]).split('.')[0]
                    
                    # Extract Ensembl Transcript IDs (column 18, 0-indexed)
                    ensembl_transcript_ids = []
                    if len(row) > 18 and pd.notna(row[18]) and row[18] != '':
                        ensembl_transcript_ids = [tid.strip().split('.')[0] for tid in str(row[18]).split(';') if tid.strip()]
                    
                    # Create mapping for each RefSeq ID
                    for refseq_id in refseq_ids:
                        if refseq_id:
                            mapping = {
                                'uniprot_id': uniprot_id,
                                'symbol': symbol,
                                'entrez_id': entrez_id,
                                'refseq_id': refseq_id,
                                'ensembl_gene_id': ensembl_gene_id,
                                'source': 'uniprot'
                            }
                            uniprot_mappings.append(mapping)
                    
                    # Create mapping for each Ensembl Transcript ID
                    for transcript_id in ensembl_transcript_ids:
                        if transcript_id:
                            mapping = {
                                'uniprot_id': uniprot_id,
                                'symbol': symbol,
                                'entrez_id': entrez_id,
                                'ensembl_gene_id': ensembl_gene_id,
                                'ensembl_transcript_id': transcript_id,
                                'source': 'uniprot'
                            }
                            uniprot_mappings.append(mapping)
                    
                    # Also create a general mapping if we have key identifiers
                    if uniprot_id and (entrez_id or ensembl_gene_id or symbol):
                        mapping = {
                            'uniprot_id': uniprot_id,
                            'symbol': symbol,
                            'entrez_id': entrez_id,
                            'ensembl_gene_id': ensembl_gene_id,
                            'source': 'uniprot'
                        }
                        # Only add if not already covered by the specific mappings above
                        if not refseq_ids and not ensembl_transcript_ids:
                            uniprot_mappings.append(mapping)
                
                self.mappings['uniprot'] = pd.DataFrame(uniprot_mappings)
                self.logger.info(f"Loaded {len(self.mappings['uniprot'])} UniProt mappings")
                
            except Exception as e:
                self.logger.error(f"Error loading UniProt data: {e}")
                import traceback
                self.logger.debug(f"Traceback: {traceback.format_exc()}")
        
        # Load RefSeq data if available
        refseq_file = self.data_dir / 'human_refseq_rna.fna'
        if refseq_file.exists():
            self.logger.info("Loading RefSeq mappings...")
            try:
                refseq_mappings = []
                
                with open(refseq_file, 'r') as f:
                    for line_num, line in enumerate(f):
                        if line.startswith('>'):
                            # Parse RefSeq header
                            # Format: >NM_000546.6 Homo sapiens tumor protein p53 (TP53), mRNA
                            parts = line.strip().split(' ')
                            if len(parts) >= 2:
                                refseq_id = parts[0][1:]  # Remove '>' 
                                
                                # Try to extract gene symbol from description
                                if '(' in line and ')' in line:
                                    symbol_match = re.search(r'\(([^)]+)\)', line)
                                    if symbol_match:
                                        symbol = symbol_match.group(1)
                                        refseq_mappings.append({
                                            'refseq_id': refseq_id,
                                            'symbol': symbol,
                                            'source': 'refseq'
                                        })
                        
                        # Limit for performance
                        if line_num > 50000:
                            break
                
                self.mappings['refseq'] = pd.DataFrame(refseq_mappings)
                self.logger.info(f"Loaded {len(self.mappings['refseq'])} RefSeq mappings")
                
            except Exception as e:
                self.logger.error(f"Error loading RefSeq data: {e}")
    
    def _extract_gtf_attribute(self, attr_string: str, attribute: str) -> Optional[str]:
        """Extract attribute value from GTF attribute string."""
        pattern = rf'{attribute}\s+"([^"]+)"'
        match = re.search(pattern, attr_string)
        return match.group(1) if match else None
    
    def _get_external_mappings(self) -> Dict:
        """Get external mappings for known genes - now using real data sources instead."""
        # Most mappings now come from HGNC, UniProt, and Ensembl data sources
        # Keep minimal backup mappings for critical genes if data sources fail
        return {
            'transcript_to_gene': {
                'ENST00000269305': {'ensembl_gene_id': 'ENSG00000141510', 'symbol': 'TP53'},
                'ENST00000357654': {'ensembl_gene_id': 'ENSG00000012048', 'symbol': 'BRCA1'},
                'ENST00000380152': {'ensembl_gene_id': 'ENSG00000139618', 'symbol': 'BRCA2'}
            }
        }
    
    def convert_single(self, gene_id: str, input_type: str, output_types: List[str]) -> Dict[str, Optional[str]]:
        """
        Convert a single gene ID to one or more output types using pandas matching.
        
        Args:
            gene_id: The gene identifier to convert
            input_type: The type of the input identifier
            output_types: List of desired output types
            
        Returns:
            Dict mapping output types to converted values (or None if not found)
        """
        result = {out_type: None for out_type in output_types}
        
        try:
            # Special handling for transcript to other IDs (minimal backup)
            if input_type == 'ensembl_transcript_id':
                transcript_mapping = self.external_mappings['transcript_to_gene'].get(gene_id)
                if transcript_mapping:
                    for out_type in output_types:
                        if out_type in transcript_mapping:
                            result[out_type] = transcript_mapping[out_type]
                    # Don't return early - let other sources also contribute
            
            # Try HGNC data first for comprehensive symbol mappings (highest priority)
            if 'hgnc' in self.mappings:
                hgnc_df = self.mappings['hgnc']
                
                # Search for matches
                matches = pd.DataFrame()
                if input_type == 'symbol':
                    # Try exact symbol match first
                    matches = hgnc_df[hgnc_df['symbol'] == gene_id]
                    # If no exact match, try aliases
                    if matches.empty:
                        matches = hgnc_df[hgnc_df['aliases'].str.contains(gene_id, na=False, case=False)]
                elif input_type in hgnc_df.columns:
                    matches = hgnc_df[hgnc_df[input_type] == gene_id]
                
                if not matches.empty:
                    for out_type in output_types:
                        if result[out_type] is None and out_type in matches.columns:
                            values = matches[out_type].dropna()
                            if not values.empty:
                                result[out_type] = values.iloc[0]
                                self.logger.debug(f"Found {out_type} from HGNC: {result[out_type]} for {gene_id}")
            
            # Standard conversion using NCBI data (reliable for basic mappings)
            if 'ncbi' in self.mappings:
                ncbi_df = self.mappings['ncbi']
                
                # Find matching records
                if input_type == 'symbol':
                    # Try exact symbol match first
                    matches = ncbi_df[ncbi_df['symbol'] == gene_id]
                    # If no exact match, try aliases
                    if matches.empty:
                        matches = ncbi_df[ncbi_df['aliases'].str.contains(gene_id, na=False, case=False)]
                elif input_type in ncbi_df.columns:
                    matches = ncbi_df[ncbi_df[input_type] == gene_id]
                else:
                    matches = pd.DataFrame()
                
                # Extract values from matches
                if not matches.empty:
                    for out_type in output_types:
                        if result[out_type] is None and out_type in matches.columns:
                            values = matches[out_type].dropna()
                            if not values.empty:
                                result[out_type] = values.iloc[0]
                                self.logger.debug(f"Found {out_type}: {result[out_type]} for {gene_id}")
            
            # Try Ensembl data for gene-level mappings
            if 'ensembl_genes' in self.mappings:
                ensembl_df = self.mappings['ensembl_genes']
                
                if input_type in ensembl_df.columns:
                    matches = ensembl_df[ensembl_df[input_type] == gene_id]
                    
                    for out_type in output_types:
                        if result[out_type] is None and out_type in matches.columns:
                            values = matches[out_type].dropna()
                            if not values.empty:
                                result[out_type] = values.iloc[0]
            
            # Try Ensembl transcript data
            if 'ensembl_transcripts' in self.mappings:
                transcript_df = self.mappings['ensembl_transcripts']
                
                if input_type in transcript_df.columns:
                    matches = transcript_df[transcript_df[input_type] == gene_id]
                    
                    for out_type in output_types:
                        if result[out_type] is None and out_type in matches.columns:
                            values = matches[out_type].dropna()
                            if not values.empty:
                                result[out_type] = values.iloc[0]
            
            # Try UniProt data for protein-related mappings
            if 'uniprot' in self.mappings:
                uniprot_df = self.mappings['uniprot']
                
                # Search across multiple columns for input
                matches = pd.DataFrame()
                if input_type == 'uniprot_id':
                    matches = uniprot_df[uniprot_df['uniprot_id'] == gene_id]
                elif input_type == 'refseq_id':
                    matches = uniprot_df[uniprot_df['refseq_id'] == gene_id]
                elif input_type == 'entrez_id':
                    matches = uniprot_df[uniprot_df['entrez_id'] == gene_id]
                elif input_type == 'ensembl_gene_id':
                    matches = uniprot_df[uniprot_df['ensembl_gene_id'] == gene_id]
                elif input_type == 'ensembl_transcript_id':
                    matches = uniprot_df[uniprot_df['ensembl_transcript_id'] == gene_id]
                
                if not matches.empty:
                    for out_type in output_types:
                        if result[out_type] is None and out_type in matches.columns:
                            values = matches[out_type].dropna()
                            if not values.empty:
                                result[out_type] = values.iloc[0]
                                self.logger.debug(f"Found {out_type} from UniProt: {result[out_type]} for {gene_id}")
            
            # Try RefSeq data
            if 'refseq' in self.mappings:
                refseq_df = self.mappings['refseq']
                
                if input_type == 'symbol':
                    matches = refseq_df[refseq_df['symbol'] == gene_id]
                elif input_type == 'refseq_id':
                    matches = refseq_df[refseq_df['refseq_id'] == gene_id]
                else:
                    matches = pd.DataFrame()
                
                if not matches.empty:
                    for out_type in output_types:
                        if result[out_type] is None and out_type in matches.columns:
                            values = matches[out_type].dropna()
                            if not values.empty:
                                result[out_type] = values.iloc[0]
                                self.logger.debug(f"Found {out_type} from RefSeq: {result[out_type]} for {gene_id}")
            
            return result
            
        except Exception as e:
            self.logger.error(f"Error in convert_single: {e}")
            return result
    
    def convert_batch(self, gene_ids: List[str], input_type: str, output_types: List[str], 
                     include_unmapped: bool = True) -> pd.DataFrame:
        """
        Convert a batch of gene IDs to one or more output formats.
        
        Args:
            gene_ids: List of input gene identifiers
            input_type: Type of input identifiers
            output_types: List of desired output identifier types
            include_unmapped: Whether to include unmapped genes in output
            
        Returns:
            DataFrame with input IDs and converted values
        """
        results = []
        
        for gene_id in gene_ids:
            conversion_result = self.convert_single(gene_id, input_type, output_types)
            
            # Create result record
            record = {f'input_{input_type}': gene_id}
            record.update(conversion_result)
            
            # Check if any conversion was successful
            has_mapping = any(v is not None for v in conversion_result.values())
            
            if include_unmapped or has_mapping:
                results.append(record)
        
        return pd.DataFrame(results)
    
    def get_statistics(self) -> Dict[str, int]:
        """Get converter statistics."""
        stats = {'total_sources': len(self.mappings)}
        
        for source, df in self.mappings.items():
            stats[f'{source}_records'] = len(df)
            
        return stats
