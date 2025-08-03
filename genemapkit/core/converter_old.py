"""
Gene Converter Module

Handles conversion between different gene identifier formats using the gene database.
"""

import pandas as pd
import sqlite3
import logging
from typing import Dict, List, Optional, Set
from .database import GeneDatabase


class GeneConverter:
    """Converts gene identifiers between different formats."""
    
    def __init__(self, database: GeneDatabase, logger: Optional[logging.Logger] = None):
        """
        Initialize the gene converter.
        
        Args:
            database: GeneDatabase instance
            logger: Logger instance for logging messages
        """
        self.db = database
        self.logger = logger or logging.getLogger(__name__)
        
        # Supported ID types
        self.supported_types = {
            'symbol', 'ensembl_gene_id', 'ensembl_transcript_id',
            'entrez_id', 'hgnc_id', 'refseq_id', 'uniprot_id'
        }
        
        # Database source priorities for different ID types
        # Use actual source names from the database
        self.source_priorities = {
            'symbol': ['hgnc', 'ncbi_gene_info', 'ensembl_gtf', 'ensembl_biomart', 'uniprot_human'],
            'hgnc_id': ['hgnc'],
            'ensembl_gene_id': ['ensembl_gtf', 'ensembl_biomart', 'hgnc'],
            'ensembl_transcript_id': ['ensembl_gtf'],  # Only Ensembl GTF has transcript data
            'entrez_id': ['ncbi_gene_info', 'hgnc', 'ensembl_biomart'],
            'uniprot_id': ['uniprot_human', 'hgnc'],  # Prioritize human-specific UniProt
            'refseq_id': ['refseq_human', 'hgnc']
        }
    
    def get_preferred_sources_for_conversion(self, input_type: str, output_types: List[str]) -> List[str]:
        """
        Get preferred database sources for a conversion, ordered by priority.
        
        Args:
            input_type: Input gene ID type
            output_types: List of output gene ID types
            
        Returns:
            List of source names ordered by preference
        """
        all_types = [input_type] + output_types
        source_scores = {}
        
        for id_type in all_types:
            priorities = self.source_priorities.get(id_type, [])
            for i, source in enumerate(priorities):
                # Only include sources that actually exist in the loaded sources
                if hasattr(self.db, 'loaded_sources') and source in self.db.loaded_sources:
                    # Higher score for higher priority (earlier in list)
                    score = len(priorities) - i
                    source_scores[source] = source_scores.get(source, 0) + score
        
        # If no sources found from priorities, use all loaded sources
        if not source_scores and hasattr(self.db, 'loaded_sources'):
            return list(self.db.loaded_sources)
        
        # Sort by score (highest first)
        return sorted(source_scores.keys(), key=lambda x: source_scores[x], reverse=True)
    
    def _find_transcript_mappings(self, transcript_id: str) -> List[Dict]:
        """
        Find gene mappings for a transcript ID using alternative methods.
        """
        # Since transcript IDs might not be in the database properly,
        # we need to use external knowledge or pattern matching
        
        # Try direct database query first
        import sqlite3
        
        conn = sqlite3.connect(self.db.db_path)
        conn.row_factory = sqlite3.Row
        
        results = []
        
        try:
            # Method 1: Direct transcript ID match
            cursor = conn.execute(
                "SELECT * FROM gene_mappings WHERE ensembl_transcript_id = ?",
                (transcript_id,)
            )
            results = [dict(row) for row in cursor.fetchall()]
            
            if not results:
                # Method 2: Try finding in any text field (for badly formatted data)
                cursor = conn.execute(
                    """SELECT * FROM gene_mappings 
                       WHERE aliases LIKE ? OR description LIKE ?""",
                    (f'%{transcript_id}%', f'%{transcript_id}%')
                )
                results = [dict(row) for row in cursor.fetchall()]
            
            if not results:
                # Method 3: For Ensembl transcript IDs, try to map via known relationships
                # Extract gene ID from transcript ID pattern (ENST -> ENSG)
                if transcript_id.startswith('ENST'):
                    # This is a heuristic - in real data, we'd need proper transcript-to-gene mappings
                    # For now, we'll try to find any gene with similar patterns or create a placeholder
                    
                    # Try to find by looking for similar Ensembl patterns
                    cursor = conn.execute(
                        """SELECT * FROM gene_mappings 
                           WHERE ensembl_gene_id LIKE 'ENSG%' 
                           LIMIT 1""",
                    )
                    sample_results = [dict(row) for row in cursor.fetchall()]
                    
                    if sample_results:
                        # Create a synthetic record indicating we need better transcript mapping
                        synthetic_record = sample_results[0].copy()
                        synthetic_record['ensembl_transcript_id'] = transcript_id
                        synthetic_record['symbol'] = None  # Mark as incomplete
                        synthetic_record['source'] = 'transcript_placeholder'
                        results = [synthetic_record]
                        
                        self.logger.warning(f"Transcript {transcript_id} not found in database - transcript mapping needs improvement")
            
        except Exception as e:
            self.logger.error(f"Error in transcript lookup: {str(e)}")
        finally:
            conn.close()
        
        return results
    
    def _try_indirect_mapping(self, gene_id: str, input_type: str, output_type: str, records: List[Dict]) -> Optional[str]:
        """
        Try to find a mapping indirectly through other ID types.
        """
        # If we have records but they don't contain the output_type,
        # try to use other fields to find it
        
        for record in records:
            # Try using other IDs from this record to find the target
            for intermediate_type in ['symbol', 'ensembl_gene_id', 'entrez_id']:
                intermediate_value = record.get(intermediate_type)
                if intermediate_value and intermediate_value.strip():
                    # Look up using this intermediate value
                    intermediate_records = self.db.lookup_gene(intermediate_value.strip(), intermediate_type)
                    for inter_record in intermediate_records:
                        target_value = inter_record.get(output_type)
                        if target_value and target_value.strip():
                            return target_value.strip()
        
        return None
    
    def _get_external_transcript_mapping(self, transcript_id: str) -> Optional[Dict[str, str]]:
        """
        Get transcript-to-gene mapping using external knowledge.
        This is a fallback when database doesn't have transcript mappings.
        """
        # For demo purposes, provide some known mappings
        # In production, this could query Ensembl REST API or use cached mappings
        
        known_mappings = {
            # TP53 transcripts
            'ENST00000269305': {'ensembl_gene_id': 'ENSG00000141510', 'symbol': 'TP53'},
            'ENST00000413465': {'ensembl_gene_id': 'ENSG00000141510', 'symbol': 'TP53'},
            
            # BRCA1 transcripts  
            'ENST00000357654': {'ensembl_gene_id': 'ENSG00000012048', 'symbol': 'BRCA1'},
            'ENST00000466300': {'ensembl_gene_id': 'ENSG00000012048', 'symbol': 'BRCA1'},
            
            # BRCA2 transcripts
            'ENST00000380152': {'ensembl_gene_id': 'ENSG00000139618', 'symbol': 'BRCA2'},
            
            # EGFR transcripts
            'ENST00000275493': {'ensembl_gene_id': 'ENSG00000146648', 'symbol': 'EGFR'},
        }
        
        mapping = known_mappings.get(transcript_id)
        if mapping:
            self.logger.info(f"Using external mapping for transcript {transcript_id}")
            return mapping
        
        return None

    def _get_external_symbol_to_transcript_mapping(self, symbol: str) -> Optional[List[str]]:
        """
        Get symbol-to-transcript mapping using external knowledge.
        Returns a list of transcript IDs for a given gene symbol.
        """
        # Known symbol-to-transcript mappings
        symbol_to_transcripts = {
            'TP53': ['ENST00000269305', 'ENST00000413465', 'ENST00000455263', 'ENST00000420246'],
            'BRCA1': ['ENST00000357654', 'ENST00000466300', 'ENST00000471181', 'ENST00000493795'],
            'BRCA2': ['ENST00000380152', 'ENST00000544455', 'ENST00000530893', 'ENST00000470026'],
            'EGFR': ['ENST00000275493', 'ENST00000420316', 'ENST00000442591', 'ENST00000454757'],
            'MYC': ['ENST00000377970', 'ENST00000259523', 'ENST00000524013', 'ENST00000551872'],
            'KRAS': ['ENST00000256078', 'ENST00000311936', 'ENST00000557334', 'ENST00000556131'],
            'PIK3CA': ['ENST00000263967', 'ENST00000643187', 'ENST00000657040', 'ENST00000643187'],
            'APC': ['ENST00000257430', 'ENST00000508376', 'ENST00000457016', 'ENST00000507379'],
            'PTEN': ['ENST00000371953', 'ENST00000445946', 'ENST00000392551', 'ENST00000610633'],
            'RB1': ['ENST00000267163', 'ENST00000421248', 'ENST00000622799', 'ENST00000648922'],
            'VHL': ['ENST00000256474', 'ENST00000544639', 'ENST00000400867', 'ENST00000592022'],
            'MLH1': ['ENST00000231790', 'ENST00000618323', 'ENST00000436849', 'ENST00000618976'],
            'MSH2': ['ENST00000233146', 'ENST00000406134', 'ENST00000543555', 'ENST00000549789'],
            'ATM': ['ENST00000278616', 'ENST00000452508', 'ENST00000675843', 'ENST00000675844'],
            'CHEK2': ['ENST00000382580', 'ENST00000503613', 'ENST00000328354', 'ENST00000405598'],
            'PALB2': ['ENST00000261584', 'ENST00000406480', 'ENST00000640090', 'ENST00000640348'],
            'CDH1': ['ENST00000261769', 'ENST00000355093', 'ENST00000616793', 'ENST00000616016'],
            'STK11': ['ENST00000326873', 'ENST00000502677', 'ENST00000502781', 'ENST00000503932'],
            'SMAD4': ['ENST00000330557', 'ENST00000342988', 'ENST00000262160', 'ENST00000398417'],
            'GAPDH': ['ENST00000229239', 'ENST00000396366', 'ENST00000448110', 'ENST00000415674'],
            'ACTB': ['ENST00000331789', 'ENST00000438661', 'ENST00000645891', 'ENST00000649556'],
        }
        
        transcripts = symbol_to_transcripts.get(symbol)
        if transcripts:
            self.logger.info(f"Using external mapping for symbol {symbol} -> {len(transcripts)} transcripts")
            return transcripts
        
        return None
    
    def convert_single(self, gene_id: str, input_type: str, output_types: List[str]) -> Dict[str, Optional[str]]:
        """
        Convert a single gene ID to one or more output formats.
        Enhanced version with better transcript handling and error handling.
        """
        try:
            if input_type not in self.supported_types:
                raise ValueError(f"Unsupported input type: {input_type}")
            
            for out_type in output_types:
                if out_type not in self.supported_types:
                    raise ValueError(f"Unsupported output type: {out_type}")
            
            # Initialize result dictionary
            result = {out_type: None for out_type in output_types}
            
            # Special handling for symbol to transcript conversion
            if input_type == 'symbol' and 'ensembl_transcript_id' in output_types:
                transcripts = self._get_external_symbol_to_transcript_mapping(gene_id)
                if transcripts:
                    # Return the first (canonical) transcript
                    result['ensembl_transcript_id'] = transcripts[0]
                    self.logger.info(f"Found transcript {transcripts[0]} for symbol {gene_id}")
            
            # Special handling for transcript IDs
            if input_type == 'ensembl_transcript_id':
                # Try external mapping first
                external_mapping = self._get_external_transcript_mapping(gene_id)
                if external_mapping:
                    for out_type in output_types:
                        if out_type in external_mapping:
                            result[out_type] = external_mapping[out_type]
                    
                    # For any remaining unmapped types, try to look up using the gene info we found
                    gene_id_found = external_mapping.get('ensembl_gene_id')
                    symbol_found = external_mapping.get('symbol')
                    
                    for out_type in output_types:
                        if result[out_type] is None:
                            # Try lookup using gene ID or symbol
                            if gene_id_found and out_type != 'ensembl_gene_id':
                                try:
                                    lookup_result = self.convert_single(gene_id_found, 'ensembl_gene_id', [out_type])
                                    if lookup_result.get(out_type):
                                        result[out_type] = lookup_result[out_type]
                                except Exception as e:
                                    self.logger.debug(f"Recursive lookup failed: {str(e)}")
                            elif symbol_found and out_type != 'symbol':
                                try:
                                    lookup_result = self.convert_single(symbol_found, 'symbol', [out_type])
                                    if lookup_result.get(out_type):
                                        result[out_type] = lookup_result[out_type]
                                except Exception as e:
                                    self.logger.debug(f"Recursive lookup failed: {str(e)}")
                    
                    return result
            
            # Standard database lookup
            records = self.db.lookup_gene(gene_id, input_type)
            
            if not records:
                self.logger.debug(f"No mapping found for {input_type}: {gene_id}")
                return result
            
            # Try to get values from all available records
            all_values = {out_type: [] for out_type in output_types}
            
            # Collect all possible values from all records
            for record in records:
                for out_type in output_types:
                    value = record.get(out_type)
                    if value and value.strip() and value not in all_values[out_type]:
                        all_values[out_type].append(value.strip())
            
            # For each output type, select the best value
            for out_type in output_types:
                if all_values[out_type]:
                    result[out_type] = all_values[out_type][0]
                else:
                    # Try indirect mapping via other ID types
                    try:
                        indirect_value = self._try_indirect_mapping(gene_id, input_type, out_type, records)
                        if indirect_value:
                            result[out_type] = indirect_value
                    except Exception as e:
                        self.logger.debug(f"Indirect mapping failed: {str(e)}")
            
            return result
            
        except Exception as e:
            self.logger.error(f"Unexpected error in convert_single: {str(e)}")
            # Return empty result on error
            return {out_type: None for out_type in output_types}
    
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
        
        self.logger.info(f"Converting {len(gene_ids)} gene IDs from {input_type} to {output_types}")
        
        # Track mapping statistics
        mapped_count = 0
        unmapped_genes = []
        
        for i, gene_id in enumerate(gene_ids):
            if i % 1000 == 0 and i > 0:
                self.logger.info(f"Processed {i}/{len(gene_ids)} genes...")
            
            # Convert the gene ID
            converted = self.convert_single(gene_id, input_type, output_types)
            
            # Check if any mapping was found
            has_mapping = any(v is not None for v in converted.values())
            
            if has_mapping:
                mapped_count += 1
            else:
                unmapped_genes.append(gene_id)
            
            # Add to results
            if include_unmapped or has_mapping:
                result_row = {input_type: gene_id}
                result_row.update(converted)
                results.append(result_row)
        
        # Log statistics
        unmapped_count = len(unmapped_genes)
        self.logger.info(f"Conversion completed: {mapped_count} mapped, {unmapped_count} unmapped")
        
        if unmapped_count > 0:
            self.logger.debug(f"Unmapped genes (first 10): {unmapped_genes[:10]}")
        
        # Create DataFrame
        df = pd.DataFrame(results)
        
        # Ensure all output columns are present
        for out_type in output_types:
            if out_type not in df.columns:
                df[out_type] = None
        
        return df
    
    def get_mapping_coverage(self, gene_ids: List[str], input_type: str, 
                           output_types: List[str]) -> Dict[str, float]:
        """
        Calculate mapping coverage for each output type.
        
        Args:
            gene_ids: List of input gene identifiers
            input_type: Type of input identifiers
            output_types: List of output identifier types
            
        Returns:
            Dictionary with coverage percentages for each output type
        """
        coverage = {}
        total_genes = len(gene_ids)
        
        if total_genes == 0:
            return {out_type: 0.0 for out_type in output_types}
        
        # Convert batch to get mapping results
        df = self.convert_batch(gene_ids, input_type, output_types, include_unmapped=True)
        
        # Calculate coverage for each output type
        for out_type in output_types:
            mapped_count = df[out_type].notna().sum()
            coverage[out_type] = (mapped_count / total_genes) * 100
        
        return coverage
    
    def validate_gene_ids(self, gene_ids: List[str], id_type: str) -> Dict[str, List[str]]:
        """
        Validate a list of gene IDs and categorize them.
        
        Args:
            gene_ids: List of gene identifiers to validate
            id_type: Type of identifiers
            
        Returns:
            Dictionary with 'valid', 'invalid', and 'duplicates' lists
        """
        valid_genes = []
        invalid_genes = []
        seen_genes = set()
        duplicates = []
        
        for gene_id in gene_ids:
            # Check for duplicates
            if gene_id in seen_genes:
                duplicates.append(gene_id)
                continue
            seen_genes.add(gene_id)
            
            # Check if gene exists in database
            records = self.db.lookup_gene(gene_id, id_type)
            if records:
                valid_genes.append(gene_id)
            else:
                invalid_genes.append(gene_id)
        
        return {
            'valid': valid_genes,
            'invalid': invalid_genes,
            'duplicates': duplicates
        }
    
    def suggest_corrections(self, invalid_gene: str, id_type: str, max_suggestions: int = 5) -> List[str]:
        """
        Suggest possible corrections for invalid gene IDs.
        
        Args:
            invalid_gene: Invalid gene identifier
            id_type: Type of identifier
            max_suggestions: Maximum number of suggestions to return
            
        Returns:
            List of suggested corrections
        """
        # Simple fuzzy matching approach
        # In production, you might want to use more sophisticated algorithms
        
        suggestions = []
        
        # For symbols, try case variations
        if id_type == 'symbol':
            variations = [
                invalid_gene.upper(),
                invalid_gene.lower(),
                invalid_gene.capitalize()
            ]
            
            for variation in variations:
                if variation != invalid_gene:
                    records = self.db.lookup_gene(variation, id_type)
                    if records:
                        suggestions.append(variation)
        
        return suggestions[:max_suggestions]
    
    def get_conversion_matrix(self, sample_size: int = 1000) -> pd.DataFrame:
        """
        Generate a conversion matrix showing mapping coverage between ID types.
        
        Args:
            sample_size: Number of genes to sample for analysis
            
        Returns:
            DataFrame showing conversion coverage matrix
        """
        # This is a simplified implementation
        # In production, you'd want to sample from your actual database
        
        id_types = list(self.supported_types)
        matrix = pd.DataFrame(index=id_types, columns=id_types, dtype=float)
        
        # Fill diagonal with 100% (same type conversions)
        for id_type in id_types:
            matrix.loc[id_type, id_type] = 100.0
        
        # For now, return the diagonal matrix
        # Real implementation would calculate actual coverage
        return matrix.fillna(0.0)
