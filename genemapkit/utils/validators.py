"""
Input Validation Utilities

Provides validation functions for user inputs and file formats.
"""

import os
import pandas as pd
from pathlib import Path
from typing import List, Optional, Tuple
import logging


def validate_input_file(file_path: str, logger: Optional[logging.Logger] = None) -> bool:
    """
    Validate that input file exists and is readable.
    
    Args:
        file_path: Path to the input file
        logger: Logger instance
        
    Returns:
        True if file is valid, raises exception otherwise
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    # Check if file exists
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Input file not found: {file_path}")
    
    # Check if file is readable
    if not os.access(file_path, os.R_OK):
        raise PermissionError(f"Cannot read input file: {file_path}")
    
    # Check file extension
    valid_extensions = ['.csv', '.tsv', '.txt']
    file_ext = Path(file_path).suffix.lower()
    
    if file_ext not in valid_extensions:
        logger.warning(f"File extension '{file_ext}' not in recommended formats: {valid_extensions}")
    
    # Try to read first few lines to validate format
    try:
        # Detect delimiter
        with open(file_path, 'r') as f:
            first_line = f.readline().strip()
        
        if '\t' in first_line:
            delimiter = '\t'
        elif ',' in first_line:
            delimiter = ','
        else:
            delimiter = ','  # Default
        
        # Try to read the file
        df = pd.read_csv(file_path, delimiter=delimiter, nrows=5)
        
        if df.empty:
            raise ValueError("Input file is empty")
        
        if len(df.columns) == 0:
            raise ValueError("Input file has no columns")
        
        logger.info(f"Input file validation passed: {len(df.columns)} columns detected")
        return True
        
    except Exception as e:
        raise ValueError(f"Error reading input file: {str(e)}")


def validate_id_type(id_type: str, supported_types: List[str]) -> bool:
    """
    Validate that the ID type is supported.
    
    Args:
        id_type: Gene ID type to validate
        supported_types: List of supported ID types
        
    Returns:
        True if valid, raises exception otherwise
    """
    if id_type not in supported_types:
        raise ValueError(f"Unsupported ID type: {id_type}. Supported types: {supported_types}")
    
    return True


def validate_gene_id_format(gene_id: str, id_type: str) -> Tuple[bool, str]:
    """
    Validate gene ID format based on its type.
    
    Args:
        gene_id: Gene identifier to validate
        id_type: Type of gene identifier
        
    Returns:
        Tuple of (is_valid, message)
    """
    import re
    
    if not gene_id or not isinstance(gene_id, str):
        return False, "Gene ID must be a non-empty string"
    
    gene_id = gene_id.strip()
    
    if not gene_id:
        return False, "Gene ID cannot be empty after trimming whitespace"
    
    # Format-specific validations
    if id_type == 'symbol':
        # Gene symbols are typically 1-15 characters, alphanumeric with some special chars
        if not re.match(r'^[A-Za-z0-9@_.-]+$', gene_id):
            return False, "Gene symbol contains invalid characters"
        if len(gene_id) > 20:
            return False, "Gene symbol is too long (>20 characters)"
            
    elif id_type == 'ensembl_gene_id':
        # Ensembl gene IDs: ENSG followed by 11 digits
        if not re.match(r'^ENSG\d{11}$', gene_id):
            return False, "Ensembl gene ID should match pattern 'ENSG' + 11 digits"
            
    elif id_type == 'ensembl_transcript_id':
        # Ensembl transcript IDs: ENST followed by 11 digits
        if not re.match(r'^ENST\d{11}$', gene_id):
            return False, "Ensembl transcript ID should match pattern 'ENST' + 11 digits"
            
    elif id_type == 'entrez_id':
        # Entrez IDs are numeric
        if not re.match(r'^\d+$', gene_id):
            return False, "Entrez ID should be numeric"
            
    elif id_type == 'hgnc_id':
        # HGNC IDs: HGNC: followed by digits
        if not re.match(r'^HGNC:\d+$', gene_id):
            return False, "HGNC ID should match pattern 'HGNC:' + digits"
            
    elif id_type == 'uniprot_id':
        # UniProt IDs: 6 characters (alphanumeric)
        if not re.match(r'^[A-Z0-9]{6}$', gene_id):
            return False, "UniProt ID should be 6 alphanumeric characters"
    
    return True, "Valid format"


def validate_output_directory(output_path: str, create_if_missing: bool = True) -> bool:
    """
    Validate output directory and optionally create it.
    
    Args:
        output_path: Path to output file or directory
        create_if_missing: Whether to create directory if it doesn't exist
        
    Returns:
        True if valid/created successfully
    """
    output_dir = Path(output_path).parent
    
    if not output_dir.exists():
        if create_if_missing:
            output_dir.mkdir(parents=True, exist_ok=True)
        else:
            raise FileNotFoundError(f"Output directory does not exist: {output_dir}")
    
    if not os.access(output_dir, os.W_OK):
        raise PermissionError(f"Cannot write to output directory: {output_dir}")
    
    return True


def sanitize_filename(filename: str) -> str:
    """
    Sanitize filename by removing/replacing invalid characters.
    
    Args:
        filename: Original filename
        
    Returns:
        Sanitized filename
    """
    import re
    
    # Remove/replace invalid characters for most file systems
    sanitized = re.sub(r'[<>:"/\\|?*]', '_', filename)
    
    # Remove leading/trailing dots and spaces
    sanitized = sanitized.strip('. ')
    
    # Ensure filename is not empty
    if not sanitized:
        sanitized = 'output'
    
    return sanitized


def detect_file_format(file_path: str) -> str:
    """
    Detect file format based on content and extension.
    
    Args:
        file_path: Path to the file
        
    Returns:
        Detected format ('csv', 'tsv', or 'unknown')
    """
    file_ext = Path(file_path).suffix.lower()
    
    # Check extension first
    if file_ext in ['.csv']:
        return 'csv'
    elif file_ext in ['.tsv', '.txt']:
        return 'tsv'
    
    # Analyze content
    try:
        with open(file_path, 'r') as f:
            first_line = f.readline().strip()
        
        tab_count = first_line.count('\t')
        comma_count = first_line.count(',')
        
        if tab_count > comma_count:
            return 'tsv'
        elif comma_count > 0:
            return 'csv'
        else:
            return 'unknown'
            
    except Exception:
        return 'unknown'
