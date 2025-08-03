#!/usr/bin/env python3
"""
External Testing Script for GeneMapKit

This script compares our gene ID conversion results with MyGene.info:
1. Takes gene symbols from our results
2. Queries MyGene.info for each symbol to get Ensembl IDs
3. Compares our Ensembl IDs with MyGene.info Ensembl IDs
4. Reports number of matches vs differences for each symbol
"""

import pandas as pd
import sys
import os
from pathlib import Path
import logging
from typing import Dict, List, Tuple, Optional

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def install_mygene():
    """Install mygene if not available."""
    try:
        import mygene
        return True
    except ImportError:
        logger.info("Installing mygene...")
        import subprocess
        result = subprocess.run([sys.executable, "-m", "pip", "install", "mygene"], 
                              capture_output=True, text=True)
        if result.returncode == 0:
            logger.info("Successfully installed mygene")
            return True
        else:
            logger.error(f"Failed to install mygene: {result.stderr}")
            return False

def get_mygene_ensembl_ids(symbols: List[str], batch_size: int = 50) -> Dict[str, str]:
    """
    Get Ensembl Gene IDs for symbols from MyGene.info.
    
    Args:
        symbols: List of gene symbols to query
        batch_size: Number of symbols to query at once
        
    Returns:
        Dictionary mapping symbols to their Ensembl Gene IDs from MyGene.info
    """
    try:
        import mygene
    except ImportError:
        if not install_mygene():
            return {}
        import mygene
    
    mg = mygene.MyGeneInfo()
    results = {}
    
    # Process in batches to avoid API limits
    for i in range(0, len(symbols), batch_size):
        batch = symbols[i:i + batch_size]
        logger.info(f"Querying MyGene.info for batch {i//batch_size + 1}/{(len(symbols)-1)//batch_size + 1}")
        
        try:
            # Query MyGene.info for gene information
            gene_data = mg.querymany(
                batch, 
                scopes='symbol', 
                fields=['ensembl.gene', 'entrezgene', 'name'],
                species='human',
                returnall=True
            )
            
            # Process results
            for gene_info in gene_data['out']:
                query_symbol = gene_info.get('query', '')
                if query_symbol and 'notfound' not in gene_info:
                    # Extract Ensembl Gene ID
                    ensembl_id = None
                    if 'ensembl' in gene_info:
                        ensembl_data = gene_info['ensembl']
                        if isinstance(ensembl_data, dict):
                            ensembl_id = ensembl_data.get('gene', '')
                        elif isinstance(ensembl_data, list) and len(ensembl_data) > 0:
                            ensembl_id = ensembl_data[0].get('gene', '') if isinstance(ensembl_data[0], dict) else ''
                    
                    if ensembl_id:
                        results[query_symbol] = ensembl_id
                    
        except Exception as e:
            logger.error(f"Error querying MyGene.info for batch: {e}")
            continue
    
    return results

def compare_ensembl_ids(our_results_file: str) -> None:
    """
    Compare our Ensembl IDs with MyGene.info Ensembl IDs for each symbol.
    
    Args:
        our_results_file: Path to our conversion results CSV
    """
    # Load our results
    try:
        our_df = pd.read_csv(our_results_file)
        logger.info(f"Loaded {len(our_df)} conversion results from {our_results_file}")
    except Exception as e:
        logger.error(f"Error loading our results: {e}")
        return
    
    # Get symbols from our results
    symbols = our_df['symbol'].tolist()
    
    # Get MyGene.info Ensembl IDs
    logger.info("Getting Ensembl IDs from MyGene.info...")
    mygene_ensembl_ids = get_mygene_ensembl_ids(symbols)
    logger.info(f"Retrieved Ensembl IDs for {len(mygene_ensembl_ids)} symbols from MyGene.info")
    
    # Compare results
    comparison_results = []
    matches = 0
    differences = 0
    our_only = 0
    mygene_only = 0
    both_missing = 0
    
    print("\n" + "="*100)
    print("ENSEMBL ID COMPARISON: Our Results vs MyGene.info")
    print("="*100)
    print(f"{'Symbol':<15} {'Our Ensembl ID':<20} {'MyGene Ensembl ID':<20} {'Status':<15}")
    print("-" * 100)
    
    for _, row in our_df.iterrows():
        symbol = row.get('symbol', '')
        our_ensembl = row.get('ensembl_gene_id', '')
        mygene_ensembl = mygene_ensembl_ids.get(symbol, '')
        
        # Determine status
        if our_ensembl and mygene_ensembl:
            if our_ensembl == mygene_ensembl:
                status = "MATCH"
                matches += 1
            else:
                status = "DIFFERENT"
                differences += 1
        elif our_ensembl and not mygene_ensembl:
            status = "OUR_ONLY"
            our_only += 1
        elif not our_ensembl and mygene_ensembl:
            status = "MYGENE_ONLY"
            mygene_only += 1
        else:
            status = "BOTH_MISSING"
            both_missing += 1
        
        # Truncate long IDs for display
        our_display = our_ensembl[:18] + ".." if len(our_ensembl) > 18 else our_ensembl
        mygene_display = mygene_ensembl[:18] + ".." if len(mygene_ensembl) > 18 else mygene_ensembl
        
        print(f"{symbol:<15} {our_display:<20} {mygene_display:<20} {status:<15}")
        
        # Store detailed results
        comparison_results.append({
            'symbol': symbol,
            'our_ensembl_gene_id': our_ensembl,
            'mygene_ensembl_gene_id': mygene_ensembl,
            'status': status,
            'match': status == "MATCH"
        })
    
    # Print summary statistics
    total_symbols = len(symbols)
    print("\n" + "="*100)
    print("SUMMARY STATISTICS")
    print("="*100)
    print(f"Total symbols tested: {total_symbols}")
    print(f"Exact matches:        {matches:4d} ({matches/total_symbols*100:5.1f}%)")
    print(f"Different IDs:        {differences:4d} ({differences/total_symbols*100:5.1f}%)")
    print(f"Only in our results:  {our_only:4d} ({our_only/total_symbols*100:5.1f}%)")
    print(f"Only in MyGene.info:  {mygene_only:4d} ({mygene_only/total_symbols*100:5.1f}%)")
    print(f"Missing in both:      {both_missing:4d} ({both_missing/total_symbols*100:5.1f}%)")
    
    # Calculate accuracy where both have results
    both_have_results = matches + differences
    if both_have_results > 0:
        accuracy = matches / both_have_results * 100
        print(f"\nAccuracy (when both have results): {matches}/{both_have_results} ({accuracy:.1f}%)")
    
    # Save detailed comparison results
    results_dir = 'results'
    os.makedirs(results_dir, exist_ok=True)
    
    comparison_df = pd.DataFrame(comparison_results)
    comparison_file = f'{results_dir}/mygene_comparison.csv'
    comparison_df.to_csv(comparison_file, index=False)
    logger.info(f"Detailed comparison results saved to {comparison_file}")
    
    # Show some examples of differences
    differences_df = comparison_df[comparison_df['status'] == 'DIFFERENT']
    if not differences_df.empty:
        print(f"\nExamples of different Ensembl IDs:")
        print("-" * 80)
        for i, (_, row) in enumerate(differences_df.head(10).iterrows()):
            print(f"{row['symbol']}: Ours={row['our_ensembl_gene_id']}, MyGene={row['mygene_ensembl_gene_id']}")
    
    # Show symbols only found by us
    our_only_df = comparison_df[comparison_df['status'] == 'OUR_ONLY']
    if not our_only_df.empty:
        print(f"\nSymbols only found by our converter:")
        print("-" * 50)
        for _, row in our_only_df.head(10).iterrows():
            print(f"{row['symbol']}: {row['our_ensembl_gene_id']}")
    
    # Show symbols only found by MyGene
    mygene_only_df = comparison_df[comparison_df['status'] == 'MYGENE_ONLY']
    if not mygene_only_df.empty:
        print(f"\nSymbols only found by MyGene.info:")
        print("-" * 50)
        for _, row in mygene_only_df.head(10).iterrows():
            print(f"{row['symbol']}: {row['mygene_ensembl_gene_id']}")

def main():
    """Main function to run the comparison."""
    
    # Check if results file exists
    results_file = 'results/symbol_to_ensembl_gene_id.csv'
    if not os.path.exists(results_file):
        logger.error(f"Results file not found: {results_file}")
        logger.info("Please run the main conversion script first to generate results")
        
        # Create a sample file for testing
        sample_data = {
            'symbol': ['TP53', 'BRCA1', 'BRCA2', 'EGFR', 'MYC', 'KRAS', 'INVALID_GENE'],
            'ensembl_gene_id': ['ENSG00000141510', 'ENSG00000012048', 'ENSG00000139618', 
                              'ENSG00000146648', 'ENSG00000136997', 'ENSG00000133703', '']
        }
        sample_df = pd.DataFrame(sample_data)
        os.makedirs('results', exist_ok=True)
        sample_df.to_csv(results_file, index=False)
        logger.info(f"Created sample results file: {results_file}")
    
    # Run the comparison
    compare_ensembl_ids(results_file)
    logger.info("Comparison complete!")

if __name__ == "__main__":
    main()
