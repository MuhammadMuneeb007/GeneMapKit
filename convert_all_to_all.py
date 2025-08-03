#!/usr/bin/env python3
"""
All-to-All Gene ID Conversion Script

This script converts sample data from each ID type to all other supported ID types,
creating comprehensive conversion matrices and sample outputs.

Usage:
    python convert_all_to_all.py
"""

import sys
import logging
import pandas as pd
from pathlib import Path
from typing import List, Dict, Optional
import os

# Add project root to path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from genemapkit.core.converter import GeneConverter
from genemapkit.utils.logger import setup_logging


class AllToAllConverter:
    """Converts sample data from each ID type to all other types."""
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        """Initialize the converter."""
        self.logger = logger or logging.getLogger(__name__)
        self.converter = GeneConverter(data_dir="data/databases")
        
        # All supported ID types
        self.all_id_types = [
            'symbol', 'ensembl_gene_id', 'ensembl_transcript_id',
            'entrez_id', 'hgnc_id', 'refseq_id', 'uniprot_id'
        ]
        
        # Sample data directory
        self.sample_dir = Path("data/sample")
        self.output_dir = Path("data/conversions")
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def convert_all_sample_files(self) -> Dict[str, Dict[str, int]]:
        """
        Convert all sample files from each ID type to all other types.
        
        Returns:
            Dictionary with conversion statistics
        """
        self.logger.info("🔄 Starting all-to-all conversion of sample files...")
        
        conversion_stats = {}
        
        for input_type in self.all_id_types:
            self.logger.info(f"\n{'='*60}")
            self.logger.info(f"Converting from {input_type.upper()}")
            self.logger.info(f"{'='*60}")
            
            # Check if sample file exists
            sample_file = self.sample_dir / f"{input_type}_sample.csv"
            if not sample_file.exists():
                self.logger.warning(f"Sample file not found: {sample_file}")
                continue
            
            # Read sample data
            try:
                df = pd.read_csv(sample_file)
                gene_ids = df[input_type].dropna().tolist()
                self.logger.info(f"Loaded {len(gene_ids)} {input_type} IDs from {sample_file}")
            except Exception as e:
                self.logger.error(f"Error reading {sample_file}: {e}")
                continue
            
            # Get output types (all except current input type)
            output_types = [t for t in self.all_id_types if t != input_type]
            
            # Perform conversion
            results_df = self.converter.convert_batch(
                gene_ids, input_type, output_types, include_unmapped=True
            )
            
            # Add mapping statistics
            mapped_count = (results_df[output_types].notna().any(axis=1)).sum()
            results_df['mapping_success'] = results_df[output_types].notna().any(axis=1)
            results_df['mapped_types_count'] = results_df[output_types].notna().sum(axis=1)
            
            # Save results
            output_file = self.output_dir / f"from_{input_type}_to_all.csv"
            results_df.to_csv(output_file, index=False)
            
            # Calculate per-type mapping statistics
            type_stats = {}
            for output_type in output_types:
                mapped_for_type = results_df[output_type].notna().sum()
                type_stats[output_type] = mapped_for_type
                self.logger.info(f"  {output_type}: {mapped_for_type}/{len(gene_ids)} ({(mapped_for_type/len(gene_ids)*100):.1f}%)")
            
            conversion_stats[input_type] = {
                'total_genes': len(gene_ids),
                'overall_mapped': mapped_count,
                'mapping_rate': (mapped_count / len(gene_ids)) * 100,
                'per_type': type_stats,
                'output_file': str(output_file)
            }
            
            self.logger.info(f"✅ Saved results to: {output_file}")
            self.logger.info(f"Overall mapping: {mapped_count}/{len(gene_ids)} ({(mapped_count/len(gene_ids)*100):.1f}%)")
        
        return conversion_stats
    
    def create_comprehensive_sample(self) -> None:
        """Create a comprehensive sample file with all ID types."""
        self.logger.info("📝 Creating comprehensive sample file...")
        
        # Start with symbol sample as base
        symbol_file = self.sample_dir / "symbol_sample.csv"
        if not symbol_file.exists():
            self.logger.error("Symbol sample file not found!")
            return
        
        df = pd.read_csv(symbol_file)
        gene_symbols = df['symbol'].tolist()
        
        # Convert to all other types
        output_types = [t for t in self.all_id_types if t != 'symbol']
        results_df = self.converter.convert_batch(
            gene_symbols, 'symbol', output_types, include_unmapped=True
        )
        
        # Add additional metadata
        results_df['total_mapped_types'] = results_df[output_types].notna().sum(axis=1)
        results_df['mapping_completeness'] = (results_df['total_mapped_types'] / len(output_types) * 100).round(1)
        
        # Save comprehensive file
        comprehensive_file = self.sample_dir / "comprehensive_all_ids.csv"
        results_df.to_csv(comprehensive_file, index=False)
        
        self.logger.info(f"✅ Created comprehensive sample: {comprehensive_file}")
        
        # Print summary
        total_genes = len(results_df)
        fully_mapped = (results_df['total_mapped_types'] == len(output_types)).sum()
        partially_mapped = (results_df['total_mapped_types'] > 0).sum()
        
        self.logger.info(f"Summary:")
        self.logger.info(f"  Total genes: {total_genes}")
        self.logger.info(f"  Fully mapped (all types): {fully_mapped}")
        self.logger.info(f"  Partially mapped: {partially_mapped}")
        self.logger.info(f"  Unmapped: {total_genes - partially_mapped}")
    
    def create_conversion_matrix(self, stats: Dict[str, Dict]) -> None:
        """Create a conversion success matrix."""
        self.logger.info("📊 Creating conversion matrix...")
        
        # Create matrix data
        matrix_data = []
        
        for input_type, input_stats in stats.items():
            row = {'input_type': input_type}
            for output_type in self.all_id_types:
                if output_type == input_type:
                    row[output_type] = 100.0  # Self-mapping is always 100%
                elif output_type in input_stats['per_type']:
                    mapped_count = input_stats['per_type'][output_type]
                    total_count = input_stats['total_genes']
                    row[output_type] = round((mapped_count / total_count) * 100, 1)
                else:
                    row[output_type] = 0.0
            matrix_data.append(row)
        
        # Create DataFrame
        matrix_df = pd.DataFrame(matrix_data)
        
        # Save matrix
        matrix_file = self.output_dir / "conversion_success_matrix.csv"
        matrix_df.to_csv(matrix_file, index=False)
        
        self.logger.info(f"✅ Conversion matrix saved: {matrix_file}")
        
        # Print matrix
        print("\n🎯 Conversion Success Matrix (% successfully mapped):")
        print("=" * 80)
        print(matrix_df.to_string(index=False, float_format='%.1f'))
    
    def generate_summary_report(self, stats: Dict[str, Dict]) -> None:
        """Generate a summary report of all conversions."""
        self.logger.info("📋 Generating summary report...")
        
        report_lines = [
            "# Gene ID Conversion Summary Report",
            f"Generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
            "",
            "## Overview",
            f"- Total ID types tested: {len(stats)}",
            f"- Conversion combinations: {len(stats) * (len(self.all_id_types) - 1)}",
            "",
            "## Conversion Results by Input Type",
            ""
        ]
        
        for input_type, input_stats in stats.items():
            report_lines.extend([
                f"### {input_type.upper()}",
                f"- Total genes: {input_stats['total_genes']}",
                f"- Overall mapping success: {input_stats['mapping_rate']:.1f}%",
                f"- Output file: {input_stats['output_file']}",
                "",
                "Per-target mappings:",
            ])
            
            for output_type, count in input_stats['per_type'].items():
                percentage = (count / input_stats['total_genes']) * 100
                report_lines.append(f"- {output_type}: {count}/{input_stats['total_genes']} ({percentage:.1f}%)")
            
            report_lines.append("")
        
        # Add recommendations
        report_lines.extend([
            "## Recommendations",
            "",
            "### Best Sources for Each ID Type:",
        ])
        
        # Find best input types for each output type
        best_sources = {}
        for output_type in self.all_id_types:
            best_rate = 0
            best_source = None
            
            for input_type, input_stats in stats.items():
                if output_type in input_stats['per_type']:
                    rate = (input_stats['per_type'][output_type] / input_stats['total_genes']) * 100
                    if rate > best_rate:
                        best_rate = rate
                        best_source = input_type
            
            if best_source:
                best_sources[output_type] = (best_source, best_rate)
                report_lines.append(f"- For {output_type}: use {best_source} ({best_rate:.1f}% success)")
        
        # Save report
        report_file = self.output_dir / "conversion_summary_report.md"
        with open(report_file, 'w') as f:
            f.write('\n'.join(report_lines))
        
        self.logger.info(f"✅ Summary report saved: {report_file}")


def main():
    """Main function."""
    print("🧬 All-to-All Gene ID Conversion Tool")
    print("=" * 60)
    
    # Setup logging
    logger = setup_logging('INFO')
    
    # Initialize converter
    converter = AllToAllConverter(logger=logger)
    
    try:
        # Check if sample files exist
        missing_files = []
        for id_type in converter.all_id_types:
            sample_file = converter.sample_dir / f"{id_type}_sample.csv"
            if not sample_file.exists():
                missing_files.append(str(sample_file))
        
        if missing_files:
            logger.warning(f"Missing sample files: {missing_files}")
            logger.info("Creating missing sample files...")
            # The files should already be created by previous code
        
        # Perform all-to-all conversions
        stats = converter.convert_all_sample_files()
        
        # Create comprehensive sample
        converter.create_comprehensive_sample()
        
        # Create conversion matrix
        converter.create_conversion_matrix(stats)
        
        # Generate summary report
        converter.generate_summary_report(stats)
        
        logger.info("\n🎉 All-to-all conversion completed successfully!")
        logger.info(f"Check the '{converter.output_dir}' directory for results.")
        
    except Exception as e:
        logger.error(f"Error during conversion: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
