#!/usr/bin/env python3
"""
GeneMapKit - Entry Point

This is the main entry point for the GeneMapKit command-line tool.
Enhanced to convert one input format to all other formats with visualization.
"""

import sys
import csv
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
import logging

# Optional visualization imports
try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    VISUALIZATION_AVAILABLE = True
except ImportError:
    VISUALIZATION_AVAILABLE = False

# Add the project root to Python path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from genemapkit.core.converter import GeneConverter
from genemapkit.utils.logger import setup_logging

def convert_to_all_formats(input_file, input_column, results_dir="results"):
    """
    Convert input gene identifiers to all other supported formats.
    
    Args:
        input_file: Path to input CSV file
        input_column: Column name containing gene identifiers
        results_dir: Directory to save results
    """
    
    # Setup logging
    logger = setup_logging(level='INFO')
    
    # Create results directory
    results_path = Path(results_dir)
    results_path.mkdir(exist_ok=True)
    
    # Initialize pandas-based converter (no database needed)
    logger.info("Initializing GeneMapKit pandas-based converter...")
    converter = GeneConverter(data_dir="data/databases")
    
    # Read input file
    logger.info(f"Reading input file: {input_file}")
    try:
        df = pd.read_csv(input_file)
        if input_column not in df.columns:
            logger.error(f"Column '{input_column}' not found in input file. Available columns: {list(df.columns)}")
            return
    except Exception as e:
        logger.error(f"Error reading input file: {e}")
        return
    
    # Detect input type based on column name or content
    input_type = detect_input_type(input_column, df[input_column].iloc[0] if len(df) > 0 else "")
    logger.info(f"Detected input type: {input_type}")
    
    # Define all possible output types
    all_types = ['symbol', 'ensembl_gene_id', 'ensembl_transcript_id', 'entrez_id', 'hgnc_id', 'refseq_id', 'uniprot_id']
    output_types = [t for t in all_types if t != input_type]
    
    logger.info(f"Converting {len(df)} identifiers from {input_type} to: {', '.join(output_types)}")
    
    # Track conversion statistics
    conversion_stats = {}
    
    # Convert to each output type
    for output_type in output_types:
        logger.info(f"Converting to {output_type}...")
        
        try:
            # Perform conversion
            results = []
            successful_conversions = 0
            
            for idx, input_id in enumerate(df[input_column]):
                if pd.isna(input_id) or str(input_id).strip() == "":
                    results.append("")
                    continue
                
                converted = converter.convert_single(str(input_id).strip(), input_type, [output_type])
                result_value = converted.get(output_type) if isinstance(converted, dict) else converted
                if result_value and result_value != "":
                    results.append(result_value)
                    successful_conversions += 1
                else:
                    results.append("")
            
            # Create output dataframe
            output_df = df.copy()
            output_df[output_type] = results
            
            # Save results
            output_file = results_path / f"{input_type}_to_{output_type}.csv"
            output_df.to_csv(output_file, index=False)
            
            # Calculate success rate
            total_non_empty = len([x for x in df[input_column] if not pd.isna(x) and str(x).strip() != ""])
            success_rate = (successful_conversions / total_non_empty * 100) if total_non_empty > 0 else 0
            
            conversion_stats[output_type] = {
                'total_input': total_non_empty,
                'successful': successful_conversions,
                'success_rate': success_rate,
                'output_file': str(output_file)
            }
            
            logger.info(f"  ✓ {output_type}: {successful_conversions}/{total_non_empty} ({success_rate:.1f}%) - saved to {output_file}")
            
        except Exception as e:
            logger.error(f"  ✗ Error converting to {output_type}: {e}")
            conversion_stats[output_type] = {
                'total_input': len(df),
                'successful': 0,
                'success_rate': 0,
                'error': str(e)
            }
    
    # Create conversion summary
    create_conversion_summary(results_path, input_file, input_type, conversion_stats)
    
    # Create visualization diagrams if matplotlib is available
    if VISUALIZATION_AVAILABLE:
        create_visualization_diagrams(results_path, input_type, conversion_stats)
    else:
        logger.info("Matplotlib not available - skipping visualization diagrams")
        logger.info("Install matplotlib with: pip install matplotlib seaborn")
    
    logger.info(f"\nConversion complete! Results saved in: {results_path}")
    logger.info("Check conversion_summary.txt for detailed statistics.")
    if VISUALIZATION_AVAILABLE:
        logger.info("Check ConversionSuccessDiagram.png for visual analysis.")

def detect_input_type(column_name, sample_value):
    """
    Detect the input gene identifier type based on column name and sample value.
    """
    column_lower = column_name.lower()
    sample_str = str(sample_value).upper()
    
    # Check column name patterns
    if 'symbol' in column_lower or column_lower in ['gene', 'gene_name', 'gene_symbol']:
        return 'symbol'
    elif 'ensembl_gene' in column_lower or 'ensg' in column_lower:
        return 'ensembl_gene_id'
    elif 'ensembl_transcript' in column_lower or 'enst' in column_lower:
        return 'ensembl_transcript_id'
    elif 'entrez' in column_lower or 'ncbi' in column_lower:
        return 'entrez_id'
    elif 'hgnc' in column_lower:
        return 'hgnc_id'
    elif 'refseq' in column_lower or 'nm_' in column_lower:
        return 'refseq_id'
    elif 'uniprot' in column_lower:
        return 'uniprot_id'
    
    # Check sample value patterns
    if sample_str.startswith('ENSG'):
        return 'ensembl_gene_id'
    elif sample_str.startswith('ENST'):
        return 'ensembl_transcript_id'
    elif sample_str.startswith('HGNC:'):
        return 'hgnc_id'
    elif sample_str.startswith(('NM_', 'NR_', 'XM_', 'XR_')):
        return 'refseq_id'
    elif sample_str.replace('.', '').isdigit():
        return 'entrez_id'
    elif len(sample_str) == 6 and sample_str[0].isalpha():
        return 'uniprot_id'
    
    # Default to symbol if uncertain
    return 'symbol'

def create_conversion_summary(results_path, input_file, input_type, conversion_stats):
    """
    Create a summary file with conversion statistics.
    """
    summary_file = results_path / "conversion_summary.txt"
    
    with open(summary_file, 'w') as f:
        f.write("GeneMapKit Conversion Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Input file: {input_file}\n")
        f.write(f"Input type: {input_type}\n\n")
        
        f.write("Conversion Results:\n")
        f.write("-" * 30 + "\n")
        
        total_success = 0
        total_attempted = 0
        
        for output_type, stats in conversion_stats.items():
            if 'error' in stats:
                f.write(f"{output_type:20}: ERROR - {stats['error']}\n")
            else:
                f.write(f"{output_type:20}: {stats['successful']:4d}/{stats['total_input']:4d} ({stats['success_rate']:5.1f}%)\n")
                total_success += stats['successful']
                total_attempted += stats['total_input']
        
        f.write("-" * 30 + "\n")
        overall_rate = (total_success / total_attempted * 100) if total_attempted > 0 else 0
        f.write(f"{'Overall':20}: {total_success:4d}/{total_attempted:4d} ({overall_rate:5.1f}%)\n\n")
        
        f.write("Output Files:\n")
        f.write("-" * 30 + "\n")
        for output_type, stats in conversion_stats.items():
            if 'output_file' in stats:
                f.write(f"{output_type}: {stats['output_file']}\n")

def create_visualization_diagrams(results_path, input_type, conversion_stats):
    """
    Create visualization diagrams for conversion success rates.
    """
    try:
        # Prepare data for visualization
        output_types = list(conversion_stats.keys())
        success_rates = [stats['success_rate'] for stats in conversion_stats.values()]
        successful_counts = [stats['successful'] for stats in conversion_stats.values()]
        total_counts = [stats['total_input'] for stats in conversion_stats.values()]
        
        # Create figure with subplots
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle(f'Gene ID Conversion Analysis: {input_type.upper()} to All Types', fontsize=16, fontweight='bold')
        
        # 1. Success Rate Bar Chart
        colors = plt.cm.RdYlGn([rate/100 for rate in success_rates])
        bars = ax1.bar(output_types, success_rates, color=colors)
        ax1.set_title('Conversion Success Rates', fontweight='bold')
        ax1.set_ylabel('Success Rate (%)')
        ax1.set_ylim(0, 100)
        
        # Add value labels on bars
        for bar, rate in zip(bars, success_rates):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height + 1,
                    f'{rate:.1f}%', ha='center', va='bottom', fontweight='bold')
        
        # Rotate x-axis labels
        ax1.tick_params(axis='x', rotation=45)
        
        # 2. Absolute Numbers Bar Chart
        ax2.bar(output_types, successful_counts, alpha=0.7, color='skyblue', label='Successful')
        ax2.bar(output_types, [total - success for total, success in zip(total_counts, successful_counts)], 
                bottom=successful_counts, alpha=0.7, color='lightcoral', label='Failed')
        ax2.set_title('Conversion Counts (Successful vs Failed)', fontweight='bold')
        ax2.set_ylabel('Number of Genes')
        ax2.legend()
        ax2.tick_params(axis='x', rotation=45)
        
        # 3. Success Rate Pie Chart (Overall Summary)
        total_success = sum(successful_counts)
        total_attempts = sum(total_counts)
        total_failure = total_attempts - total_success
        
        if total_attempts > 0:
            overall_rate = (total_success / total_attempts) * 100
            pie_data = [total_success, total_failure]
            pie_labels = [f'Successful\n({overall_rate:.1f}%)', f'Failed\n({100-overall_rate:.1f}%)']
            colors_pie = ['lightgreen', 'lightcoral']
            
            ax3.pie(pie_data, labels=pie_labels, colors=colors_pie, autopct='%1.0f', startangle=90)
            ax3.set_title('Overall Conversion Success', fontweight='bold')
        
        # 4. Performance Heatmap Style
        # Create a simple matrix showing the performance
        matrix_data = np.array(success_rates).reshape(1, -1)
        im = ax4.imshow(matrix_data, cmap='RdYlGn', aspect='auto', vmin=0, vmax=100)
        ax4.set_title('Conversion Success Heatmap', fontweight='bold')
        ax4.set_xticks(range(len(output_types)))
        ax4.set_xticklabels(output_types, rotation=45)
        ax4.set_yticks([0])
        ax4.set_yticklabels([input_type.upper()])
        
        # Add text annotations
        for i, rate in enumerate(success_rates):
            ax4.text(i, 0, f'{rate:.1f}%', ha='center', va='center', 
                    color='white' if rate < 50 else 'black', fontweight='bold')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax4, orientation='horizontal', pad=0.1)
        cbar.set_label('Success Rate (%)')
        
        plt.tight_layout()
        
        # Save the plot
        diagram_file = results_path / "ConversionSuccessDiagram.png"
        plt.savefig(diagram_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        # Create a detailed text-based diagram as well
        create_text_diagram(results_path, input_type, conversion_stats)
        
        return str(diagram_file)
        
    except Exception as e:
        print(f"Error creating visualization: {e}")
        return None

def create_text_diagram(results_path, input_type, conversion_stats):
    """
    Create a text-based conversion success diagram.
    """
    diagram_file = results_path / "ConversionSuccessTextDiagram.txt"
    
    with open(diagram_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write(f"Gene ID Conversion Success Analysis\n")
        f.write(f"Input Type: {input_type.upper()}\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("=" * 80 + "\n\n")
        
        f.write("Conversion Results Summary:\n")
        f.write("-" * 40 + "\n")
        
        for output_type, stats in conversion_stats.items():
            success_rate = stats['success_rate']
            successful = stats['successful']
            total = stats['total_input']
            
            # Create visual bar
            bar_length = 20
            filled_length = int(bar_length * success_rate / 100)
            bar = "█" * filled_length + "░" * (bar_length - filled_length)
            
            # Status emoji
            if success_rate >= 90:
                emoji = "🟢"
            elif success_rate >= 70:
                emoji = "🟡"
            elif success_rate >= 50:
                emoji = "🟠"
            else:
                emoji = "🔴"
            
            f.write(f"{output_type:20} {emoji} [{bar}] {success_rate:5.1f}% ({successful:3d}/{total:3d})\n")
        
        f.write("\n" + "-" * 40 + "\n")
        
        # Overall statistics
        total_success = sum(stats['successful'] for stats in conversion_stats.values())
        total_attempts = sum(stats['total_input'] for stats in conversion_stats.values())
        overall_rate = (total_success / total_attempts * 100) if total_attempts > 0 else 0
        
        f.write(f"Overall Success Rate: {overall_rate:.1f}% ({total_success}/{total_attempts})\n")
        
        # Best and worst performers
        if conversion_stats:
            best = max(conversion_stats.items(), key=lambda x: x[1]['success_rate'])
            worst = min(conversion_stats.items(), key=lambda x: x[1]['success_rate'])
            
            f.write(f"Best Performance:  {best[0]} ({best[1]['success_rate']:.1f}%)\n")
            f.write(f"Worst Performance: {worst[0]} ({worst[1]['success_rate']:.1f}%)\n")
        
        f.write("\n" + "=" * 80 + "\n")
        f.write("Legend: 🟢 Excellent (≥90%) 🟡 Good (≥70%) 🟠 Fair (≥50%) 🔴 Poor (<50%)\n")
        f.write("Bar: ████████████████████ = 100% success rate\n")

def main():
    """
    Main function with command line interface.
    """
    if len(sys.argv) != 3:
        print("Usage: python main.py <input_file> <input_column>")
        print("\nExample:")
        print("  python main.py data/sample/gene_symbols.csv gene_symbol")
        print("  python main.py my_genes.csv symbol")
        print("  python main.py data/sample/ensembl_gene_id_100.csv ensembl_gene_id")
        print("\nThis will convert your input identifiers to ALL other supported formats.")
        print("Results will be saved in the 'results/' directory.")
        sys.exit(1)
    
    input_file = sys.argv[1]
    input_column = sys.argv[2]
    
    # Check if input file exists
    if not Path(input_file).exists():
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)
    
    # Run conversion
    convert_to_all_formats(input_file, input_column)

if __name__ == '__main__':
    main()
