#!/usr/bin/env python3
"""
Comprehensive Sample File Conversion with Performance Analysis

This script converts all sample ID files to all other formats and provides
detailed performance analysis and comparison.

Usage:
    python comprehensive_sample_conversion.py
"""

import sys
import logging
import pandas as pd
from pathlib import Path
from typing import List, Dict, Optional, Tuple
import json
from datetime import datetime

# Optional visualization dependencies
try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    VISUALIZATION_AVAILABLE = True
except ImportError:
    VISUALIZATION_AVAILABLE = False
    plt = sns = np = None

# Add project root to path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from genemapkit.core.converter import GeneConverter
from genemapkit.utils.logger import setup_logging


class ComprehensiveSampleConverter:
    """Converts all sample files and analyzes conversion performance."""
    
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
        self.results_dir = Path("Results")
        self.results_dir.mkdir(parents=True, exist_ok=True)
        
        # Performance tracking
        self.performance_data = {}
    
    def convert_all_sample_files(self) -> Dict[str, Dict]:
        """
        Convert all sample files from each ID type to all other types.
        
        Returns:
            Dictionary with detailed conversion results and performance metrics
        """
        self.logger.info("🧬 Starting comprehensive sample file conversion...")
        self.logger.info(f"Results will be saved to: {self.results_dir}")
        
        conversion_results = {}
        performance_summary = {}
        
        for input_type in self.all_id_types:
            self.logger.info(f"\n{'='*70}")
            self.logger.info(f"🔄 Converting from {input_type.upper()}")
            self.logger.info(f"{'='*70}")
            
            # Check if sample file exists
            sample_file = self.sample_dir / f"{input_type}_sample.csv"
            if not sample_file.exists():
                self.logger.warning(f"❌ Sample file not found: {sample_file}")
                continue
            
            # Process this input type
            result = self._process_input_type(input_type, sample_file)
            if result:
                conversion_results[input_type] = result
                performance_summary[input_type] = result['summary']
        
        # Create overall performance analysis
        self._create_performance_analysis(performance_summary)
        
        # Create conversion matrix
        self._create_conversion_matrix(performance_summary)
        
        # Create visualization diagrams
        self._create_visualization_diagrams(performance_summary)
        
        # Generate detailed report
        self._generate_detailed_report(conversion_results)
        
        return conversion_results
    
    def _process_input_type(self, input_type: str, sample_file: Path) -> Optional[Dict]:
        """Process conversion for a specific input type."""
        try:
            # Read sample data
            df = pd.read_csv(sample_file)
            gene_ids = df[input_type].dropna().unique().tolist()
            total_genes = len(gene_ids)
            
            self.logger.info(f"📋 Loaded {total_genes} {input_type} IDs from {sample_file.name}")
            
            if total_genes == 0:
                self.logger.warning(f"No valid IDs found in {sample_file}")
                return None
            
            # Get output types (all except current input type)
            output_types = [t for t in self.all_id_types if t != input_type]
            
            # Perform conversion
            start_time = datetime.now()
            results_df = self.converter.convert_batch(
                gene_ids, input_type, output_types, include_unmapped=True
            )
            conversion_time = (datetime.now() - start_time).total_seconds()
            
            # Calculate detailed statistics
            stats = self._calculate_detailed_stats(results_df, input_type, output_types, total_genes)
            
            # Add timing information
            stats['conversion_time_seconds'] = conversion_time
            stats['genes_per_second'] = total_genes / conversion_time if conversion_time > 0 else 0
            
            # Save results to CSV
            output_file = self.results_dir / f"{input_type}ToAll.csv"
            results_df.to_csv(output_file, index=False)
            
            # Add mapping success and statistics columns
            results_df['mapping_success'] = results_df[output_types].notna().any(axis=1)
            results_df['mapped_types_count'] = results_df[output_types].notna().sum(axis=1)
            results_df['mapping_completeness_percent'] = (results_df['mapped_types_count'] / len(output_types) * 100).round(1)
            
            # Save enhanced results
            enhanced_output_file = self.results_dir / f"{input_type}ToAll_detailed.csv"
            results_df.to_csv(enhanced_output_file, index=False)
            
            self.logger.info(f"💾 Results saved to: {output_file}")
            self.logger.info(f"💾 Detailed results saved to: {enhanced_output_file}")
            
            # Log performance summary
            overall_success_rate = stats['overall_mapping_rate']
            best_target = max(stats['per_target_stats'].items(), key=lambda x: x[1]['success_rate'])
            worst_target = min(stats['per_target_stats'].items(), key=lambda x: x[1]['success_rate'])
            
            self.logger.info(f"📊 Performance Summary:")
            self.logger.info(f"   Overall success rate: {overall_success_rate:.1f}%")
            self.logger.info(f"   Best target: {best_target[0]} ({best_target[1]['success_rate']:.1f}%)")
            self.logger.info(f"   Worst target: {worst_target[0]} ({worst_target[1]['success_rate']:.1f}%)")
            self.logger.info(f"   Conversion time: {conversion_time:.2f}s ({stats['genes_per_second']:.1f} genes/sec)")
            
            return {
                'input_type': input_type,
                'sample_file': str(sample_file),
                'output_file': str(output_file),
                'results_dataframe': results_df,
                'summary': stats
            }
            
        except Exception as e:
            self.logger.error(f"❌ Error processing {input_type}: {str(e)}")
            import traceback
            traceback.print_exc()
            return None
    
    def _calculate_detailed_stats(self, results_df: pd.DataFrame, input_type: str, 
                                 output_types: List[str], total_genes: int) -> Dict:
        """Calculate detailed conversion statistics."""
        
        # Overall statistics
        mapped_any = (results_df[output_types].notna().any(axis=1)).sum()
        overall_mapping_rate = (mapped_any / total_genes) * 100
        
        # Per-target statistics
        per_target_stats = {}
        for target_type in output_types:
            mapped_count = results_df[target_type].notna().sum()
            success_rate = (mapped_count / total_genes) * 100
            
            per_target_stats[target_type] = {
                'mapped_count': mapped_count,
                'total_count': total_genes,
                'success_rate': success_rate,
                'unmapped_count': total_genes - mapped_count
            }
        
        # Completeness analysis
        completeness_scores = results_df[output_types].notna().sum(axis=1)
        fully_mapped = (completeness_scores == len(output_types)).sum()
        partially_mapped = ((completeness_scores > 0) & (completeness_scores < len(output_types))).sum()
        unmapped = (completeness_scores == 0).sum()
        
        # Quality metrics
        avg_completeness = (completeness_scores / len(output_types) * 100).mean()
        
        return {
            'input_type': input_type,
            'total_genes': total_genes,
            'overall_mapped': mapped_any,
            'overall_mapping_rate': overall_mapping_rate,
            'fully_mapped': fully_mapped,
            'partially_mapped': partially_mapped,
            'unmapped': unmapped,
            'average_completeness_percent': avg_completeness,
            'per_target_stats': per_target_stats
        }
    
    def _create_performance_analysis(self, performance_summary: Dict[str, Dict]) -> None:
        """Create comprehensive performance analysis."""
        self.logger.info("📊 Creating performance analysis...")
        
        # Create performance comparison DataFrame
        comparison_data = []
        
        for input_type, stats in performance_summary.items():
            row = {
                'input_type': input_type,
                'total_genes': stats['total_genes'],
                'overall_success_rate': stats['overall_mapping_rate'],
                'fully_mapped_count': stats['fully_mapped'],
                'partially_mapped_count': stats['partially_mapped'],
                'unmapped_count': stats['unmapped'],
                'avg_completeness_percent': stats['average_completeness_percent'],
                'conversion_time_seconds': stats.get('conversion_time_seconds', 0),
                'genes_per_second': stats.get('genes_per_second', 0)
            }
            
            # Add per-target success rates
            for target_type in self.all_id_types:
                if target_type != input_type and target_type in [t for t in self.all_id_types if t != input_type]:
                    target_stats = stats['per_target_stats'].get(target_type, {})
                    row[f'to_{target_type}_success_rate'] = target_stats.get('success_rate', 0)
            
            comparison_data.append(row)
        
        comparison_df = pd.DataFrame(comparison_data)
        
        # Save performance comparison
        performance_file = self.results_dir / "PerformanceComparison.csv"
        comparison_df.to_csv(performance_file, index=False)
        
        self.logger.info(f"💾 Performance comparison saved to: {performance_file}")
        
        # Print top performers
        print("\n🏆 TOP PERFORMERS:")
        print("=" * 60)
        
        # Best overall success rates
        best_overall = comparison_df.nlargest(3, 'overall_success_rate')
        print("\n📈 Best Overall Success Rates:")
        for _, row in best_overall.iterrows():
            print(f"  {row['input_type']:15}: {row['overall_success_rate']:.1f}% overall success")
        
        # Best completeness
        best_completeness = comparison_df.nlargest(3, 'avg_completeness_percent')
        print("\n🎯 Best Average Completeness:")
        for _, row in best_completeness.iterrows():
            print(f"  {row['input_type']:15}: {row['avg_completeness_percent']:.1f}% average completeness")
        
        # Fastest conversions
        fastest = comparison_df.nlargest(3, 'genes_per_second')
        print("\n⚡ Fastest Conversions:")
        for _, row in fastest.iterrows():
            print(f"  {row['input_type']:15}: {row['genes_per_second']:.1f} genes/second")
    
    def _create_conversion_matrix(self, performance_summary: Dict[str, Dict]) -> None:
        """Create a conversion success matrix."""
        self.logger.info("🎯 Creating conversion success matrix...")
        
        # Create matrix data
        matrix_data = []
        
        for input_type in self.all_id_types:
            if input_type not in performance_summary:
                continue
                
            row = {'input_type': input_type}
            stats = performance_summary[input_type]
            
            for target_type in self.all_id_types:
                if target_type == input_type:
                    row[target_type] = 100.0  # Self-mapping is always 100%
                elif target_type in stats['per_target_stats']:
                    row[target_type] = stats['per_target_stats'][target_type]['success_rate']
                else:
                    row[target_type] = 0.0
                    
            matrix_data.append(row)
        
        matrix_df = pd.DataFrame(matrix_data)
        
        # Save matrix
        matrix_file = self.results_dir / "ConversionSuccessMatrix.csv"
        matrix_df.to_csv(matrix_file, index=False)
        
        self.logger.info(f"💾 Conversion matrix saved to: {matrix_file}")
        
        # Print matrix
        print("\n🎯 Conversion Success Matrix (% successfully mapped):")
        print("=" * 100)
        print(matrix_df.to_string(index=False, float_format='%.1f'))
    
    def _create_visualization_diagrams(self, performance_summary: Dict[str, Dict]) -> None:
        """Create visualization diagrams for conversion success rates."""
        self.logger.info("📊 Creating visualization diagrams...")
        
        if not VISUALIZATION_AVAILABLE:
            self.logger.warning("⚠️  Matplotlib/Seaborn not available. Skipping visualization diagrams.")
            self.logger.info("💡 Install with: pip install matplotlib seaborn")
            # Create text-based diagram as fallback
            self._create_text_diagram(performance_summary)
            return
        
        try:
            # Set up the plotting style
            plt.style.use('default')
            sns.set_palette("husl")
            
            # Create figure with multiple subplots
            fig = plt.figure(figsize=(20, 15))
            
            # 1. Conversion Success Heatmap
            self._create_heatmap(performance_summary, fig)
            
            # 2. Overall Success Rates Bar Chart
            self._create_overall_success_chart(performance_summary, fig)
            
            # 3. Target-specific Performance Chart
            self._create_target_performance_chart(performance_summary, fig)
            
            # 4. Network-style Conversion Diagram
            self._create_network_diagram(performance_summary, fig)
            
            # Save the combined figure
            diagram_file = self.results_dir / "ConversionSuccessDiagrams.png"
            plt.tight_layout()
            plt.savefig(diagram_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            self.logger.info(f"📊 Visualization diagrams saved to: {diagram_file}")
            
            # Create individual specialized diagrams
            self._create_detailed_heatmap(performance_summary)
            self._create_success_flow_diagram(performance_summary)
            
        except Exception as e:
            self.logger.error(f"❌ Error creating visualization diagrams: {str(e)}")
            import traceback
            traceback.print_exc()
    
    def _create_heatmap(self, performance_summary: Dict[str, Dict], fig) -> None:
        """Create a heatmap showing conversion success rates."""
        # Prepare data for heatmap
        matrix_data = []
        id_types = sorted(self.all_id_types)
        
        for input_type in id_types:
            row = []
            for target_type in id_types:
                if input_type == target_type:
                    row.append(100.0)  # Self-mapping
                elif input_type in performance_summary:
                    stats = performance_summary[input_type].get('per_target_stats', {})
                    if target_type in stats:
                        row.append(stats[target_type]['success_rate'])
                    else:
                        row.append(0.0)
                else:
                    row.append(0.0)
            matrix_data.append(row)
        
        # Create heatmap
        ax1 = plt.subplot(2, 2, 1)
        heatmap_data = np.array(matrix_data)
        
        im = ax1.imshow(heatmap_data, cmap='RdYlGn', aspect='auto', vmin=0, vmax=100)
        
        # Set labels
        ax1.set_xticks(range(len(id_types)))
        ax1.set_yticks(range(len(id_types)))
        ax1.set_xticklabels([t.replace('_', '\n') for t in id_types], rotation=45, ha='right')
        ax1.set_yticklabels([t.replace('_', '\n') for t in id_types])
        
        # Add text annotations
        for i in range(len(id_types)):
            for j in range(len(id_types)):
                text = ax1.text(j, i, f'{heatmap_data[i, j]:.0f}%',
                               ha="center", va="center", color="black", fontsize=8)
        
        ax1.set_title('Conversion Success Rate Matrix\n(From Row → To Column)', fontsize=12, fontweight='bold')
        ax1.set_xlabel('Target ID Type')
        ax1.set_ylabel('Source ID Type')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax1, shrink=0.8)
        cbar.set_label('Success Rate (%)', rotation=270, labelpad=15)
    
    def _create_overall_success_chart(self, performance_summary: Dict[str, Dict], fig) -> None:
        """Create bar chart of overall success rates."""
        ax2 = plt.subplot(2, 2, 2)
        
        # Prepare data
        id_types = []
        success_rates = []
        colors = []
        
        for input_type, stats in performance_summary.items():
            id_types.append(input_type.replace('_', '\n'))
            rate = stats['overall_mapping_rate']
            success_rates.append(rate)
            
            # Color coding based on performance
            if rate >= 95:
                colors.append('#2E8B57')  # Green for excellent
            elif rate >= 80:
                colors.append('#FFD700')  # Gold for good
            elif rate >= 60:
                colors.append('#FFA500')  # Orange for moderate
            else:
                colors.append('#DC143C')  # Red for poor
        
        # Create bar chart
        bars = ax2.bar(id_types, success_rates, color=colors, alpha=0.8, edgecolor='black', linewidth=1)
        
        # Add value labels on bars
        for bar, rate in zip(bars, success_rates):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height + 1,
                    f'{rate:.1f}%', ha='center', va='bottom', fontweight='bold')
        
        ax2.set_title('Overall Conversion Success Rates by Source ID Type', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Success Rate (%)')
        ax2.set_xlabel('Source ID Type')
        ax2.set_ylim(0, 110)
        ax2.grid(axis='y', alpha=0.3)
        
        # Add legend
        legend_elements = [
            plt.Rectangle((0,0),1,1, facecolor='#2E8B57', label='Excellent (95-100%)'),
            plt.Rectangle((0,0),1,1, facecolor='#FFD700', label='Good (80-94%)'),
            plt.Rectangle((0,0),1,1, facecolor='#FFA500', label='Moderate (60-79%)'),
            plt.Rectangle((0,0),1,1, facecolor='#DC143C', label='Poor (<60%)')
        ]
        ax2.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1, 1))
    
    def _create_target_performance_chart(self, performance_summary: Dict[str, Dict], fig) -> None:
        """Create chart showing average performance for each target type."""
        ax3 = plt.subplot(2, 2, 3)
        
        # Calculate average success rates for each target
        target_success = {}
        for result in performance_summary.values():
            for target_type, target_stats in result['per_target_stats'].items():
                if target_type not in target_success:
                    target_success[target_type] = []
                target_success[target_type].append(target_stats['success_rate'])
        
        # Calculate averages
        target_averages = {}
        for target, rates in target_success.items():
            target_averages[target] = sum(rates) / len(rates)
        
        # Sort by performance
        sorted_targets = sorted(target_averages.items(), key=lambda x: x[1], reverse=True)
        
        target_names = [t[0].replace('_', '\n') for t, _ in sorted_targets]
        avg_rates = [rate for _, rate in sorted_targets]
        
        # Create horizontal bar chart
        bars = ax3.barh(target_names, avg_rates, color='skyblue', alpha=0.8, edgecolor='navy', linewidth=1)
        
        # Add value labels
        for bar, rate in zip(bars, avg_rates):
            width = bar.get_width()
            ax3.text(width + 1, bar.get_y() + bar.get_height()/2,
                    f'{rate:.1f}%', ha='left', va='center', fontweight='bold')
        
        ax3.set_title('Average Success Rates by Target ID Type', fontsize=12, fontweight='bold')
        ax3.set_xlabel('Average Success Rate (%)')
        ax3.set_ylabel('Target ID Type')
        ax3.set_xlim(0, max(avg_rates) + 10)
        ax3.grid(axis='x', alpha=0.3)
    
    def _create_network_diagram(self, performance_summary: Dict[str, Dict], fig) -> None:
        """Create a network-style diagram showing conversion flows."""
        ax4 = plt.subplot(2, 2, 4)
        
        # Create circular layout for ID types
        n_types = len(self.all_id_types)
        angles = np.linspace(0, 2*np.pi, n_types, endpoint=False)
        
        # Position ID types in a circle
        positions = {}
        radius = 1.0
        
        for i, id_type in enumerate(self.all_id_types):
            x = radius * np.cos(angles[i])
            y = radius * np.sin(angles[i])
            positions[id_type] = (x, y)
        
        # Draw connections based on success rates
        for input_type, stats in performance_summary.items():
            if input_type not in positions:
                continue
                
            start_pos = positions[input_type]
            
            for target_type, target_stats in stats['per_target_stats'].items():
                if target_type not in positions:
                    continue
                    
                end_pos = positions[target_type]
                success_rate = target_stats['success_rate']
                
                # Only draw lines for successful conversions
                if success_rate > 0:
                    # Line thickness based on success rate
                    linewidth = max(0.5, success_rate / 20)
                    
                    # Color based on success rate
                    if success_rate >= 95:
                        color = '#00FF00'  # Green
                        alpha = 0.8
                    elif success_rate >= 80:
                        color = '#FFFF00'  # Yellow
                        alpha = 0.6
                    elif success_rate >= 60:
                        color = '#FFA500'  # Orange
                        alpha = 0.4
                    else:
                        color = '#FF0000'  # Red
                        alpha = 0.3
                    
                    ax4.plot([start_pos[0], end_pos[0]], [start_pos[1], end_pos[1]],
                            color=color, linewidth=linewidth, alpha=alpha)
        
        # Draw ID type nodes
        for i, id_type in enumerate(self.all_id_types):
            pos = positions[id_type]
            
            # Node color based on overall performance as source
            if id_type in performance_summary:
                overall_rate = performance_summary[id_type]['overall_mapping_rate']
                if overall_rate >= 95:
                    node_color = '#2E8B57'
                elif overall_rate >= 80:
                    node_color = '#FFD700'
                elif overall_rate >= 60:
                    node_color = '#FFA500'
                else:
                    node_color = '#DC143C'
            else:
                node_color = '#808080'
            
            circle = plt.Circle(pos, 0.15, color=node_color, alpha=0.8, zorder=10)
            ax4.add_patch(circle)
            
            # Add labels
            label = id_type.replace('_', '\n')
            ax4.text(pos[0], pos[1], label, ha='center', va='center', 
                    fontsize=8, fontweight='bold', zorder=15)
        
        ax4.set_xlim(-1.5, 1.5)
        ax4.set_ylim(-1.5, 1.5)
        ax4.set_aspect('equal')
        ax4.set_title('Conversion Success Network\n(Line thickness = Success rate)', fontsize=12, fontweight='bold')
        ax4.axis('off')
    
    def _create_detailed_heatmap(self, performance_summary: Dict[str, Dict]) -> None:
        """Create a detailed standalone heatmap."""
        if not VISUALIZATION_AVAILABLE:
            return
            
        try:
            fig, ax = plt.subplots(figsize=(12, 10))
            
            # Prepare matrix data
            id_types = sorted(self.all_id_types)
            matrix_data = []
            
            for input_type in id_types:
                row = []
                for target_type in id_types:
                    if input_type == target_type:
                        row.append(100.0)
                    elif input_type in performance_summary:
                        stats = performance_summary[input_type].get('per_target_stats', {})
                        row.append(stats.get(target_type, {}).get('success_rate', 0.0))
                    else:
                        row.append(0.0)
                matrix_data.append(row)
            
            matrix_data = np.array(matrix_data)
            
            # Create heatmap with annotations
            sns.heatmap(matrix_data, 
                       xticklabels=[t.replace('_', ' ').title() for t in id_types],
                       yticklabels=[t.replace('_', ' ').title() for t in id_types],
                       annot=True, 
                       fmt='.1f',
                       cmap='RdYlGn',
                       vmin=0, vmax=100,
                       square=True,
                       ax=ax,
                       cbar_kws={'label': 'Success Rate (%)'})
            
            plt.title('Gene ID Conversion Success Rate Matrix\n(From Row → To Column)', 
                     fontsize=16, fontweight='bold', pad=20)
            plt.xlabel('Target ID Type', fontsize=12, fontweight='bold')
            plt.ylabel('Source ID Type', fontsize=12, fontweight='bold')
            
            plt.tight_layout()
            
            heatmap_file = self.results_dir / "DetailedConversionHeatmap.png"
            plt.savefig(heatmap_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            self.logger.info(f"📊 Detailed heatmap saved to: {heatmap_file}")
            
        except Exception as e:
            self.logger.error(f"❌ Error creating detailed heatmap: {str(e)}")
    
    def _create_success_flow_diagram(self, performance_summary: Dict[str, Dict]) -> None:
        """Create a flow diagram showing successful conversion paths."""
        if not VISUALIZATION_AVAILABLE:
            return
            
        try:
            fig, ax = plt.subplots(figsize=(15, 10))
            
            # Create a sankey-style flow diagram
            # Group conversions by success rate ranges
            excellent_conversions = []  # 95-100%
            good_conversions = []      # 80-94%
            moderate_conversions = []  # 60-79%
            poor_conversions = []      # <60%
            
            for input_type, stats in performance_summary.items():
                for target_type, target_stats in stats['per_target_stats'].items():
                    success_rate = target_stats['success_rate']
                    conversion = f"{input_type} → {target_type}"
                    
                    if success_rate >= 95:
                        excellent_conversions.append((conversion, success_rate))
                    elif success_rate >= 80:
                        good_conversions.append((conversion, success_rate))
                    elif success_rate >= 60:
                        moderate_conversions.append((conversion, success_rate))
                    elif success_rate > 0:
                        poor_conversions.append((conversion, success_rate))
            
            # Create text summary
            y_pos = 0.9
            ax.text(0.05, y_pos, "🧬 Gene ID Conversion Success Summary", 
                   fontsize=18, fontweight='bold', transform=ax.transAxes)
            
            y_pos -= 0.08
            ax.text(0.05, y_pos, f"Total Conversion Combinations: {len(excellent_conversions + good_conversions + moderate_conversions + poor_conversions)}", 
                   fontsize=14, transform=ax.transAxes)
            
            # Excellent conversions
            y_pos -= 0.1
            ax.text(0.05, y_pos, "🟢 Excellent Conversions (95-100% success):", 
                   fontsize=14, fontweight='bold', color='green', transform=ax.transAxes)
            
            y_pos -= 0.04
            for conversion, rate in sorted(excellent_conversions, key=lambda x: x[1], reverse=True):
                if y_pos > 0.1:
                    ax.text(0.1, y_pos, f"• {conversion}: {rate:.1f}%", 
                           fontsize=10, transform=ax.transAxes)
                    y_pos -= 0.03
            
            # Good conversions
            if good_conversions and y_pos > 0.2:
                y_pos -= 0.05
                ax.text(0.05, y_pos, "🟡 Good Conversions (80-94% success):", 
                       fontsize=14, fontweight='bold', color='orange', transform=ax.transAxes)
                
                y_pos -= 0.04
                for conversion, rate in sorted(good_conversions, key=lambda x: x[1], reverse=True):
                    if y_pos > 0.1:
                        ax.text(0.1, y_pos, f"• {conversion}: {rate:.1f}%", 
                               fontsize=10, transform=ax.transAxes)
                        y_pos -= 0.03
            
            # Add summary statistics on the right side
            stats_x = 0.55
            y_pos = 0.85
            
            ax.text(stats_x, y_pos, "📊 Performance Statistics", 
                   fontsize=16, fontweight='bold', transform=ax.transAxes)
            
            y_pos -= 0.08
            total_conversions = len(excellent_conversions + good_conversions + moderate_conversions + poor_conversions)
            if total_conversions > 0:
                excellent_pct = len(excellent_conversions) / total_conversions * 100
                good_pct = len(good_conversions) / total_conversions * 100
                moderate_pct = len(moderate_conversions) / total_conversions * 100
                poor_pct = len(poor_conversions) / total_conversions * 100
                
                ax.text(stats_x, y_pos, f"Excellent: {len(excellent_conversions)} ({excellent_pct:.1f}%)", 
                       fontsize=12, color='green', transform=ax.transAxes)
                y_pos -= 0.05
                ax.text(stats_x, y_pos, f"Good: {len(good_conversions)} ({good_pct:.1f}%)", 
                       fontsize=12, color='orange', transform=ax.transAxes)
                y_pos -= 0.05
                ax.text(stats_x, y_pos, f"Moderate: {len(moderate_conversions)} ({moderate_pct:.1f}%)", 
                       fontsize=12, color='red', transform=ax.transAxes)
                y_pos -= 0.05
                ax.text(stats_x, y_pos, f"Poor: {len(poor_conversions)} ({poor_pct:.1f}%)", 
                       fontsize=12, color='darkred', transform=ax.transAxes)
            
            # Best and worst performers
            y_pos -= 0.1
            if performance_summary:
                best_source = max(performance_summary.items(), key=lambda x: x[1]['overall_mapping_rate'])
                worst_source = min(performance_summary.items(), key=lambda x: x[1]['overall_mapping_rate'])
                
                ax.text(stats_x, y_pos, "🏆 Best Source ID Type:", 
                       fontsize=12, fontweight='bold', transform=ax.transAxes)
                y_pos -= 0.04
                ax.text(stats_x, y_pos, f"{best_source[0]} ({best_source[1]['overall_mapping_rate']:.1f}%)", 
                       fontsize=11, transform=ax.transAxes)
                
                y_pos -= 0.08
                ax.text(stats_x, y_pos, "⚠️ Worst Source ID Type:", 
                       fontsize=12, fontweight='bold', transform=ax.transAxes)
                y_pos -= 0.04
                ax.text(stats_x, y_pos, f"{worst_source[0]} ({worst_source[1]['overall_mapping_rate']:.1f}%)", 
                       fontsize=11, transform=ax.transAxes)
            
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.axis('off')
            
            plt.tight_layout()
            
            flow_file = self.results_dir / "ConversionSuccessFlow.png"
            plt.savefig(flow_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            self.logger.info(f"📊 Success flow diagram saved to: {flow_file}")
            
        except Exception as e:
            self.logger.error(f"❌ Error creating success flow diagram: {str(e)}")
    
    def _create_text_diagram(self, performance_summary: Dict[str, Dict]) -> None:
        """Create a text-based conversion success diagram."""
        self.logger.info("📊 Creating text-based conversion success diagram...")
        
        diagram_lines = [
            "# Gene ID Conversion Success Diagram",
            f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            "",
            "## Conversion Success Matrix (Text Format)",
            "",
            "Legend: ✅ Excellent (95-100%), 🟡 Good (80-94%), 🟠 Moderate (60-79%), ❌ Poor (<60%), ⚫ Failed (0%)",
            "",
        ]
        
        # Create text-based matrix
        id_types = sorted(self.all_id_types)
        
        # Header
        header = "FROM\\TO".ljust(15)
        for target in id_types:
            header += target[:8].ljust(10)
        diagram_lines.append(header)
        diagram_lines.append("=" * len(header))
        
        # Matrix rows
        for input_type in id_types:
            row = input_type[:12].ljust(15)
            
            for target_type in id_types:
                if input_type == target_type:
                    symbol = "🔄"  # Self-conversion
                    rate_str = "100%"
                elif input_type in performance_summary:
                    stats = performance_summary[input_type].get('per_target_stats', {})
                    if target_type in stats:
                        rate = stats[target_type]['success_rate']
                        rate_str = f"{rate:.0f}%"
                        
                        if rate >= 95:
                            symbol = "✅"
                        elif rate >= 80:
                            symbol = "🟡"
                        elif rate >= 60:
                            symbol = "🟠"
                        elif rate > 0:
                            symbol = "❌"
                        else:
                            symbol = "⚫"
                    else:
                        symbol = "⚫"
                        rate_str = "0%"
                else:
                    symbol = "⚫"
                    rate_str = "0%"
                
                cell = f"{symbol}{rate_str}".ljust(10)
                row += cell
            
            diagram_lines.append(row)
        
        diagram_lines.extend([
            "",
            "## Success Rate Summary by Source ID Type",
            ""
        ])
        
        # Source performance summary
        for input_type, stats in sorted(performance_summary.items(), 
                                      key=lambda x: x[1]['overall_mapping_rate'], 
                                      reverse=True):
            rate = stats['overall_mapping_rate']
            if rate >= 95:
                emoji = "🏆"
            elif rate >= 80:
                emoji = "🥈"
            elif rate >= 60:
                emoji = "🥉"
            else:
                emoji = "⚠️"
            
            diagram_lines.append(f"{emoji} {input_type.upper()}: {rate:.1f}% overall success")
        
        diagram_lines.extend([
            "",
            "## Best Conversion Paths",
            ""
        ])
        
        # Find best conversions
        best_conversions = []
        for input_type, stats in performance_summary.items():
            for target_type, target_stats in stats['per_target_stats'].items():
                if target_stats['success_rate'] >= 95:
                    best_conversions.append((
                        f"{input_type} → {target_type}",
                        target_stats['success_rate']
                    ))
        
        best_conversions.sort(key=lambda x: x[1], reverse=True)
        
        for conversion, rate in best_conversions:
            diagram_lines.append(f"✅ {conversion}: {rate:.1f}%")
        
        diagram_lines.extend([
            "",
            "## Problem Areas (Low Success Rates)",
            ""
        ])
        
        # Find problematic conversions
        problem_conversions = []
        for input_type, stats in performance_summary.items():
            for target_type, target_stats in stats['per_target_stats'].items():
                if target_stats['success_rate'] < 60:
                    problem_conversions.append((
                        f"{input_type} → {target_type}",
                        target_stats['success_rate']
                    ))
        
        problem_conversions.sort(key=lambda x: x[1])
        
        for conversion, rate in problem_conversions:
            if rate == 0:
                diagram_lines.append(f"⚫ {conversion}: {rate:.1f}% (No mappings found)")
            else:
                diagram_lines.append(f"❌ {conversion}: {rate:.1f}%")
        
        # Save text diagram
        text_diagram_file = self.results_dir / "ConversionSuccessDiagram.txt"
        with open(text_diagram_file, 'w') as f:
            f.write('\n'.join(diagram_lines))
        
        self.logger.info(f"📊 Text-based conversion diagram saved to: {text_diagram_file}")
        
        # Print summary to console
        print("\n" + "="*80)
        print("🧬 GENE ID CONVERSION SUCCESS DIAGRAM")
        print("="*80)
        
        print("\n🏆 TOP PERFORMING SOURCE ID TYPES:")
        for input_type, stats in sorted(performance_summary.items(), 
                                      key=lambda x: x[1]['overall_mapping_rate'], 
                                      reverse=True)[:3]:
            rate = stats['overall_mapping_rate']
            print(f"   {input_type.upper()}: {rate:.1f}% overall success")
        
        print(f"\n📊 Best Conversion Paths ({len(best_conversions)} excellent conversions):")
        for conversion, rate in best_conversions[:10]:  # Show top 10
            print(f"   ✅ {conversion}: {rate:.1f}%")
        
        if problem_conversions:
            print(f"\n⚠️  Problem Areas ({len(problem_conversions)} low-success conversions):")
            for conversion, rate in problem_conversions[:5]:  # Show worst 5
                if rate == 0:
                    print(f"   ⚫ {conversion}: No mappings found")
                else:
                    print(f"   ❌ {conversion}: {rate:.1f}%")
        
        print(f"\n💾 Full text diagram saved to: {text_diagram_file}")
        print("="*80)
    
    def _generate_detailed_report(self, conversion_results: Dict[str, Dict]) -> None:
        """Generate a detailed markdown report."""
        self.logger.info("📋 Generating detailed report...")
        
        report_lines = [
            "# Comprehensive Gene ID Conversion Analysis Report",
            f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"**Total ID Types Analyzed:** {len(conversion_results)}",
            "",
            "## Executive Summary",
            "",
        ]
        
        # Calculate overall statistics
        total_conversions = sum(len(result['summary']['per_target_stats']) for result in conversion_results.values())
        avg_success_rate = sum(result['summary']['overall_mapping_rate'] for result in conversion_results.values()) / len(conversion_results)
        
        report_lines.extend([
            f"- **Total Conversion Combinations Tested:** {total_conversions}",
            f"- **Average Success Rate:** {avg_success_rate:.1f}%",
            f"- **Sample Size per ID Type:** 20 genes",
            "",
            "## Detailed Results by Input Type",
            ""
        ])
        
        # Add detailed results for each input type
        for input_type, result in conversion_results.items():
            stats = result['summary']
            report_lines.extend([
                f"### {input_type.upper()}",
                f"- **Overall Success Rate:** {stats['overall_mapping_rate']:.1f}%",
                f"- **Fully Mapped Genes:** {stats['fully_mapped']}/{stats['total_genes']} ({stats['fully_mapped']/stats['total_genes']*100:.1f}%)",
                f"- **Partially Mapped Genes:** {stats['partially_mapped']}/{stats['total_genes']} ({stats['partially_mapped']/stats['total_genes']*100:.1f}%)",
                f"- **Unmapped Genes:** {stats['unmapped']}/{stats['total_genes']} ({stats['unmapped']/stats['total_genes']*100:.1f}%)",
                f"- **Average Completeness:** {stats['average_completeness_percent']:.1f}%",
                "",
                "**Target-Specific Success Rates:**",
            ])
            
            # Sort targets by success rate
            sorted_targets = sorted(stats['per_target_stats'].items(), 
                                  key=lambda x: x[1]['success_rate'], reverse=True)
            
            for target_type, target_stats in sorted_targets:
                report_lines.append(f"- {target_type}: {target_stats['success_rate']:.1f}% ({target_stats['mapped_count']}/{target_stats['total_count']})")
            
            report_lines.append("")
        
        # Add recommendations
        report_lines.extend([
            "## Recommendations",
            "",
            "### Best Source ID Types (Highest Overall Success):",
        ])
        
        # Rank input types by overall success
        ranked_inputs = sorted(conversion_results.items(), 
                             key=lambda x: x[1]['summary']['overall_mapping_rate'], 
                             reverse=True)
        
        for i, (input_type, result) in enumerate(ranked_inputs[:3], 1):
            rate = result['summary']['overall_mapping_rate']
            report_lines.append(f"{i}. **{input_type}** - {rate:.1f}% overall success rate")
        
        report_lines.extend([
            "",
            "### Best Target Mappings (Across All Sources):",
        ])
        
        # Find best target mappings
        target_success = {}
        for result in conversion_results.values():
            for target_type, target_stats in result['summary']['per_target_stats'].items():
                if target_type not in target_success:
                    target_success[target_type] = []
                target_success[target_type].append(target_stats['success_rate'])
        
        # Calculate average success rates for each target
        target_averages = {target: sum(rates)/len(rates) for target, rates in target_success.items()}
        best_targets = sorted(target_averages.items(), key=lambda x: x[1], reverse=True)
        
        for i, (target_type, avg_rate) in enumerate(best_targets[:3], 1):
            report_lines.append(f"{i}. **{target_type}** - {avg_rate:.1f}% average success rate")
        
        # Save report
        report_file = self.results_dir / "ComprehensiveAnalysisReport.md"
        with open(report_file, 'w') as f:
            f.write('\n'.join(report_lines))
        
        self.logger.info(f"💾 Detailed report saved to: {report_file}")


def main():
    """Main function."""
    print("🧬 Comprehensive Sample File Conversion & Performance Analysis")
    print("=" * 80)
    
    # Setup logging
    logger = setup_logging('INFO')
    
    # Initialize converter
    converter = ComprehensiveSampleConverter(logger=logger)
    
    try:
        # Check if sample files exist
        missing_files = []
        available_files = []
        
        for id_type in converter.all_id_types:
            sample_file = converter.sample_dir / f"{id_type}_sample.csv"
            if sample_file.exists():
                available_files.append(id_type)
            else:
                missing_files.append(str(sample_file))
        
        logger.info(f"📋 Found {len(available_files)} sample files: {available_files}")
        
        if missing_files:
            logger.warning(f"⚠️  Missing sample files: {missing_files}")
        
        if not available_files:
            logger.error("❌ No sample files found! Please create sample files first.")
            return
        
        # Perform comprehensive conversions
        results = converter.convert_all_sample_files()
        
        logger.info("\n🎉 Comprehensive conversion analysis completed!")
        logger.info(f"📁 Check the '{converter.results_dir}' directory for all results:")
        logger.info(f"   - Individual conversion files: *ToAll.csv")
        logger.info(f"   - Enhanced results with statistics: *ToAll_detailed.csv")
        logger.info(f"   - Performance comparison: PerformanceComparison.csv")
        logger.info(f"   - Conversion matrix: ConversionSuccessMatrix.csv")
        logger.info(f"   - Detailed report: ComprehensiveAnalysisReport.md")
        logger.info(f"   - Text-based diagram: ConversionSuccessDiagram.txt")
        
        if VISUALIZATION_AVAILABLE:
            logger.info(f"   - Visualization diagrams: ConversionSuccessDiagrams.png")
            logger.info(f"   - Detailed heatmap: DetailedConversionHeatmap.png")
            logger.info(f"   - Success flow diagram: ConversionSuccessFlow.png")
            logger.info("\n📊 Visualization Features:")
            logger.info(f"   - Conversion success rate heatmap")
            logger.info(f"   - Overall success rates bar chart")
            logger.info(f"   - Target performance analysis")
            logger.info(f"   - Network-style conversion flow diagram")
        else:
            logger.info(f"\n💡 For visual diagrams, install: pip install matplotlib seaborn")
        
    except Exception as e:
        logger.error(f"❌ Error during comprehensive conversion: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
