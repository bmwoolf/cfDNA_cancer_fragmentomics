#!/usr/bin/env python3
"""
Batch cfDNA Fragmentomics Analysis Script

This script processes multiple BAM files and generates comparative
fragment length analysis reports for cancer research applications.
"""

import argparse
import sys
import os
from pathlib import Path
from typing import Dict, List, Tuple
import logging
import glob

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

# Import our main analyzer
from parse_cfDNA_frag_lengths import cfDNAFragmentAnalyzer

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class BatchFragmentAnalyzer:
    """Batch analysis of multiple BAM files for fragmentomics."""
    
    def __init__(self, input_pattern: str, output_dir: str = "outputs"):
        """
        Initialize batch analyzer.
        
        Args:
            input_pattern: Glob pattern for BAM files (e.g., "data/*.bam")
            output_dir: Output directory
        """
        self.input_pattern = input_pattern
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Find BAM files
        self.bam_files = glob.glob(input_pattern)
        if not self.bam_files:
            raise FileNotFoundError(f"No BAM files found matching pattern: {input_pattern}")
        
        logger.info(f"Found {len(self.bam_files)} BAM files to process")
        
        # Results storage
        self.results = {}
        self.combined_stats = []
    
    def process_all_files(self, max_reads_per_file: int = None) -> Dict:
        """
        Process all BAM files in the batch.
        
        Args:
            max_reads_per_file: Maximum reads to process per file
            
        Returns:
            Dict: Combined results
        """
        logger.info(f"Starting batch analysis of {len(self.bam_files)} files")
        
        for bam_file in tqdm(self.bam_files, desc="Processing BAM files"):
            try:
                # Extract sample name from filename
                sample_name = Path(bam_file).stem
                
                # Create sample-specific output directory
                sample_output_dir = self.output_dir / sample_name
                sample_output_dir.mkdir(exist_ok=True)
                
                # Analyze this file
                analyzer = cfDNAFragmentAnalyzer(bam_file, str(sample_output_dir))
                stats = analyzer.analyze_fragments(max_reads=max_reads_per_file)
                
                if stats:
                    # Store results
                    self.results[sample_name] = {
                        'analyzer': analyzer,
                        'stats': stats,
                        'bam_file': bam_file
                    }
                    
                    # Save individual results
                    analyzer.save_fragment_lengths_csv(f"{sample_name}_fragment_lengths.csv")
                    analyzer.create_histogram(f"{sample_name}_histogram.png")
                    
                    # Add to combined stats
                    stats_copy = stats.copy()
                    stats_copy['sample_name'] = sample_name
                    stats_copy['bam_file'] = bam_file
                    self.combined_stats.append(stats_copy)
                    
                    logger.info(f"Successfully processed {sample_name}")
                else:
                    logger.warning(f"No valid fragments found in {sample_name}")
                    
            except Exception as e:
                logger.error(f"Error processing {bam_file}: {e}")
                continue
        
        # Save combined results
        if self.combined_stats:
            self.save_combined_results()
            self.create_comparative_plots()
        
        return self.results
    
    def save_combined_results(self):
        """Save combined statistics from all samples."""
        if not self.combined_stats:
            logger.warning("No results to save!")
            return
        
        # Create combined DataFrame
        combined_df = pd.DataFrame(self.combined_stats)
        
        # Reorder columns for better readability
        column_order = [
            'sample_name', 'bam_file', 'total_reads_processed', 'filtered_read_pairs',
            'mean_fragment_length', 'median_fragment_length', 'std_fragment_length',
            'min_fragment_length', 'max_fragment_length',
            'fragment_length_25th_percentile', 'fragment_length_75th_percentile'
        ]
        
        # Only include columns that exist
        existing_columns = [col for col in column_order if col in combined_df.columns]
        combined_df = combined_df[existing_columns]
        
        # Save combined statistics
        output_path = self.output_dir / "combined_fragment_statistics.csv"
        combined_df.to_csv(output_path, index=False)
        logger.info(f"Combined statistics saved to: {output_path}")
        
        # Create summary report
        self.create_summary_report(combined_df)
    
    def create_summary_report(self, combined_df: pd.DataFrame):
        """Create a summary report of all samples."""
        report_path = self.output_dir / "analysis_summary_report.txt"
        
        with open(report_path, 'w') as f:
            f.write("cfDNA Fragmentomics Batch Analysis Summary Report\n")
            f.write("=" * 60 + "\n\n")
            
            f.write(f"Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Total Samples Processed: {len(combined_df)}\n")
            f.write(f"Total Read Pairs Analyzed: {combined_df['filtered_read_pairs'].sum():,}\n\n")
            
            f.write("Sample Summary:\n")
            f.write("-" * 40 + "\n")
            for _, row in combined_df.iterrows():
                f.write(f"Sample: {row['sample_name']}\n")
                f.write(f"  Total Reads: {row['total_reads_processed']:,}\n")
                f.write(f"  Valid Fragments: {row['filtered_read_pairs']:,}\n")
                f.write(f"  Mean Length: {row['mean_fragment_length']:.1f} bp\n")
                f.write(f"  Median Length: {row['median_fragment_length']:.1f} bp\n")
                f.write(f"  Std Dev: {row['std_fragment_length']:.1f} bp\n\n")
            
            f.write("Overall Statistics:\n")
            f.write("-" * 40 + "\n")
            f.write(f"Mean Fragment Length (across samples): {combined_df['mean_fragment_length'].mean():.1f} bp\n")
            f.write(f"Median Fragment Length (across samples): {combined_df['median_fragment_length'].median():.1f} bp\n")
            f.write(f"Fragment Length Range: {combined_df['min_fragment_length'].min()} - {combined_df['max_fragment_length'].max()} bp\n")
            f.write(f"Average Filtering Rate: {combined_df['filtered_read_pairs'].sum() / combined_df['total_reads_processed'].sum() * 100:.2f}%\n")
        
        logger.info(f"Summary report saved to: {report_path}")
    
    def create_comparative_plots(self):
        """Create comparative plots across all samples."""
        if not self.combined_stats:
            logger.warning("No data for comparative plots!")
            return
        
        # Create comparative histogram
        plt.figure(figsize=(15, 10))
        
        # Plot fragment length distributions for each sample
        for sample_name, result in self.results.items():
            fragment_lengths = result['analyzer'].fragment_lengths
            plt.hist(fragment_lengths, bins=30, alpha=0.6, label=sample_name, density=True)
        
        plt.xlabel('Fragment Length (bp)')
        plt.ylabel('Density')
        plt.title('Comparative cfDNA Fragment Length Distributions')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Save comparative plot
        output_path = self.output_dir / "comparative_fragment_distributions.png"
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Comparative plot saved to: {output_path}")
        
        # Create box plot of fragment length statistics
        self.create_box_plot()
    
    def create_box_plot(self):
        """Create box plot comparing fragment length statistics across samples."""
        if not self.combined_stats:
            return
        
        # Prepare data for box plot
        box_data = []
        labels = []
        
        for sample_name, result in self.results.items():
            fragment_lengths = result['analyzer'].fragment_lengths
            box_data.append(fragment_lengths)
            labels.append(sample_name)
        
        # Create box plot
        plt.figure(figsize=(12, 8))
        plt.boxplot(box_data, labels=labels)
        plt.xlabel('Sample')
        plt.ylabel('Fragment Length (bp)')
        plt.title('Fragment Length Distribution Comparison')
        plt.xticks(rotation=45)
        plt.grid(True, alpha=0.3)
        
        # Save box plot
        output_path = self.output_dir / "fragment_length_boxplot.png"
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Box plot saved to: {output_path}")
    
    def print_batch_summary(self):
        """Print batch analysis summary to console."""
        if not self.combined_stats:
            logger.warning("No results to summarize!")
            return
        
        print("\n" + "="*80)
        print("BATCH cfDNA Fragmentomics Analysis Summary")
        print("="*80)
        print(f"Total Samples Processed: {len(self.combined_stats)}")
        print(f"Total Fragment Pairs Analyzed: {sum(s['filtered_read_pairs'] for s in self.combined_stats):,}")
        
        # Calculate overall statistics
        mean_lengths = [s['mean_fragment_length'] for s in self.combined_stats]
        median_lengths = [s['median_fragment_length'] for s in self.combined_stats]
        
        print(f"\nOverall Fragment Length Statistics:")
        print(f"  Mean Length (across samples): {np.mean(mean_lengths):.1f} ± {np.std(mean_lengths):.1f} bp")
        print(f"  Median Length (across samples): {np.mean(median_lengths):.1f} ± {np.std(median_lengths):.1f} bp")
        
        print(f"\nSample Details:")
        for stats in self.combined_stats:
            print(f"  {stats['sample_name']}: {stats['filtered_read_pairs']:,} fragments, "
                  f"mean={stats['mean_fragment_length']:.1f}bp")
        
        print("="*80)


def main():
    """Main function for batch fragment analysis."""
    parser = argparse.ArgumentParser(
        description="Batch analysis of cfDNA fragment lengths from multiple BAM files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
            Examples:
            python batch_fragment_analysis.py "data/*.bam"
            python batch_fragment_analysis.py "samples/*.bam" --max-reads 5000 --output-dir results
            python batch_fragment_analysis.py "*.bam" --verbose
        """
    )
    
    parser.add_argument(
        'input_pattern',
        help='Glob pattern for BAM files (e.g., "data/*.bam")'
    )
    
    parser.add_argument(
        '--output-dir',
        default='outputs',
        help='Output directory (default: outputs)'
    )
    
    parser.add_argument(
        '--max-reads',
        type=int,
        help='Maximum read pairs to process per file (for testing)'
    )
    
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )
    
    args = parser.parse_args()
    
    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    try:
        # Initialize batch analyzer
        batch_analyzer = BatchFragmentAnalyzer(args.input_pattern, args.output_dir)
        
        # Process all files
        results = batch_analyzer.process_all_files(max_reads_per_file=args.max_reads)
        
        if results:
            # Print summary
            batch_analyzer.print_batch_summary()
            logger.info("Batch analysis completed successfully!")
        else:
            logger.error("Batch analysis failed - no valid results!")
            sys.exit(1)
            
    except Exception as e:
        logger.error(f"Error during batch analysis: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main() 