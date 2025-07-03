#!/usr/bin/env python3
"""
cfDNA Fragmentomics Analysis Script

This script extracts cfDNA fragment length distributions from BAM files
for cancer research applications, specifically for optimizing amplicon sizes
in qPCR assays.
"""

import argparse
import sys
import os
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import logging

import pysam
import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class cfDNAFragmentAnalyzer:
    """Analyzes cfDNA fragment lengths from BAM files."""
    
    def __init__(self, bam_path: str, output_dir: str = "outputs"):
        """
        Initialize the analyzer.
        
        Args:
            bam_path: Path to the BAM file
            output_dir: Directory to save outputs
        """
        self.bam_path = bam_path
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Fragment length statistics
        self.fragment_lengths = []
        self.stats = {}
        
        # Validate BAM file
        if not os.path.exists(bam_path):
            raise FileNotFoundError(f"BAM file not found: {bam_path}")
        
        # Check if BAM is indexed
        bai_path = bam_path + ".bai"
        if not os.path.exists(bai_path):
            logger.warning(f"BAM index not found: {bai_path}")
            logger.info("Creating BAM index...")
            pysam.index(bam_path)
    
    def filter_read_pair(self, read1, read2) -> bool:
        """
        Filter read pairs based on quality criteria.
        
        Args:
            read1: First read in pair
            read2: Second read in pair
            
        Returns:
            bool: True if read pair passes filters
        """
        # Check if reads are properly paired
        if not (read1.is_proper_pair and read2.is_proper_pair):
            return False
        
        # Check mapping quality (MAPQ >= 30)
        if read1.mapping_quality < 30 or read2.mapping_quality < 30:
            return False
        
        # Ignore supplementary or secondary alignments
        if read1.is_supplementary or read1.is_secondary:
            return False
        if read2.is_supplementary or read2.is_secondary:
            return False
        
        # Check if reads are on the same chromosome
        if read1.reference_name != read2.reference_name:
            return False
        
        # For testing, limit to chr1
        if read1.reference_name != "chr1":
            return False
        
        # Check if reads are not too far apart (reasonable fragment size)
        insert_size = abs(read1.template_length)
        if insert_size > 1000 or insert_size < 50:  # cfDNA typically 50-1000bp
            return False
        
        return True
    
    def calculate_fragment_length(self, read1, read2) -> int:
        """
        Calculate fragment length (insert size) for a read pair.
        
        Args:
            read1: First read in pair
            read2: Second read in pair
            
        Returns:
            int: Fragment length
        """
        return abs(read1.template_length)
    
    def analyze_fragments(self, max_reads: Optional[int] = None) -> Dict:
        """
        Analyze fragment lengths from BAM file.
        
        Args:
            max_reads: Maximum number of read pairs to process (for testing)
            
        Returns:
            Dict: Analysis statistics
        """
        logger.info(f"Starting fragment analysis of {self.bam_path}")
        
        with pysam.AlignmentFile(self.bam_path, "rb") as bam:
            total_reads = 0
            filtered_reads = 0
            
            # Process only first read in each pair
            for read in tqdm(bam.fetch(until_eof=True), desc="Processing reads"):
                total_reads += 1
                
                # Only process first read in pair
                if not read.is_read1:
                    continue
                
                # Apply filters (mimic filter_read_pair logic, but for single read)
                if not read.is_proper_pair:
                    continue
                if read.mapping_quality < 30:
                    continue
                if read.is_supplementary or read.is_secondary:
                    continue
                if read.reference_name != "chr1":
                    continue
                insert_size = abs(read.template_length)
                if insert_size > 1000 or insert_size < 50:
                    continue
                
                # Calculate fragment length
                fragment_length = insert_size
                self.fragment_lengths.append(fragment_length)
                filtered_reads += 1
                
                # Limit reads for testing
                if max_reads and filtered_reads >= max_reads:
                    break
        
        # Calculate statistics
        if self.fragment_lengths:
            self.stats = self.calculate_stats_from_list(self.fragment_lengths)
            self.stats['filtered_read_pairs'] = len(self.fragment_lengths)
            self.stats['total_reads_processed'] = total_reads
        else:
            logger.warning("No valid fragment lengths found!")
            self.stats = {}
        
        logger.info(f"Analysis complete. Found {len(self.fragment_lengths)} valid fragment lengths.")
        return self.stats
    
    def save_fragment_lengths_csv(self, filename: str = "fragment_lengths.csv") -> str:
        """
        Save fragment lengths to CSV file.
        
        Args:
            filename: Output filename
            
        Returns:
            str: Path to saved file
        """
        if not self.fragment_lengths:
            logger.warning("No fragment lengths to save!")
            return ""
        
        # Create DataFrame
        df = pd.DataFrame({
            'fragment_length': self.fragment_lengths
        })
        
        # Add statistics
        stats_df = pd.DataFrame([self.stats])
        
        # Save fragment lengths
        output_path = self.output_dir / filename
        df.to_csv(output_path, index=False)
        logger.info(f"Fragment lengths saved to: {output_path}")
        
        # Save statistics
        stats_path = self.output_dir / "fragment_statistics.csv"
        stats_df.to_csv(stats_path, index=False)
        logger.info(f"Statistics saved to: {stats_path}")
        
        return str(output_path)
    
    def create_histogram(self, filename: str = "fragment_length_histogram.png") -> str:
        """
        Create and save fragment length histogram.
        
        Args:
            filename: Output filename
            
        Returns:
            str: Path to saved file
        """
        if not self.fragment_lengths:
            logger.warning("No fragment lengths to plot!")
            return ""
        
        # Set up the plot
        plt.figure(figsize=(12, 8))
        
        # Create histogram
        plt.hist(self.fragment_lengths, bins=50, alpha=0.7, color='steelblue', edgecolor='black')
        
        # Add vertical lines for mean and median
        plt.axvline(self.stats['mean_fragment_length'], color='red', linestyle='--', 
                   label=f'Mean: {self.stats["mean_fragment_length"]:.1f}bp')
        plt.axvline(self.stats['median_fragment_length'], color='orange', linestyle='--', 
                   label=f'Median: {self.stats["median_fragment_length"]:.1f}bp')
        
        # Customize plot
        plt.xlabel('Fragment Length (bp)')
        plt.ylabel('Frequency')
        plt.title('cfDNA Fragment Length Distribution')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Save plot
        output_path = self.output_dir / filename
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Histogram saved to: {output_path}")
        return str(output_path)
    
    def print_summary(self):
        """Print analysis summary to console."""
        if not self.stats:
            logger.warning("No statistics available!")
            return
        
        print("\n" + "="*60)
        print("cfDNA Fragmentomics Analysis Summary")
        print("="*60)
        print(f"BAM File: {self.bam_path}")
        print(f"Total Reads Processed: {self.stats['total_reads_processed']:,}")
        print(f"Valid Fragment Pairs: {self.stats['filtered_read_pairs']:,}")
        print(f"Filtering Rate: {self.stats['filtered_read_pairs']/self.stats['total_reads_processed']*100:.2f}%")
        print("\nFragment Length Statistics:")
        print(f"  Mean: {self.stats['mean_fragment_length']:.1f} bp")
        print(f"  Median: {self.stats['median_fragment_length']:.1f} bp")
        print(f"  Std Dev: {self.stats['std_fragment_length']:.1f} bp")
        print(f"  Range: {self.stats['min_fragment_length']} - {self.stats['max_fragment_length']} bp")
        print(f"  25th Percentile: {self.stats['fragment_length_25th_percentile']:.1f} bp")
        print(f"  75th Percentile: {self.stats['fragment_length_75th_percentile']:.1f} bp")
        print("="*60)

    @staticmethod
    def calculate_stats_from_list(fragment_lengths):
        if not fragment_lengths:
            return {}
        return {
            'mean_fragment_length': float(np.mean(fragment_lengths)),
            'median_fragment_length': float(np.median(fragment_lengths)),
            'std_fragment_length': float(np.std(fragment_lengths)),
            'min_fragment_length': int(np.min(fragment_lengths)),
            'max_fragment_length': int(np.max(fragment_lengths)),
            'fragment_length_25th_percentile': float(np.percentile(fragment_lengths, 25)),
            'fragment_length_75th_percentile': float(np.percentile(fragment_lengths, 75))
        }


def main():
    """Main function to run the cfDNA fragment analysis."""
    parser = argparse.ArgumentParser(
        description="Analyze cfDNA fragment lengths from BAM files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python parse_cfDNA_frag_lengths.py sample.bam
  python parse_cfDNA_frag_lengths.py sample.bam --max-reads 10000 --output-dir results
  python parse_cfDNA_frag_lengths.py sample.bam --no-plot
        """
    )
    
    parser.add_argument(
        'bam_file',
        help='Path to input BAM file'
    )
    
    parser.add_argument(
        '--output-dir',
        default='outputs',
        help='Output directory (default: outputs)'
    )
    
    parser.add_argument(
        '--max-reads',
        type=int,
        help='Maximum number of read pairs to process (for testing)'
    )
    
    parser.add_argument(
        '--no-plot',
        action='store_true',
        help='Skip creating histogram plot'
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
        # Initialize analyzer
        analyzer = cfDNAFragmentAnalyzer(args.bam_file, args.output_dir)
        
        # Run analysis
        stats = analyzer.analyze_fragments(max_reads=args.max_reads)
        
        if stats:
            # Save results
            analyzer.save_fragment_lengths_csv()
            
            # Create plot
            if not args.no_plot:
                analyzer.create_histogram()
            
            # Print summary
            analyzer.print_summary()
        else:
            logger.error("Analysis failed - no valid fragments found!")
            sys.exit(1)
            
    except Exception as e:
        logger.error(f"Error during analysis: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main() 