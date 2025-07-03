#!/usr/bin/env python3
"""
Generate Test Data for cfDNA Fragmentomics Analysis

This script creates mock BAM files with simulated cfDNA fragment data
for testing the fragmentomics analysis pipeline.
"""

import argparse
import sys
import os
from pathlib import Path
import random
import logging

import pysam
import numpy as np

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class MockBAMGenerator:
    """Generates mock BAM files with simulated cfDNA data."""
    
    def __init__(self, output_dir: str = "data"):
        """
        Initialize the mock data generator.
        
        Args:
            output_dir: Directory to save generated BAM files
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # cfDNA fragment length distribution parameters
        # Based on typical cfDNA characteristics
        self.fragment_params = {
            'normal': {
                'mean': 167,  # Typical cfDNA fragment size
                'std': 50,
                'min': 50,
                'max': 500
            },
            'cancer': {
                'mean': 145,  # Shorter fragments in cancer
                'std': 45,
                'min': 50,
                'max': 400
            }
        }
    
    def generate_fragment_lengths(self, num_fragments: int, sample_type: str = 'normal') -> list:
        """
        Generate fragment lengths following cfDNA distribution.
        
        Args:
            num_fragments: Number of fragments to generate
            sample_type: 'normal' or 'cancer'
            
        Returns:
            list: Fragment lengths
        """
        params = self.fragment_params[sample_type]
        
        # Generate lengths using truncated normal distribution
        lengths = []
        while len(lengths) < num_fragments:
            # Generate from normal distribution
            new_lengths = np.random.normal(params['mean'], params['std'], num_fragments * 2)
            
            # Filter to valid range
            valid_lengths = new_lengths[
                (new_lengths >= params['min']) & 
                (new_lengths <= params['max'])
            ]
            
            lengths.extend(valid_lengths[:num_fragments - len(lengths)])
        
        return [int(length) for length in lengths[:num_fragments]]
    
    def create_mock_bam(self, filename: str, num_fragments: int = 10000, 
                       sample_type: str = 'normal', chromosome: str = 'chr1') -> str:
        """
        Create a mock BAM file with simulated cfDNA data.
        
        Args:
            filename: Output BAM filename
            num_fragments: Number of fragment pairs to generate
            sample_type: 'normal' or 'cancer'
            chromosome: Chromosome to simulate
            
        Returns:
            str: Path to created BAM file
        """
        output_path = self.output_dir / filename
        
        # Generate fragment lengths
        fragment_lengths = self.generate_fragment_lengths(num_fragments, sample_type)
        
        # Create header
        header = {
            'HD': {'VN': '1.6', 'SO': 'coordinate'},
            'SQ': [
                {'LN': 248956422, 'SN': 'chr1'},  # Human chr1 length
                {'LN': 242193529, 'SN': 'chr2'},
                {'LN': 198295559, 'SN': 'chr3'},
                {'LN': 190214555, 'SN': 'chr4'},
                {'LN': 181538259, 'SN': 'chr5'},
                {'LN': 170805979, 'SN': 'chr6'},
                {'LN': 159345973, 'SN': 'chr7'},
                {'LN': 145138636, 'SN': 'chr8'},
                {'LN': 138394717, 'SN': 'chr9'},
                {'LN': 133797422, 'SN': 'chr10'},
                {'LN': 135086622, 'SN': 'chr11'},
                {'LN': 133275309, 'SN': 'chr12'},
                {'LN': 114364328, 'SN': 'chr13'},
                {'LN': 107043718, 'SN': 'chr14'},
                {'LN': 101991189, 'SN': 'chr15'},
                {'LN': 90338345, 'SN': 'chr16'},
                {'LN': 83257441, 'SN': 'chr17'},
                {'LN': 80373285, 'SN': 'chr18'},
                {'LN': 58617616, 'SN': 'chr19'},
                {'LN': 64444167, 'SN': 'chr20'},
                {'LN': 46709983, 'SN': 'chr21'},
                {'LN': 50818468, 'SN': 'chr22'},
                {'LN': 156040895, 'SN': 'chrX'},
                {'LN': 57227415, 'SN': 'chrY'}
            ]
        }
        
        # Create temporary unsorted BAM file
        temp_bam_path = str(output_path) + ".unsorted"
        with pysam.AlignmentFile(temp_bam_path, "wb", header=header) as bam:
            read_id = 1
            
            for fragment_length in fragment_lengths:
                # Random start position on chromosome
                start_pos = random.randint(1000000, 10000000)  # Avoid telomeres
                
                # Generate read pair
                read1 = self._create_read_pair(
                    bam, read_id, chromosome, start_pos, fragment_length, is_first=True
                )
                read2 = self._create_read_pair(
                    bam, read_id, chromosome, start_pos, fragment_length, is_first=False
                )
                
                # Write reads
                bam.write(read1)
                bam.write(read2)
                
                read_id += 1
        
        # Sort the BAM file
        pysam.sort("-o", str(output_path), temp_bam_path)
        
        # Remove temporary unsorted file
        os.remove(temp_bam_path)
        
        # Create index
        pysam.index(str(output_path))
        
        logger.info(f"Created mock BAM file: {output_path}")
        logger.info(f"  Sample type: {sample_type}")
        logger.info(f"  Fragments: {num_fragments}")
        logger.info(f"  Mean fragment length: {np.mean(fragment_lengths):.1f} bp")
        
        return str(output_path)
    
    def _create_read_pair(self, bam, read_id: int, chromosome: str, start_pos: int, 
                         fragment_length: int, is_first: bool) -> pysam.AlignedSegment:
        """
        Create a single read in a pair.
        
        Args:
            bam: BAM file object
            read_id: Read pair ID
            chromosome: Chromosome name
            start_pos: Start position
            fragment_length: Fragment length
            is_first: Whether this is the first read in pair
            
        Returns:
            pysam.AlignedSegment: Read object
        """
        read = pysam.AlignedSegment()
        
        # Basic read properties
        read.query_name = f"READ_{read_id:08d}"
        read.reference_id = bam.get_tid(chromosome)
        read.reference_start = start_pos if is_first else start_pos + fragment_length - 150
        read.mapping_quality = 60  # High quality
        
        # Read sequence (mock)
        read.query_sequence = "A" * 150
        read.query_qualities = pysam.qualitystring_to_array("I" * 150)  # High quality
        
        # CIGAR string (150M = 150 matches)
        read.cigarstring = "150M"
        
        # Template length
        if is_first:
            read.template_length = fragment_length
        else:
            read.template_length = -fragment_length
        
        # Flags
        if is_first:
            read.flag = 99  # Proper pair, first in pair, forward strand
        else:
            read.flag = 147  # Proper pair, second in pair, reverse strand
        
        # Tags
        read.set_tag("NM", 0)  # Edit distance
        read.set_tag("AS", 150)  # Alignment score
        read.set_tag("XS", 0)  # Suboptimal alignment score
        
        return read
    
    def generate_test_dataset(self, num_samples: int = 3) -> list:
        """
        Generate a complete test dataset with multiple samples.
        
        Args:
            num_samples: Number of samples to generate
            
        Returns:
            list: Paths to generated BAM files
        """
        generated_files = []
        
        for i in range(num_samples):
            # Generate normal sample
            normal_filename = f"normal_sample_{i+1:02d}.bam"
            normal_path = self.create_mock_bam(
                normal_filename, 
                num_fragments=5000, 
                sample_type='normal'
            )
            generated_files.append(normal_path)
            
            # Generate cancer sample
            cancer_filename = f"cancer_sample_{i+1:02d}.bam"
            cancer_path = self.create_mock_bam(
                cancer_filename, 
                num_fragments=5000, 
                sample_type='cancer'
            )
            generated_files.append(cancer_path)
        
        logger.info(f"Generated {len(generated_files)} test BAM files")
        return generated_files


def main():
    """Main function to generate test data."""
    parser = argparse.ArgumentParser(
        description="Generate mock BAM files for cfDNA fragmentomics testing",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python generate_test_data.py --samples 2
  python generate_test_data.py --output-dir test_data --fragments 1000
  python generate_test_data.py --single normal_sample.bam --type normal
        """
    )
    
    parser.add_argument(
        '--output-dir',
        default='data',
        help='Output directory for BAM files (default: data)'
    )
    
    parser.add_argument(
        '--samples',
        type=int,
        default=3,
        help='Number of sample pairs to generate (default: 3)'
    )
    
    parser.add_argument(
        '--fragments',
        type=int,
        default=5000,
        help='Number of fragments per sample (default: 5000)'
    )
    
    parser.add_argument(
        '--single',
        help='Generate single BAM file with specified name'
    )
    
    parser.add_argument(
        '--type',
        choices=['normal', 'cancer'],
        default='normal',
        help='Sample type for single file generation (default: normal)'
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
        # Initialize generator
        generator = MockBAMGenerator(args.output_dir)
        
        if args.single:
            # Generate single file
            generator.create_mock_bam(
                args.single, 
                num_fragments=args.fragments, 
                sample_type=args.type
            )
        else:
            # Generate test dataset
            generator.generate_test_dataset(args.samples)
        
        logger.info("Test data generation completed successfully!")
        
    except Exception as e:
        logger.error(f"Error generating test data: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main() 