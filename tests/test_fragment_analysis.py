#!/usr/bin/env python3
"""
Unit Tests for cfDNA Fragmentomics Analysis

Tests for the fragment analysis pipeline including BAM parsing,
filtering, statistics calculation, and visualization.
"""

import unittest
import tempfile
import os
import shutil
from pathlib import Path
import pandas as pd
import numpy as np

# Add src directory to path for imports
import sys
sys.path.append(str(Path(__file__).parent.parent / 'src'))

from parse_cfDNA_frag_lengths import cfDNAFragmentAnalyzer
from generate_test_data import MockBAMGenerator


class TestMockBAMGenerator(unittest.TestCase):
    """Test the mock BAM file generator."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.generator = MockBAMGenerator(self.temp_dir)
    
    def tearDown(self):
        """Clean up test environment."""
        shutil.rmtree(self.temp_dir)
    
    def test_fragment_length_generation_normal(self):
        """Test normal cfDNA fragment length generation."""
        lengths = self.generator.generate_fragment_lengths(1000, 'normal')
        
        self.assertEqual(len(lengths), 1000)
        self.assertTrue(all(50 <= length <= 500 for length in lengths))
        
        # Check distribution characteristics
        mean_length = np.mean(lengths)
        self.assertTrue(150 <= mean_length <= 185)  # Expected range
    
    def test_fragment_length_generation_cancer(self):
        """Test cancer cfDNA fragment length generation."""
        lengths = self.generator.generate_fragment_lengths(1000, 'cancer')
        
        self.assertEqual(len(lengths), 1000)
        self.assertTrue(all(50 <= length <= 400 for length in lengths))
        
        # Check distribution characteristics
        mean_length = np.mean(lengths)
        self.assertTrue(130 <= mean_length <= 160)  # Expected range
    
    def test_mock_bam_creation(self):
        """Test mock BAM file creation."""
        bam_path = self.generator.create_mock_bam(
            "test.bam", 
            num_fragments=100, 
            sample_type='normal'
        )
        
        self.assertTrue(os.path.exists(bam_path))
        self.assertTrue(os.path.exists(bam_path + ".bai"))  # Index file
        
        # Test BAM file can be opened
        import pysam
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            read_count = 0
            for read in bam.fetch():
                read_count += 1
            
            # Should have 200 reads (100 pairs)
            self.assertEqual(read_count, 200)


class TestcfDNAFragmentAnalyzer(unittest.TestCase):
    """Test the main fragment analyzer."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.generator = MockBAMGenerator(self.temp_dir)
        
        # Create test BAM file
        self.test_bam = self.generator.create_mock_bam(
            "test_sample.bam",
            num_fragments=500,
            sample_type='normal'
        )
    
    def tearDown(self):
        """Clean up test environment."""
        shutil.rmtree(self.temp_dir)
    
    def test_analyzer_initialization(self):
        """Test analyzer initialization."""
        analyzer = cfDNAFragmentAnalyzer(self.test_bam, self.temp_dir)
        
        self.assertEqual(analyzer.bam_path, self.test_bam)
        self.assertEqual(str(analyzer.output_dir), self.temp_dir)
        self.assertEqual(len(analyzer.fragment_lengths), 0)
        self.assertEqual(len(analyzer.stats), 0)
    
    def test_analyzer_with_nonexistent_file(self):
        """Test analyzer with non-existent BAM file."""
        with self.assertRaises(FileNotFoundError):
            cfDNAFragmentAnalyzer("nonexistent.bam")
    
    def test_fragment_analysis(self):
        """Test complete fragment analysis."""
        analyzer = cfDNAFragmentAnalyzer(self.test_bam, self.temp_dir)
        stats = analyzer.analyze_fragments(max_reads=100)
        
        # Check that analysis produced results
        self.assertIsNotNone(stats)
        self.assertGreater(len(analyzer.fragment_lengths), 0)
        self.assertGreater(stats['filtered_read_pairs'], 0)
        
        # Check statistics structure
        required_keys = [
            'total_reads_processed', 'filtered_read_pairs',
            'mean_fragment_length', 'median_fragment_length',
            'std_fragment_length', 'min_fragment_length',
            'max_fragment_length'
        ]
        
        for key in required_keys:
            self.assertIn(key, stats)
    
    def test_fragment_length_calculation(self):
        """Test fragment length calculation."""
        analyzer = cfDNAFragmentAnalyzer(self.test_bam, self.temp_dir)
        
        # Test with mock read objects
        class MockRead:
            def __init__(self, template_length):
                self.template_length = template_length
        
        read1 = MockRead(200)
        read2 = MockRead(-200)
        
        length = analyzer.calculate_fragment_length(read1, read2)
        self.assertEqual(length, 200)
    
    def test_read_pair_filtering(self):
        """Test read pair filtering logic."""
        analyzer = cfDNAFragmentAnalyzer(self.test_bam, self.temp_dir)
        
        # Create mock reads for testing
        class MockRead:
            def __init__(self, **kwargs):
                for key, value in kwargs.items():
                    setattr(self, key, value)
        
        # Test valid read pair
        read1 = MockRead(
            is_proper_pair=True,
            mapping_quality=60,
            is_supplementary=False,
            is_secondary=False,
            reference_name="chr1",
            template_length=200
        )
        read2 = MockRead(
            is_proper_pair=True,
            mapping_quality=60,
            is_supplementary=False,
            is_secondary=False,
            reference_name="chr1",
            template_length=-200
        )
        
        self.assertTrue(analyzer.filter_read_pair(read1, read2))
        
        # Test invalid read pair (low mapping quality)
        read1.mapping_quality = 20
        self.assertFalse(analyzer.filter_read_pair(read1, read2))
        
        # Test invalid read pair (wrong chromosome)
        read1.mapping_quality = 60
        read1.reference_name = "chr2"
        self.assertFalse(analyzer.filter_read_pair(read1, read2))
        
        # Test invalid read pair (too long fragment)
        read1.reference_name = "chr1"
        read1.template_length = 1500
        self.assertFalse(analyzer.filter_read_pair(read1, read2))
    
    def test_csv_output(self):
        """Test CSV file output."""
        analyzer = cfDNAFragmentAnalyzer(self.test_bam, self.temp_dir)
        analyzer.analyze_fragments(max_reads=50)
        
        # Test fragment lengths CSV
        csv_path = analyzer.save_fragment_lengths_csv("test_fragments.csv")
        self.assertTrue(os.path.exists(csv_path))
        
        # Check CSV content
        df = pd.read_csv(csv_path)
        self.assertIn('fragment_length', df.columns)
        self.assertEqual(len(df), len(analyzer.fragment_lengths))
        
        # Test statistics CSV
        stats_path = os.path.join(self.temp_dir, "fragment_statistics.csv")
        self.assertTrue(os.path.exists(stats_path))
        
        stats_df = pd.read_csv(stats_path)
        self.assertGreater(len(stats_df), 0)
    
    def test_histogram_creation(self):
        """Test histogram creation."""
        analyzer = cfDNAFragmentAnalyzer(self.test_bam, self.temp_dir)
        analyzer.analyze_fragments(max_reads=50)
        
        # Test histogram creation
        plot_path = analyzer.create_histogram("test_histogram.png")
        self.assertTrue(os.path.exists(plot_path))
        
        # Check file size (should be non-zero)
        self.assertGreater(os.path.getsize(plot_path), 0)
    
    def test_empty_results_handling(self):
        """Test handling of empty results."""
        # Create analyzer but don't run analysis
        analyzer = cfDNAFragmentAnalyzer(self.test_bam, self.temp_dir)
        
        # Test CSV output with no data
        csv_path = analyzer.save_fragment_lengths_csv("empty.csv")
        self.assertEqual(csv_path, "")
        
        # Test histogram with no data
        plot_path = analyzer.create_histogram("empty.png")
        self.assertEqual(plot_path, "")
    
    def test_statistics_calculation(self):
        """Test statistics calculation accuracy."""
        analyzer = cfDNAFragmentAnalyzer(self.test_bam, self.temp_dir)
        
        # Add known fragment lengths
        test_lengths = [100, 150, 200, 250, 300]
        analyzer.fragment_lengths = test_lengths
        
        # Calculate statistics
        stats = analyzer.analyze_fragments(max_reads=5)
        
        # Check calculated values
        self.assertAlmostEqual(stats['mean_fragment_length'], 200.0, places=1)
        self.assertAlmostEqual(stats['median_fragment_length'], 200.0, places=1)
        self.assertEqual(stats['min_fragment_length'], 100)
        self.assertEqual(stats['max_fragment_length'], 300)
        self.assertAlmostEqual(stats['std_fragment_length'], 70.71, places=1)  # Approximate


class TestBatchAnalysis(unittest.TestCase):
    """Test batch analysis functionality."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.generator = MockBAMGenerator(self.temp_dir)
        
        # Create multiple test BAM files
        self.bam_files = []
        for i in range(3):
            bam_path = self.generator.create_mock_bam(
                f"sample_{i+1}.bam",
                num_fragments=100,
                sample_type='normal' if i % 2 == 0 else 'cancer'
            )
            self.bam_files.append(bam_path)
    
    def tearDown(self):
        """Clean up test environment."""
        shutil.rmtree(self.temp_dir)
    
    def test_batch_analyzer_initialization(self):
        """Test batch analyzer initialization."""
        from batch_fragment_analysis import BatchFragmentAnalyzer
        
        pattern = os.path.join(self.temp_dir, "*.bam")
        batch_analyzer = BatchFragmentAnalyzer(pattern, self.temp_dir)
        
        self.assertEqual(len(batch_analyzer.bam_files), 3)
        self.assertEqual(len(batch_analyzer.results), 0)
        self.assertEqual(len(batch_analyzer.combined_stats), 0)


def run_tests():
    """Run all tests."""
    # Create tests directory if it doesn't exist
    tests_dir = Path(__file__).parent
    tests_dir.mkdir(exist_ok=True)
    
    # Run tests
    unittest.main(verbosity=2)


if __name__ == "__main__":
    run_tests() 