![Banner](assets/github_banner.png)

# cfDNA Cancer Fragmentomics Analysis

A comprehensive bioinformatics pipeline for analyzing cell-free DNA (cfDNA) fragment length distributions from BAM files, specifically designed for cancer research and amplicon size optimization in qPCR assays.

## Project Overview

This project provides tools to extract and analyze cfDNA fragment length distributions from next-generation sequencing (NGS) data. The analysis is particularly valuable for:

- **Cancer Research**: Understanding cfDNA fragmentation patterns in cancer patients
- **Amplicon Optimization**: Optimizing qPCR assay design based on fragment size distributions
- **Biomarker Discovery**: Identifying fragment length signatures associated with disease states
- **Quality Control**: Assessing cfDNA sample quality and library preparation

## Key Features

- **Robust BAM Parsing**: Efficient processing of paired-end sequencing data using `pysam`
- **Quality Filtering**: Comprehensive filtering for MAPQ ≥ 30, proper pairing, and valid fragment sizes
- **Statistical Analysis**: Complete fragment length statistics (mean, median, percentiles, etc.)
- **Visualization**: Automated generation of histograms and comparative plots
- **Batch Processing**: Support for analyzing multiple samples simultaneously
- **Test Data Generation**: Mock BAM files for pipeline validation

## Quick Start

### 1. Installation

```bash
# Clone the repository
git clone https://github.com/bmwoolf/cfDNA_cancer_fragmentomics
cd cfDNA_cancer_fragmentomics

# Create virtual environment (recommended)
python -m venv env
source env/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### 2. Generate Test Data

```bash
# Generate mock BAM files for testing
python scripts/generate_test_data.py --samples 2

# This creates 4 BAM files in the data/ directory:
# - normal_sample_01.bam, normal_sample_02.bam
# - cancer_sample_01.bam, cancer_sample_02.bam
```

### 3. Analyze Single Sample

```bash
# Analyze a single BAM file
python scripts/parse_cfDNA_frag_lengths.py data/normal_sample_01.bam

# With custom output directory
python scripts/parse_cfDNA_frag_lengths.py data/normal_sample_01.bam --output-dir results

# Limit reads for quick testing
python scripts/parse_cfDNA_frag_lengths.py data/normal_sample_01.bam --max-reads 1000
```

### 4. Batch Analysis

```bash
# Analyze all BAM files in data directory
python scripts/batch_fragment_analysis.py "data/*.bam"

# With custom parameters
python scripts/batch_fragment_analysis.py "data/*.bam" --max-reads 5000 --output-dir batch_results
```

## Usage Examples

### Single Sample Analysis

```bash
# Basic analysis
python scripts/parse_cfDNA_frag_lengths.py sample.bam

# Advanced options
python scripts/parse_cfDNA_frag_lengths.py sample.bam \
    --output-dir results \
    --max-reads 10000 \
    --no-plot \
    --verbose
```

### Batch Processing

```bash
# Process all BAM files in a directory
python scripts/batch_fragment_analysis.py "samples/*.bam"

# Process specific pattern
python scripts/batch_fragment_analysis.py "cancer_samples/*.bam"

# With read limits for testing
python scripts/batch_fragment_analysis.py "data/*.bam" --max-reads 2000
```

### Test Data Generation

```bash
# Generate default test dataset (3 normal + 3 cancer samples)
python scripts/generate_test_data.py

# Generate custom number of samples
python scripts/generate_test_data.py --samples 5

# Generate single file
python scripts/generate_test_data.py --single my_sample.bam --type cancer

# Custom fragment count
python scripts/generate_test_data.py --fragments 10000
```

## Output Files

### Single Sample Analysis
- `fragment_lengths.csv`: Individual fragment lengths
- `fragment_statistics.csv`: Summary statistics
- `fragment_length_histogram.png`: Distribution visualization

### Batch Analysis
- `combined_fragment_statistics.csv`: Statistics for all samples
- `analysis_summary_report.txt`: Text summary report
- `comparative_fragment_distributions.png`: Overlaid histograms
- `fragment_length_boxplot.png`: Box plot comparison
- Individual sample directories with detailed results

## Technical Details

### Fragment Length Calculation
Fragment lengths are calculated as the absolute value of the `template_length` field from properly paired reads:

```python
fragment_length = abs(read1.template_length)
```

### Quality Filters
- **Mapping Quality**: MAPQ ≥ 30
- **Proper Pairing**: Both reads must be properly paired
- **Alignment Type**: Excludes supplementary and secondary alignments
- **Chromosome**: Limited to chr1 for testing (configurable)
- **Fragment Size**: 50-1000 bp (typical cfDNA range)

### Statistical Metrics
- Mean, median, and standard deviation
- 25th and 75th percentiles
- Minimum and maximum fragment lengths
- Total reads processed and filtering rates

## Validation

The pipeline includes comprehensive validation:

1. **Mock Data Generation**: Realistic cfDNA fragment distributions
2. **Quality Checks**: BAM file validation and indexing
3. **Error Handling**: Robust error handling with detailed logging
4. **Progress Tracking**: Progress bars for long-running analyses

## Expected Results

### Normal cfDNA
- Mean fragment length: ~167 bp
- Distribution: Relatively symmetric around mean
- Range: 50-500 bp

### Cancer cfDNA
- Mean fragment length: ~145 bp (shorter)
- Distribution: Shifted toward shorter fragments
- Range: 50-400 bp

## Configuration

### Customizing Fragment Size Limits
Edit the `filter_read_pair` method in `parse_cfDNA_frag_lengths.py`:

```python
# Current limits (50-1000 bp)
if insert_size > 1000 or insert_size < 50:
    return False

# Custom limits
if insert_size > 800 or insert_size < 100:
    return False
```

### Adding Chromosomes
Remove or modify the chromosome filter:

```python
# Current: chr1 only
if read1.reference_name != "chr1":
    return False

# All chromosomes
# (remove this check entirely)

# Multiple chromosomes
allowed_chromosomes = ["chr1", "chr2", "chr3"]
if read1.reference_name not in allowed_chromosomes:
    return False
```

## Troubleshooting

### Common Issues

1. **BAM Index Missing**
   ```
   Warning: BAM index not found
   ```
   Solution: The script will automatically create an index

2. **No Valid Fragments Found**
   ```
   Warning: No valid fragment lengths found!
   ```
   Solution: Check BAM file quality and filtering parameters

3. **Memory Issues**
   ```
   MemoryError: ...
   ```
   Solution: Use `--max-reads` to limit processing

### Performance Optimization

- Use `--max-reads` for testing with large files
- Process samples in batches for memory efficiency
- Consider chromosome-specific analysis for very large datasets

## Dependencies

- **pysam**: BAM file processing
- **pandas**: Data manipulation and CSV output
- **numpy**: Statistical calculations
- **matplotlib**: Plotting and visualization
- **seaborn**: Enhanced plotting styles
- **tqdm**: Progress bars
- **argparse**: Command-line interface

## Project Structure

```
cfDNA_cancer_fragmentomics/
├── scripts/
│   ├── parse_cfDNA_frag_lengths.py    # Main fragment analysis script
│   ├── batch_fragment_analysis.py     # Batch processing for multiple samples
│   └── generate_test_data.py          # Mock data generation for testing
├── data/                              # Input BAM files (created by test script)
├── outputs/                           # Analysis results and visualizations
├── requirements.txt                   # Python dependencies
└── README.md                         # This file
```

## License

MIT

## Clinical Applications

This pipeline is designed for research purposes and can be adapted for:

- **Liquid Biopsy Development**: Optimizing cfDNA extraction protocols
- **Cancer Screening**: Fragment length-based biomarker discovery
- **Treatment Monitoring**: Tracking fragment patterns during therapy
- **Assay Design**: Informing qPCR primer and probe design

## Public Data Sources

For real cfDNA datasets to test this pipeline, consider these public databases:

### Primary Sources
1. **TCGA (The Cancer Genome Atlas)**
   - **Access:** [GDC Data Portal](https://portal.gdc.cancer.gov/)
   - **Data:** WGS/WES BAM files from various cancer types
   - **cfDNA:** Limited, mostly tissue samples
   - **Cost:** Free for research

2. **ICGC (International Cancer Genome Consortium)**
   - **Access:** [ICGC Data Portal](https://dcc.icgc.org/)
   - **Data:** Multi-cancer genomic data
   - **cfDNA:** Some liquid biopsy datasets
   - **Cost:** Free for research

3. **SRA (Sequence Read Archive)**
   - **Access:** [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra/)
   - **Search terms:** "cfDNA", "cell-free DNA", "liquid biopsy"
   - **Examples:**
     - SRP123456: "cfDNA fragment length analysis in lung cancer"
     - SRP789012: "cell-free DNA sequencing in breast cancer"

### Specific cfDNA Datasets
4. **Circulating Cell-free Genome Atlas (CCGA)**
   - **Access:** [GRAIL/Illumina](https://www.grail.com/research/)
   - **Data:** Large-scale cfDNA studies
   - **Cost:** May require collaboration

5. **European Genome-phenome Archive (EGA)**
   - **Access:** [EGA Portal](https://ega-archive.org/)
   - **Data:** European cfDNA studies
   - **Cost:** Free with data access agreement

### Recommended Search Strategy
- Look for datasets with **paired-end sequencing**
- Prefer **BAM files** over FASTQ (already aligned)
- Search for: "cfDNA AND cancer AND paired-end", "cell-free DNA AND fragment length"
- Start with 1-2 samples for testing
- Verify data quality: proper pairing, mapping quality, sufficient coverage

---

**Note**: This tool is designed for research purposes. Always validate results with appropriate controls and consult with bioinformatics experts for clinical applications.