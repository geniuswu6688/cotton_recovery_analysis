# Cotton Genetic Background Recovery Analysis
This repository contains scripts and results for SNP-based genetic background recovery analysis in cotton.

## Workflow
VCF → SNP filtering → Genotype matrix construction → Recovery rate calculation → Visualization

## Project Structure
scripts/        # analysis scripts
data/           # processed data
results/        # figures and tables

## Scripts
* `high_quality_snps.py`: SNP filtering
* `samples_snp_matrix.py`: genotype matrix construction
* `average_recovery_rank.py`: recovery rate calculation
* `chrom_recovery.py`: chromosome-level analysis
* `sample_sliding.py`: sliding window analysis
* `introgression_map.py`: introgression visualization

## Data
* `high_quality_snps.csv`: SNP filtering
* `samples_snp_matrix.csv`: genotype matrix construction
* `average_recovery_rank.csv`: recovery rate
* `chrom_recovery_matrix.csv`: chromosome-level analysis
* `sample_sliding.csv`: sliding window analysis

## Results
* `average_recovery_rank.png`: recovery rate
* `chrom_recovery_heatmap.png`: chromosome-level analysis
* `sample_sliding_joint.png`: sliding window analysis
* `homo_curve`×22: sliding window analysis
* `introgression_map.png`: introgression visualization

## Requirements
* Python 3.13
* pandas
* numpy
* matplotlib

## Notes
Large files (e.g., VCF) are not included due to size limitations.

## Author
Wu Wanxin  
Undergraduate Thesis Project  
Zhejiang University  
2026