# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

DIVIA is a Nextflow-based bioinformatics pipeline for pediatric leukemia diagnosis and classification. It processes RNA-seq data to detect fusions, variants, and gene expression patterns for BALL, TALL, and AML lineage classification.

## Development Commands

### Running the Pipeline

```bash
# Show help and available parameters
nextflow run DIVIA.nf --help

# Run with local execution profile
nextflow run DIVIA.nf -profile local

# Run with HPC cluster (LSF) profile
nextflow run DIVIA.nf -profile cluster

# Run with Singularity containers
nextflow run DIVIA.nf -profile singularity

# Specify work directory and output location
nextflow run DIVIA.nf -profile cluster -w ~/Nextflow_work --outdir ./results

# Resume a previous run
nextflow run DIVIA.nf -profile cluster -resume
```

### Configuration Management

```bash
# View current configuration
nextflow config DIVIA.nf

# View configuration for specific profile
nextflow config DIVIA.nf -profile cluster

# Validate pipeline syntax
nextflow run DIVIA.nf -syntax
```

### Container Development

```bash
# Build Singularity container from definition file
sudo singularity build divia.sif divia_arun.def

# Test container interactively
singularity shell divia.sif
```

## Architecture Overview

### Pipeline Structure

DIVIA implements a **modular directed acyclic graph (DAG)** with 42 interconnected processes organized into main branches:

1. **Input Processing**: `Zcat_MergeFastq` → consolidates multi-lane FASTQ files
2. **Quality Control**: `Trim_Galore` (optional)
3. **Alignment**: `STAR_Mapping` → central hub for downstream analyses
4. **Fusion Detection**: `Fusion_Catcher`, `STAR_Fusion`, `Arriba_Fusion` (ensemble approach)
5. **Gene Expression**: `RSEM` → expression matrices → ML classification
6. **Variant Calling**: GATK pipeline, VarDict (ITD detection), Pindel (structural variants)
7. **Specialized Analyses**: RNASeqCNV, RNAIndel, iAdmix

### Key Design Patterns

- **Conditional Execution**: Each analysis branch controlled by Y/N parameters in `nextflow.config`
- **Channel Architecture**: Named channels (`star_bam_ch`, `rsem_ch`, etc.) enable data flow between processes
- **Adapter Processes**: Lightweight format converters (e.g., `STARMapping_Adaptor_*`)
- **Memory Scaling**: Error handling with retry logic and increasing memory allocation
- **Output Organization**: Structured by sample and process type under `${params.outdir}/${params.project}/`

### Configuration System

The `nextflow.config` file contains:
- **39 Feature Toggles**: Enable/disable specific analyses (Select_*)
- **Reference Data Paths**: Tool-specific databases and indices
- **Resource Allocation**: Memory and CPU per process
- **Execution Profiles**: `local`, `cluster` (LSF), `singularity`

Example key parameters:
- `params.fastq_filelist`: TSV file with format `SampleName \t Fastq_Block \t Strandness`
- `params.Lineage_Type`: "BALL", "TALL", or "AML" for appropriate ML models
- `params.Select_*`: Y/N toggles for each analysis module

### Sample Input Format

Input is a TSV file where each row contains:
```
SampleName \t Fastq_Block \t Strandness
```

- **Fastq_Block**: Space-separated R1/R2 pairs, comma-separated for multiple lanes
- **Strandness**: Used for stranded RNA-seq library protocols

### Container Environment

The `divia_arun.def` Singularity definition provides:
- Ubuntu 20.04 base with R 4.2.2
- 50+ bioinformatics tools pre-installed
- 40+ R packages (Bioconductor/CRAN)
- Custom `Rscript` environment via AppRun

### Machine Learning Components

The pipeline includes pediatric leukemia classifiers:
- **RANK Classification**: Distance-based classification
- **Random Forest**: Ensemble classification
- **t-SNE Visualization**: Dimensionality reduction for visualization
- Lineage-specific training models stored as RDS files

## Important Notes

- This is a **production bioinformatics pipeline** requiring significant computational resources
- Designed for **pediatric leukemia research** at St. Jude Children's Research Hospital
- Uses **ensemble approaches** for fusion detection and variant calling
- **Memory requirements**: Processes range from 8GB to 100GB RAM
- **Reference data dependencies**: Requires multiple genome references, databases, and pre-trained ML models
- **HPC Integration**: Optimized for LSF cluster execution with appropriate queue selection