# Sugarcane–RKN RNA-seq Pipeline (R570 Reference)

A reproducible, HPC-friendly bioinformatics pipeline for analyzing **Sugarcane (*Saccharum* spp. hybrid)** RNA-seq
data focused on **Root-Knot Nematode (RKN)** interactions using the **R570 reference**.

A key feature is a dedicated **smut decontamination step**: reads from *Sporisorium scitamineum*
(sugarcane smut pathogen) are removed **before** quantification to reduce fungal noise and improve accuracy
in plant–nematode expression signals.

---

## 🧬 Project Overview

This workflow is designed to run on SLURM-based HPC clusters and supports parallel processing via array jobs.
It includes a fully defined Conda environment and an integrated R workflow for differential expression and
functional enrichment.

### Key Features

- **HPC-Optimized**: Designed for SLURM-based clusters using array jobs for efficient parallel processing.
- **Reproducible Environment**: Fully defined Conda environment (`environment.yml`).
- **Smut Decontamination**: Filters out *S. scitamineum* reads using Bowtie2 alignment.
- **Accurate Quantification**: Uses **Salmon** for transcript-level quantification against the R570 genome.
- **Integrated R Workflow**: Includes scripts for **DESeq2** (Differential Expression) and
  **ClusterProfiler** (Enrichment).

---

## 📊 Pipeline Workflow

```mermaid
graph TD
    subgraph Pre-processing
    Raw[Raw FastQ] --> QC1[FastQC / MultiQC]
    QC1 --> Trim[Fastp: Trimming & Adapter Removal]
    Trim --> QC2[Clean Data QC]
    end

    subgraph Decontamination [Key Step: Smut Removal]
    QC2 --> Filter{Bowtie2: Smut Index}
    Filter -- Aligned Reads --> Smut[Discard Smut Contamination]
    Filter -- Unaligned Reads --> Clean[Clean Sugarcane Reads]
    end

    subgraph Quantification
    Clean --> Salmon[Salmon Quant (R570)]
    Salmon --> Counts[Quant.sf / Count Matrix]
    end

    subgraph Analysis
    Counts --> R[R: DESeq2 Analysis]
    R --> DEGs[Differentially Expressed Genes]
    DEGs --> Enrich[GO/KEGG Enrichment]
    end
```

---

## 📂 Directory Structure

```text
.
├── config/
│   └── config.env                 # Global configuration (paths to R570/Smut genomes, output dirs)
├── scripts/
│   ├── 00_init_dirs.sh            # Initialize directory structure
│   ├── 02_prepare_samples...      # Generate sample ID lists
│   ├── 03_merge_lanes...          # Merge lane files (if applicable)
│   ├── 04-06_qc_trim...           # QC and Fastp trimming wrappers
│   ├── 07_build_bowtie2...        # Build Smut genome index
│   ├── 08_build_salmon...         # Build R570 transcriptome index
│   ├── 09_salmon_quant...         # Quantification wrapper
│   └── ...                        # Helper scripts for metadata/tx2gene
├── slurm/
│   └── rm_smut_array.sbatch       # SLURM array script for Smut removal
├── 08_deseq2_R570/
│   └── samples_rkn_only.txt       # Input sample ID list
├── metadata.csv                   # Experimental design for R analysis
└── environment.yml                # Conda environment specification
```

---

## 🛠️ Installation & Setup

### 1) Clone Repository

```bash
git clone https://github.com/YClong015/Sugarcane-RNK-RNAseq-pipeline.git
cd Sugarcane-RNK-RNAseq-pipeline
```

### 2) Create Conda Environment

Ensure `conda` (or `mamba`) is installed.

```bash
conda env create -f environment.yml
conda activate sugarcane-rkn
```

Dependencies typically include: `fastqc`, `multiqc`, `fastp`, `bowtie2`, `salmon`, `samtools`,
`python=3.10`, `R=4.4`.

### 3) Configuration

Copy the template and edit `config.env` to set file paths (genomes, indices, output directories).

```bash
bash scripts/99_copy_config_example.sh
nano config/config.env
```

---

## 📝 Input Preparation

### 1) Sample List (`samples_rkn_only.txt`)

Create a plain text file listing sample IDs (one per line). IDs must match FASTQ prefixes used by the scripts.

**File location:** `08_deseq2_R570/samples_rkn_only.txt`


### 2) Metadata File (`metadata.csv`)

Create `metadata.csv` in the project root. This defines the experimental design for DESeq2.

Example:

```csv
sample_id,genotype,time_point,treatment,replicate
Q208_12w_C_1,Q208,12w,Control,1
Q208_12w_C_2,Q208,12w,Control,2
Q208_12w_RKN_1,Q208,12w,RKN,1
Q208_7d_C_2,Q208,1w,Control,2
```

---

## 🚀 Execution Guide

Run scripts in numerical order.

### Phase 1: QC & Trimming

```bash
# Setup directories & verify samples
bash scripts/00_init_dirs.sh
bash scripts/02_prepare_samples_rkn_only.sh

# Run QC and Trimming
bash scripts/04_qc_raw.sh
bash scripts/05_trim_fastp.sh
bash scripts/06_qc_clean.sh
```

### Phase 2: Smut Removal (Decontamination)

Reads are mapped against the *Sporisorium scitamineum* genome using Bowtie2.
Unmapped reads (Sugarcane/RKN) are preserved for downstream analysis.

```bash
# 1) Build Bowtie2 index for Smut genome
bash scripts/07_build_bowtie2_index.sh

# 2) Submit SLURM Array Job (array size N computed from sample list)
N=$(wc -l < 08_deseq2_R570/samples_rkn_only.txt)
mkdir -p slurm_logs
sbatch --array=1-"$N" slurm/rm_smut_array.sbatch
```

### Phase 3: Quantification

Quantify cleaned reads against the Sugarcane R570 transcriptome using Salmon.

```bash
# 1) Build Salmon index for R570 Transcriptome
bash scripts/08_build_salmon_index.sh

# 2) Run Quantification
bash scripts/09_salmon_quant_rkn_only.sh
```

### Phase 4: Prepare R Inputs

```bash
# Create transcript-to-gene map from GFF3
bash scripts/10_make_tx2gene_from_gff3.sh

# Generate R metadata object (if applicable)
bash scripts/11_make_sample_metadata_from_list.sh
```

---

## 📊 Downstream Analysis (R)

The R analysis workflow performs:

- **Data Import**: Load Salmon `quant.sf` using `tximport`
- **PCA**: Visualize sample clustering and detect outliers
- **DESeq2**: Differential expression analysis (e.g., RKN vs Control) across genotypes/time points
- **Enrichment**: GO and KEGG enrichment using `clusterProfiler`

---

## 📌 Typical Outputs

- QC reports: FastQC results + MultiQC summary
- Trimmed reads (fastp outputs)
- Decontaminated reads (smut-aligned discarded; unmapped retained)
- Salmon quantification (`quant.sf`) and derived count matrices
- DESeq2 results tables (log2FC, padj), PCA plots, volcano plots
- GO/KEGG enrichment tables and visualizations

---

## 📜 Citation & Credits

Maintainer: **Yanchen Zheng (Ethan)**  
Affiliation: **The University of Queensland**

If you use this pipeline, please cite the underlying tools:

- **Fastp**: Chen et al., 2018
- **Bowtie2**: Langmead et al., 2012
- **Salmon**: Patro et al., 2017
- **DESeq2**: Love et al., 2014
- **ClusterProfiler**: Wu et al., 2021
