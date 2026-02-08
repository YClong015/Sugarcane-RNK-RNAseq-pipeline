# Sugarcane RKN RNA-seq pipeline (R570 reference; RKN-only sample list)

Original scripts are preserved verbatim for provenance:
- `scripts/original/`
- `slurm/original/`

## Setup

```bash
bash scripts/99_copy_config_example.sh
# edit config/config.env (set BASE_DIR and reference paths)
```

## Prepare the sample list (required)

1) Create/edit the file:

- `08_deseq2_R570/samples_rkn_only.txt`

You can initialize it from the template and sanitize it:

```bash
bash scripts/00_init_dirs.sh
bash scripts/02_prepare_samples_rkn_only.sh
```

## Run (step-by-step)

```bash
bash scripts/03_merge_lanes_rkn_only.sh
bash scripts/04_qc_raw.sh
bash scripts/05_trim_fastp.sh
bash scripts/06_qc_clean.sh
bash scripts/07_build_bowtie2_index.sh
```

### Smut removal (Slurm array)

Compute N (number of samples) and submit with an explicit array range:

```bash
N=$(wc -l < 08_deseq2_R570/samples_rkn_only.txt)
mkdir -p slurm_logs
sbatch --array=1-"$N" slurm/rm_smut_array.sbatch
```

### Salmon

```bash
bash scripts/08_build_salmon_index.sh
bash scripts/09_salmon_quant_rkn_only.sh
```

### tx2gene + metadata

```bash
bash scripts/10_make_tx2gene_from_gff3.sh
bash scripts/11_make_sample_metadata_from_list.sh
```

## Downstream R analysis

You can run your original R scripts under `scripts/original/` using the
paths configured in `config/config.env`. The modular wrappers are intentionally
minimal; keep your exact analysis logic in the original scripts if you want
bit-for-bit reproduction.

## Notes

- If `fastp` rejects `--cut_window_size`, replace with the compatible option
  for your installed fastp version.
