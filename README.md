# Single-Cell RNA-seq Analysis of Peripheral Blood in HR+ Breast Cancer Patients Treated with Immunotherapy

[![Data: GSE300475](https://img.shields.io/badge/data-GSE300475-blue)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE300475)
[![Platform: 10x Chromium](https://img.shields.io/badge/platform-10x%20Chromium-lightgrey)](https://www.10xgenomics.com/)
[![Framework: Seurat v5](https://img.shields.io/badge/framework-Seurat%20v5-orange)](https://satijalab.org/seurat/)
[![Data versioning: DVC](https://img.shields.io/badge/data%20versioning-DVC-green)](https://dvc.org/)
[![Status: Work in progress](https://img.shields.io/badge/status-work%20in%20progress-yellow)](#status--reproducibility)

> **Status: Work in progress.** This repository is under active development. Notebooks, scripts, and results are being iterated on and may change. See [Known limitations](#known-limitations--reproducibility).

---

## Table of Contents

- [Overview](#overview)
- [Study Summary](#study-summary)
- [Data & Results Versioning](#data--results-versioning)
- [Project Structure](#project-structure)
- [Installation](#installation)
- [Usage](#usage)
- [Workflow](#workflow)
- [Requirements](#requirements)
- [Known limitations & reproducibility](#known-limitations--reproducibility)
- [Citation](#citation)
- [Contact](#contact)

---

## Overview

This project analyzes a public single-cell RNA sequencing (scRNA-seq) dataset investigating peripheral immune dynamics in hormone receptor-positive (HR+) breast cancer patients treated with neoadjuvant **nab-paclitaxel + pembrolizumab**. The data combine **CITE-seq** (gene expression + hashtag/antibody capture) and **TCR-seq** (V(D)J) libraries generated on the 10x Genomics Chromium platform from longitudinal peripheral blood mononuclear cell (PBMC) samples across multiple patients and treatment timepoints.

The overall aim is to reproduce and extend the published analysis by performing sample demultiplexing, quality control, integration, clustering, and cell-type annotation, in order to explore peripheral blood correlates of immunotherapy response.

- **GEO Accession:** [GSE300475](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE300475)
- **BioProject:** [PRJNA1281033](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1281033)
- **Organism:** *Homo sapiens*
- **Platform:** Illumina NovaSeq 6000 (GPL24676)
- **Assay types:** 5' scRNA-seq, Feature Barcode (HTO), V(D)J (TCR)
- **Samples:** 32 GEO records covering multiple patients (PT1, PT5, PT6, PT7, PT11, PT13, PT15) across timepoints (week1, week3, week7, week15), each hashtag-multiplexed (CITE-seq HTO)
- **Original study contact:** Balko Lab, Vanderbilt University
- **Publication:** [PubMed 40610460](https://www.ncbi.nlm.nih.gov/pubmed/40610460)

## Study Summary

The limited clinical benefit of immune checkpoint inhibitors in breast cancer highlights the need for predictive biomarkers that can minimize risk and maximize benefit for patients. The original study performed single-cell RNA and TCR sequencing on PBMCs to monitor peripheral immune dynamics in an exploratory cohort of HR+ breast cancer patients undergoing neoadjuvant immunotherapy, aiming to identify candidate peripheral blood biomarkers of treatment response.

## Data & Results Versioning

Raw data, intermediate files, figures, and result tables are **not stored in this Git repository**. They are version-controlled with **DVC** and hosted on DagsHub:

🔗 **https://dagshub.com/mehranbeyki/scrna-immune-checkpoint-response**

To pull the tracked data/results locally, clone the repository and run the standard DVC pull workflow (`dvc pull`) after configuring the DagsHub remote — see the DagsHub project page for setup instructions.

The original raw sequencing/count files are also available directly from GEO:

- Supplementary raw data: `GSE300475_RAW.tar` (per-sample CSV/MTX/TSV files, 10x Genomics format)
- Feature reference: `GSE300475_feature_ref.xlsx`
- Raw FASTQ reads: via [SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA1281033)

## Project Structure

The actual on-disk layout of the repository (data files are not in Git — they are DVC-tracked or created at runtime):

```
.
├── README.md
├── environment.yml                  # Conda environment (R + Jupyter + Python helpers)
├── .gitignore                      # Ignores data/, *.rds, secrets, caches
├── .dvc/                            # DVC core config (remote = DagsHub)
├── docs/
│   └── report.md                   # Reserved for the final analysis report (placeholder)
├── notebooks/                       # R-kernel Jupyter notebooks (analysis pipeline)
│   ├── 01_data_download_preprocessing.ipynb   # GEO download, HTO demultiplexing, per-sample Seurat objects, merge
│   ├── 02_qc_clustering.R.ipynb              # QC, patient/timepoint mapping, SCTransform, Harmony integration, clustering
│   ├── 03_integration_clustering.ipynb        # Marker panels (DotPlot), FindAllMarkers, marker table export
│   ├── 04_celltype_annotation.ipynb          # Manual PBMC lineage annotation, T-cell exhaustion markers
│   ├── 05_exhaustion_quantification.ipynb    # Exhaustion scoring + non-parametric statistics (work in progress)
│   └── README.md                             # ASCII pipeline diagram
├── scripts/                        # Python helper scripts (raw-file organization)
│   ├── file_org.ipy                # Reorganize flat GEO archive into per-sample 10x folders (gex/ + vdj/)
│   └── check gzip features.ipy     # Diagnostic: inspect a gzipped features.tsv.gz file
├── src/
│   └── README.md                   # Reserved for refactored R/Python modules (placeholder)
└── results/                         # DVC pointers only — actual outputs live in DVC, not Git
    ├── figures.dvc                 # 19 figures tracked via DVC
    └── tables.dvc                  # 4 result tables tracked via DVC
```

At runtime, a `data/` directory is created locally (gitignored) holding:

```
data/
├── raw/
│   ├── GSE300475_RAW/              # Original files extracted from GSE300475_RAW.tar
│   ├── gex/                        # Per-sample gene-expression + HTO folders (organized by scripts/file_org.ipy)
│   │   └── scRNAseq PT1/           #   barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
│   └── vdj/                        # Per-sample V(D)J contig annotation files
├── mapping_table.csv               # GSM accession → sample title mapping (from GEOquery metadata)
└── processed/                     # Intermediate & final Seurat objects (.rds), tracked via DVC
```

> **Note on script extensions:** the helper scripts currently use the `.ipy` (IPython) extension, so they are run as `ipython scripts/file_org.ipy` rather than `python scripts/file_org.py`. Renaming to `.py` and refactoring into importable modules is on the roadmap.

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/mehranbeyki/scrna-immune-checkpoint-response.git
cd scrna-immune-checkpoint-response
```

### 2. Create the conda environment

```bash
conda env create -f environment.yml
conda activate scrna-immune-checkpoint-response
```

This installs R (≥ 4.5), the R kernel for Jupyter, **Seurat v5**, `harmony`, `GEOquery`, `glmGamPoi`, `dplyr`, `ggplot2`, `future`, `SingleR`, `celldex`, `RColorBrewer`, `patchwork`, plus Python + Jupyter and `pandas` for the helper scripts.

### 3. Register the R kernel

```bash
R -e 'IRkernel::installspec(name = "ir", displayname = "R")'
```

### 4. (Optional) Pull versioned data and results with DVC

```bash
pip install dvc
dvc remote modify origin url https://dagshub.com/mehranbeyki/scrna-immune-checkpoint-response.dvc
dvc pull
```

Alternatively, download the raw files directly from GEO (see [Data & Results Versioning](#data--results-versioning)) and place them under `data/raw/`.

## Usage

The notebooks embed hardcoded absolute paths during development. To run them portably, set the project root as an environment variable and run from the repository root:

```bash
export SCRNA_PROJ_ROOT="$PWD"
```

Then run the pipeline in order:

1. Organize the raw GEO archive into per-sample 10x folders:
   ```bash
   ipython scripts/file_org.ipy
   ```
2. Run the notebooks in sequence (R kernel):
   1. `notebooks/01_data_download_preprocessing.ipynb`
   2. `notebooks/02_qc_clustering.R.ipynb`
   3. `notebooks/03_integration_clustering.ipynb`
   4. `notebooks/04_celltype_annotation.ipynb`
   5. `notebooks/05_exhaustion_quantification.ipynb` *(work in progress)*

   ```bash
   jupyter lab
   ```

To run a notebook headless from the command line:

```bash
jupyter nbconvert --to notebook --execute notebooks/01_data_download_preprocessing.ipynb
```

## Workflow

The pipeline combines a Python preprocessing step with an R-based analysis pipeline (run as Jupyter notebooks with an R kernel).

### 1. Raw data organization (Python)

After downloading and extracting `GSE300475_RAW.tar` from GEO, `scripts/file_org.ipy` reorganizes the flat archive of per-GSM files into per-sample folders following the 10x Genomics convention (`barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz` for gene-expression samples, and `<sample>_all_contig_annotations.csv.gz` for V(D)J samples). Sample folder names are resolved via `data/mapping_table.csv` (GSM accession → sample title), which is generated from GEO metadata in notebook `01`.

```bash
ipython scripts/file_org.ipy
```

`scripts/check gzip features.ipy` is a small diagnostic utility used to inspect the contents of a sample's `features.tsv.gz` file (e.g., to confirm feature counts/formatting) before loading it into R.

### 2. Data download, metadata mapping & demultiplexing — `01_data_download_preprocessing.ipynb` (R)

- Downloads GEO supplementary files and series metadata via `GEOquery` (`getGEOSuppFiles`, `getGEO`)
- Builds `mapping_table.csv` linking GEO accessions, sample titles, and source names
- Loads each sample with `Seurat::Read10X`, creates a `Seurat` object (gene expression) with an additional `HTO` assay (antibody capture)
- Normalizes HTO counts with CLR (Centered Log-Ratio) and demultiplexes hashtags with `HTODemux`
- Visualizes HTO signal (ridge plots, heatmaps, violin plots) and retains **singlets** only
- Loops the full pipeline across all sample folders and merges all per-sample singlet objects into a single Seurat object (`merged_singlets.rds`)

### 3. QC, patient/timepoint mapping, integration & clustering — `02_qc_clustering.R.ipynb` (R)

- Computes mitochondrial content (`percent.mt`) and inspects QC metrics (`nFeature_RNA`, `nCount_RNA`, `percent.mt`) by sample
- Builds a manual **hashtag → patient/timepoint mapping table** and joins it to cell-level metadata to annotate each cell with `patient_id` and `timepoint`
- Flags the PT15 replicate run (`sequencing_run` column) and performs a QC comparison between the original and replicate runs
- Applies QC filtering (`200 < nFeature_RNA < 6000`, `nCount_RNA < 30000`, `percent.mt < 15`)
- Runs `SCTransform` per patient (regressing out `percent.mt`), then merges and selects integration features
- Integrates samples with **Harmony** (`IntegrateLayers`), computes PCA/UMAP, and clusters cells (`FindNeighbors` + `FindClusters`)
- Saves diagnostic UMAPs (by patient, sequencing run, cluster, timepoint) and the final integrated/clustered object (`merged_integrated_clustered.rds`)

### 4. Marker scoring & integration diagnostics — `03_integration_clustering.ipynb` (R)

- Loads the integrated/clustered Seurat object
- Scores clusters against a curated panel of PBMC lineage markers (T cell, CD4/CD8, Treg, exhaustion markers, NK, monocyte, B cell, plasma cell, dendritic cell, platelet, red blood cell) via `DotPlot`
- Identifies cluster marker genes with `FindAllMarkers` (SCT assay) and exports full and top-5-per-cluster marker tables
- Saves the prepared object for downstream annotation

### 5. Cell-type annotation & exhaustion markers — `04_celltype_annotation.ipynb` (R)

- Loads the integrated/clustered Seurat object
- Performs manual PBMC lineage annotation by mapping clusters to cell types based on marker panels
- Subsets T cells and visualizes exhaustion markers (`PDCD1`, `CTLA4`, `LAG3`, `HAVCR2`, `TOX`, `TIGIT`) via `FeaturePlot` and `DotPlot`
- Saves the annotated object (`merged_annotated_final.rds`) and the T-cell subset (`tcells_subset.rds`)

### 6. Exhaustion quantification — `05_exhaustion_quantification.ipynb` (R, work in progress)

- Quantifies T-cell exhaustion scores across patients and timepoints
- Applies non-parametric statistics (Friedman test) for longitudinal comparison

## Requirements

**R** (primary analysis environment, R ≥ 4.5). Key packages used across the notebooks:

`Seurat` (v5), `harmony`, `GEOquery`, `glmGamPoi`, `dplyr`, `ggplot2`, `future`, `SingleR`, `celldex`, `RColorBrewer`, `patchwork`

**Python** (used only for raw-data organization):

`pandas`, plus standard-library modules `os`, `re`, `shutil`, `gzip`

**Data/results versioning:**

`dvc` (with a DagsHub remote — see [Data & Results Versioning](#data--results-versioning))

## Known limitations & reproducibility

This repository is a work in progress. The following are known gaps being addressed:

- **Hardcoded paths:** notebooks and scripts currently contain absolute `/Users/...` paths. Set the `SCRNA_PROJ_ROOT` environment variable (see [Usage](#usage)); full path externalization is in progress.
- **Random seeds:** `set.seed()` is not yet set in all notebooks, so clustering/integration results may vary between runs.
- **Notebook 05** is incomplete (empty placeholder) and not yet executable.
- **No `dvc.yaml` pipeline yet:** a DVC DAG is planned but not yet defined; data/results are tracked as DVC pointers.
- **Empty placeholders:** `docs/report.md`, `src/README.md`, and `results/README.md` are reserved placeholders.
- **No LICENSE / tests / CI yet:** licensing, unit tests, and continuous integration are on the roadmap.

## Citation

If you use this dataset, please cite the associated publication:

- PubMed: [40610460](https://www.ncbi.nlm.nih.gov/pubmed/40610460)

Please also cite the GEO series:

> Sun X, Axelrod ML, et al. Single cell RNA sequencing for longitudinal human peripheral blood from HR+ breast cancer patients treated with immunotherapy. GEO accession GSE300475.

## Contact

For questions about the original dataset, contact the submitting lab:

- **Justin Balko** — Vanderbilt University — justin.balko@vumc.org

For questions about this analysis repository, please open an issue on GitHub.
