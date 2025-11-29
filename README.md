# A curated single-cell transcriptomic dataset for human inflammatory skin diseases

This repository contains the source code for the web portal and the comprehensive data analysis pipeline associated with the paper **"A curated single-cell transcriptomic dataset for human inflammatory skin diseases."**

## 📁 Data Analysis Scripts

The `data_analysis_scripts` directory contains the specific scripts used for processing raw data, performing quality control, and conducting downstream analysis. The workflow and corresponding scripts are detailed below:

### 1. Data Preprocessing & QC
*   **`10xre2countsm.R`**
    *   **Function:** Raw Data Conversion.
    *   **Description:** Used to convert the output from Cell Ranger into a standard counts matrix for downstream processing.
*   **`doublet_rate.py`**
    *   **Function:** Doublet Removal.
    *   **Description:** Takes the counts matrix as input to identify and remove potential doublet cells to ensure data quality.

### 2. Single Dataset Analysis
*   **`seurat_multisample.R`**
    *   **Function:** Integration & Clustering.
    *   **Description:** Primarily used for the integration and clustering analysis within individual datasets.
*   **`seurat_multisample_de.R`**
    *   **Function:** Annotation & DE Analysis.
    *   **Description:** Performs cell type annotation for identified clusters and conducts differential expression gene (DEG) analysis between groups.

### 3. Integrated Atlas Analysis
*   **`combined_seurat_analysis_by_subject.r`**
    *   **Function:** Global Integration (Seurat v5).
    *   **Description:** Implements **Seurat v5** workflows to perform a comprehensive integrated analysis across all datasets, organized by subject.

### 4. Functional Enrichment
*   **`cluster_profile_hsa.R`**
    *   **Function:** Pathway Analysis.
    *   **Description:** Executes Gene Ontology (GO) and KEGG pathway enrichment analyses to interpret biological functions.

