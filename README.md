# RNAseq-Analysis
Comprehensive RNA-Seq Analysis and Differential Gene Expression Profiling of Human Cardiopulmonary Samples (GSE236251)

**Overview:**

This project presents a fully automated and reproducible workflow for the analysis of RNA-Seq gene expression data using HTSeq-count output files, focusing specifically on the publicly available human dataset GSE236251 from NCBI’s GEO repository.

The workflow performs key steps in RNA-Seq analysis including data import, preprocessing, normalization, gene annotation, and differential expression analysis using the DESeq2 framework. Additionally, the pipeline generates high-quality visualizations such as volcano plots, heatmaps, and PCA plots to aid in biological interpretation of results.

**The pipeline is designed to identify differentially expressed genes (DEGs) between three main clinical conditions:**

✔️Healthy controls

✔️Left Heart Disease without Pulmonary Hypertension (LHD w/o PH)

✔️Pulmonary Hypertension with Left Heart Disease (PH-LHD)

**Main Features:**

✅ Batch processing of multiple sample-level count files

✅ Sample metadata handling from custom CSV files

✅ Data normalization and transformation using DESeq2

✅ Gene annotation using biomaRt and org.Hs.eg.db

✅ Filtering and exporting of significantly up/downregulated genes

✅ Publication-ready DE plots (MA, Volcano, PCA, Heatmaps)

**Goals of the Project:**

✔️Integrate and clean high-throughput RNA-Seq datasets

✔️Identify differentially expressed genes (DEGs) associated with disease phenotypes

✔️Provide gene-level biological annotations and enrichment potential

✔️Offer visualizations for downstream interpretation or machine learning tasks

**Who Can Use This?**
This workflow is particularly useful for bioinformaticians, molecular biologists, and clinical researchers working in:

✔️Transcriptomics

✔️Cardio-pulmonary disease genomics

✔️Biomarker discovery

✔️Differential expression analysis

**Technologies Used:**
R, DESeq2, data.table, biomaRt, AnnotationDbi, ComplexHeatmap, ggplot2, pheatmap
