# Explore Visualization - RNA Sequencing and ChIP-Sequencing

An interactive bioinformatics data visualization project that involves analyzing RNA-seq and  ChIP-Seq data using existing packages like ChIpSeeker, DESeq2, and GProfiler. The goal is for students to extract insights from the data, draw conclusions about proteins, and defend their findings using visualizations. This can be performed on a specific protein where we have control data, protein-knockout data.


**KEY IDEAS:**
This project is VERY open-ended by design. We want to give students some creative ability here to make a conclusion using figures and results. A student will easily be able to look up the function of a given protein BUT can they use their data to reaffirm what they are seeing online? That is the challenge here!


## Data Preprocessing

Explain the importance of data preprocessing in bioinformatics analysis. Showcase how to perform basic preprocessing steps using existing packages (e.g., quality control, read alignment, normalization). This can be done using a series of bioinformatics tools that we can embed into a docker container for the students. 

**Important Note:**
This process is VERY computationally expensive. We will assert that students are able to run this pipeline on a very simple dataset and then will provide them with the full counts tables that are inputs for DESeq2! 

## Differential Expression Analysis (RNA-seq)

Introduce DESeq2 and its role in identifying differentially expressed genes.
Provide sample code demonstrating how to use DESeq2 for differential expression analysis.
Explain the concept of fold change and adjusted p-values.
Showcase visualization techniques for displaying differentially expressed genes (e.g., volcano plots, heatmaps).

## Functional Enrichment Analysis (RNA-seq)

Introduce gprofiler2 and its role in functional enrichment analysis. Using the genes identified from DESeq2, the students can perform functional enrichment analysis to see which 
Demonstrate how to perform functional enrichment analysis using gprofiler2.
Show how to visualize enriched Gene Ontology terms and pathways.

## ChIP-Seq Analysis

Introduce ChIP-Seq data and its application in studying protein-DNA interactions.
Provide an overview of peak calling and its significance.
Showcase how to use existing tools for ChIP-Seq analysis (e.g., MACS2).
Display visualization techniques for ChIP-Seq data, such as peak density plots and motif analysis.