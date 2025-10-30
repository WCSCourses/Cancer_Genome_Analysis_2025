# Group Mini Projects: Integrated Cancer Genomics Analysis

## Overview

In these group mini projects, you will apply the analytical skills acquired throughout this course to investigate real cancer genomics data. Working with TCGA (The Cancer Genome Atlas) data from distinct cancer types, your group will develop and test scientific hypotheses integrating concepts from driver gene identification, mutational signatures, and statistical association testing.

## Learning Objectives

By completing this mini project, you will:

1. **Apply course methods** to real genomic data in an integrated analysis
2. **Develop scientific hypotheses** based on exploratory analysis of mutation and clinical data
3. **Test hypotheses** using appropriate statistical and bioinformatic approaches
4. **Interpret results** in a biological and clinical context
5. **Communicate findings** in a scientific format with reproducible analysis

## Project Overview

Each group will work with one of four TCGA cancer types:

| Cancer Type | Directory | Abbreviation |
|-------------|-----------|--------------|
| Liver Hepatocellular Carcinoma | `lihc_tcga_pan_can_atlas_2018/` | LIHC |
| Lung Squamous Cell Carcinoma | `lusc_tcga_pan_can_atlas_2018/` | LUSC |
| Pancreatic Adenocarcinoma | `paad_tcga_pan_can_atlas_2018/` | PAAD |
| Skin Cutaneous Melanoma | `skcm_tcga_pan_can_atlas_2018/` | SKCM |

## Data Available

Each cancer type directory contains:

### 1. Mutation Data (`data_mutations.txt`)
MAF-formatted file with somatic mutation information including:
- Genomic coordinates and gene annotations
- Variant classifications and functional predictions
- Sample identifiers for linking to clinical data

### 2. Clinical Data (`*_clinical_data.tsv`)
Clinical and demographic information including:
- Patient characteristics (age, sex, stage, grade)
- Computed genomic metrics (mutation count, TMB, aneuploidy score)
- Survival outcomes (Overall Survival, Disease-Free Survival, Progression-Free Survival)
- Additional clinical variables

The example data for all groups is available online, hosted by Phileal. You can download the file using the following command:

```bash
wget https://public.phileal.com/tutorials/cancer_genome_analysis_2025/Data_projects_WCS_LATAM_2023.zip
```

Once downloaded, find your assigned group's cancer type and proceed with your analysis.

## Getting Started

### Suggested Exploration Areas

Consider investigating research questions such as:

- **Which genes are recurrently mutated** in your cancer type, and do mutations in specific genes associate with clinical outcomes?
- **What mutational processes** shaped these tumors, and do different samples exhibit distinct mutational patterns?
- **How do mutation patterns correlate** with patient demographics, disease stage, or survival?
- **Are there patient subgroups** identifiable by their genomic or mutational signature profiles?
- **What combinations of mutations** appear in aggressive vs. indolent tumors?

### Analytical Considerations

Think about:
- How to characterize the mutation burden and patterns across samples
- What clinical variables might be most informative for hypothesis testing
- How to identify and validate associations you discover
- Whether observed patterns are statistically robust and biologically plausible
- What limitations or confounders might affect your conclusions

### Tools and Methods from Course

Recall the analytical approaches covered in:
- **Module 2-3**: Variant analysis, visualization, and statistical testing
- **Module 3**: Driver gene identification and mutational burden analysis
- **Module 4**: Mutational signature analysis and interpretation

Choose tools and methods appropriate for your research questions.

## Presenting Your Findings

Prepare a short (5-10 minutes) presentation of your findings.
