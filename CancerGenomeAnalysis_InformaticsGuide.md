# <img src="https://coursesandconferences.wellcomeconnectingscience.org/wp-content/themes/wcc_courses_and_conferences/dist/assets/svg/logo.svg" width="300" height="50"> 
# Cancer Genome Analysis Informatics Guide

## **Software used during the course**

| Software | Version (if not latest) | Module | Notes |
|-------------|--------------|----------|-------------|
| [R](https://www.r-project.org/) | 4.3.3 | Statistical Computing | Required for many bioinformatics tools |
| [RStudio](https://posit.co/download/rstudio-desktop/) | 4.3.3 | IDE for R | Recommended for R users |
| [Samtools](http://www.htslib.org/) | 1.9 | NGS Analysis | Used for BAM file processing |
| [BCFtools](http://www.htslib.org/) | 1.9 | Variant Calling | Works with VCF/BCF files |
| [GATK](https://gatk.broadinstitute.org/hc/en-us) | 4.6.1.0 | Variant Calling | Genome Analysis Toolkit |
| [dndscv](https://github.com/im3sanger/dndscv) | Latest | Mutation Analysis | dN/dS ratio calculations |
| [CIBERSORT](https://cibersortx.stanford.edu/) | Latest | Deconvolution | Immune cell type estimation |
| [maftools](https://bioconductor.org/packages/release/bioc/html/maftools.html) | Latest | Mutation Analysis | Visualization of somatic mutations |
| [VEP (Variant Effect Predictor)](https://www.ensembl.org/info/docs/tools/vep/index.html) | Latest | Variant Annotation | Provided by Ensembl |
| [SnpEff](https://pcingola.github.io/SnpEff/) | 5.2 | Variant Annotation | Predicts variant effects |
| [BWA](http://bio-bwa.sourceforge.net/) | 7.18 | NGS Alignment | Aligns sequences to reference genome |
| [SigProfilerExtractor](https://github.com/AlexandrovLab/SigProfilerExtractor) | Latest | Mutation Signatures | Extracts mutational signatures |
| [SigProfilerMatrixGenerator](https://github.com/AlexandrovLab/SigProfilerMatrixGenerator) | Latest | Mutation Signatures | Generates matrices for signature analysis |
| [SigProfilerAssignment](https://github.com/AlexandrovLab/SigProfilerAssignment) | Latest | Mutation Signatures | Assigns signatures to samples|


## **To run the SigProfiler software, you will have to activate an environment (mutational_signatures), What is an Environment?**
An environment is an isolated space where specific software, dependencies, and libraries are installed. It ensures that all required tools run in a controlled and reproducible way, avoiding conflicts with other system applications.

## For this course, we use the **mutational_signatures** and **samtools_env** environment.

⚠️ **Note:** Software will **not work outside the environment.**

## **How to Activate the mutational_signatures Environment**
Before using any software, activate the environment with:

```bash
conda activate mutational_signatures
```

## **How to Activate the samtools_env Environment**
Before using any software, activate the environment with:

```bash
conda activate samtools_env
```


## **To deactivate the environment when you're done:**

```bash
conda deactivate
```

Link to [mutational_signatures.yml](https://github.com/WCSCourses/Cancer_Genome_Analysis_2025/blob/main/course_data_2025/env_mutational_signatures.yml)
Installation Guide here [SigProfilerHelper](https://github.com/edawson/SigProfilerHelper)

## Informatics Set-Up
For installation and setup, please refer to the following guides:

- **[Oracle VM VirtualBox Installation Guide](https://github.com/WCSCourses/index/blob/main/VM_Guide.md)** – Detailed instructions for installing and configuring VirtualBox on different operating systems. *(Note: Separate installations are needed for Intel-based and ARM-based Macs, and the VDI files will differ.)*

The Host Operating System Requirements are: <br />
- RAM requirement: 8GB (preferably 12GB) <br />
- Processor requirement: 4 processors (preferably 8) <br />
- Hard disk space: 200GB <br />
- Admin rights to the computer <br />

## Citing and Reusing Course Material

The course data are free to reuse and adapt with appropriate attribution. All course data in these repositories are licensed under the <a rel="license" href="https://creativecommons.org/licenses/by-nc-sa/4.0/">Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)</a>. <a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br /> 

Each course landing page is assigned a DOI via Zenodo, providing a stable and citable reference. These DOIs can be found on the respective course landing pages and can be included in CVs or research publications, offering a professional record of the course contributions.

## Interested in attending a course?

Take a look at what courses are coming up at [Wellcome Connecting Science Courses & Conference Website](https://coursesandconferences.wellcomeconnectingscience.org/our-events/).

---

[Wellcome Connecting Science GitHub Home Page](https://github.com/WCSCourses) 

For more information or queries, feel free to contact us via the [Wellcome Connecting Science website](https://coursesandconferences.wellcomeconnectingscience.org).<br /> 
Please find us on socials [Wellcome Connecting Science Linktr](https://linktr.ee/eventswcs)

---
