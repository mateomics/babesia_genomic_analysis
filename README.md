# Comparative and Functional Genomics of *Babesia divergens*: Invasion and Egress Processes

This repository contains the bioinformatics pipeline and analysis for the study of gene expression patterns in *Babesia divergens*, focusing on the molecular processes of **invasion** and **egress**. This project reproduces and extends the analysis from González et al. (2019) through RNA-seq differential expression analysis.

## Authors

- Héctor López Ordaz <hectorjl@lcg.unam.mx>
- Ismael Maximiliano De Los Santos Huesca <ismadlsh@lcg.unam.mx>
- Mateo Jiménez Sotelo <majiso@lcg.unam.mx>

## Project Overview

Babesiosis is an emerging zoonotic disease caused by the protozoan parasite *Babesia divergens*, whose incidence has increased significantly over the last 30 years. This project analyzes gene expression patterns associated with invasion and egress processes using a comprehensive bioinformatics pipeline.

## Data Sources

- **RNA-seq Data**: [Fuente](link)
- **Reference Genome**: [Fuente](link)
- **Original Study**: González LM, et al. (2019) PLoS Negl Trop Dis 13(8): e0007680

## Project Structure

babesia_genomic_analysis/  
├── bin/ # Executable scripts  
├── data/ # Raw and processed data
├── docs/ # Documentation  
│ ├──Comparative_and_functional_genomics_of_the_protozoan_Babesia_divergens_highlighting_the_invasion_and_egress_processes.pdf  
│ └── Project_Description_Babesia.md
├── results/ # Analysis results  
├── src/ # Source code  
├── lib/ # Libraries and dependencies  
├── test/ # Test files  
├── LICENSE  
└── README.md


## Analysis Pipeline

1.  **Data Acquisition**: Download FASTQ files and reference genome from NCBI
2.  **Quality Control**: FastQC and MultiQC for read quality assessment
3.  **Alignment**: HISAT2 alignment to *B. divergens* genome
4.  **Quantification**: featureCounts for read counting
5.  **Differential Expression**: DESeq2/edgeR analysis
6.  **Functional Annotation**: BLAST analysis of DEGs
7.  **Validation**: Comparison with original study results

## Installation and Usage

[*Installation instructions and usage examples will be added here as the project develops*]

## Expected Results

- Gene expression matrix of *B. divergens*
- Identification of differentially expressed genes between free and intracellular merozoites
- Functional annotation of invasion and egress-related genes
- Validation of original study findings
- Reproducible analysis pipeline

## Terms of Use

This project is available under the MIT license. See the LICENSE file for more details.

## Citation

If you use this work in your research, please cite:  
López Ordaz H, De Los Santos Huesca IM, Jiménez Sotelo M. (2025). Comparative and functional genomics of Babesia divergens highlighting the invasion and egress processes. [GitHub Repository]. https://github.com/mateomics/babesia_genomic_analysis

## Contact

For questions or issues regarding this project, please open an issue in this repository or contact:
- Héctor López Ordaz: hectorjl@lcg.unam.mx
- Ismael Maximiliano De Los Santos Huesca: ismadlsh@lcg.unam.mx
- Mateo Jiménez Sotelo: majiso@lcg.unam.mx

## References

González LM, Estrada K, Grande R, et al. (2019) Comparative and functional genomics of the protozoan parasite Babesia divergens highlighting the invasion and egress processes. PLoS Negl Trop Dis 13(8): e0007680. https://doi.org/10.1371/journal.pntd.0007680
