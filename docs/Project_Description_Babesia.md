# Comparative and functional genomics of the protozoan parasite Babesia divergens highlighting the invasion and egress processes

Fecha: 26/08/2025

**Authors**:

- Héctor López Ordaz	<hectorjl@lcg.unam.mx>
- Ismael Maximiliano De Los Santos Huesca	<ismadlsh@lcg.unam.mx>
- Mateo Jiménez Sotelo	<majiso@lcg.unam.mx>

## Statement of the problem

Babesiosis is an emerging zoonotic disease caused by the protozoan parasite *Babesia divergens*, whose incidence has increased significantly over the last 30 years. The parasite's molecular processes of **invasion** and **egress** are not fully characterized, which limits the development of effective diagnostic and treatment methods. This project seeks to analyze the gene expression patterns associated with these processes using a bioinformatics pipeline.

## Calendar

[Pending to define. General activities required for the project]

| Activity | Date   | Responsible  | Submit |
|----------|----------|----------|----------|
| Project description    | August 30  | Mateo, Héctor   | Markdown document |
| README    | August 30  | Mateo, Héctor   | Markdown document |
| Requirements specification    | September 6   | Max   | Markdown document   |
| Analysis and design   | September 13  |  Mateo, Héctor, Max   | Markdown document |
| Construction   | September 20, November  |  Mateo, Héctor, Max   | Scripts |
| Testing   | October 4  |  Mateo, Héctor, Max   | Markdown document |
| Results report  | October 12  |  Mateo, Héctor, Max   | Markdown documents |
| Project presentation   | [pending]  |  Mateo, Héctor, Max   | GitHub repository (release)|




## Methodology

Pipeline:

1.  Download the fastq files

2.  Download the Babesia divergens genome from NCBI

3.  Align the sequences of each replicate of each condition

4.  Generate the count table
    
5.  Perform differential expression analysis
    
6.  For differentially expressed genes, search for their annotation via blast

7.  Analyze the output

## Expected results

- Gene expression matrix of *B. divergens*
- Identification of differentially expressed genes between free and intracellular merozoites
- Functional annotation of genes associated with invasion and egress
- Validation of results reported in the original article
- Reproducible repository of the complete analysis


## Requirements Specification

- Automated data download from NCBI SRA
- Preprocessing and quality control pipeline
- Efficient alignment of reads to the genome
- Generation of expression count matrix
- Statistical analysis of differential expression
- Functional annotation using BLAST
- Generation of reports and visualizations


### Non-functional requirements
- Complete pipeline documentation
- Environment management with Conda
- Version control with Git
- Efficient management of computational resources

## Analysis and design

### Pipeline's architecture 

## Directory structure
babesia_genomic_analysis
├── bin
│   └── REMOVEME.txt
├── data
│   └── REMOVEME.txt
├── docs
│   ├── Comparative and functional genomics of the protozoan Babesia divergens highlighting the invasion and egress processes.pdf
│   └── Descripcion_proyecto_Babesia.md.md
├── lib
│   └── REMOVEME.txt
├── LICENSE
├── raramuri
├── README.md
├── results
│   └── REMOVEME.txt
├── src
│   └── REMOVEME.txt
└── test
    └── REMOVEME.txt

### Usage case: Differential expression analysis

**Actor**: Bioinformático  
**Description**: Ejecución del pipeline completo de análisis de RNA-seq

**Principal flux**:

[pending]
    

**Alternative fluxes**:

[pending]
    

### Computational resources needed

-   Conexión a internet estable
- [pending]
