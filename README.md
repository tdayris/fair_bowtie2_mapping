[![Snakemake](https://img.shields.io/badge/snakemake-≥7.29.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/tdayris/fair_genome_indexer/workflows/Tests/badge.svg?branch=main)](https://github.com/tdayris/fair_genome_indexer/actions?query=branch%3Amain+workflow%3ATests)

Snakemake workflow used to align ungapped reads to the genome with Bowtie2.

## Usage

The usage of this workflow is described in the [Snakemake workflow catalog](https://snakemake.github.io/snakemake-workflow-catalog?usage=tdayris/fair_bowtie2_mapping) it is also available [locally](https://github.com/tdayris/fair_bowtie2_mapping/blob/main/workflow/report/usage.rst) on a single page.
 
## Results

A complete description of the results can be found here in [workflow reports](https://github.com/tdayris/fair_bowtie2_mapping/blob/main/workflow/report/results.rst).

## Material and Methods

The tools used in this pipeline are described [here](https://github.com/tdayris/fair_bowtie2_mapping/blob/main/workflow/report/material_methods.rst) textually. Web-links are available below:

```
┌─────────────────────────────────┐                     ┌───────────────────────────────┐
│  fair_genome_indexer_pipeline   │                     │  fair_fastqc_multiqc_pipeline │
└───────────┬─────────────────────┘                     └───────────────────────────┬───┘
            │                                                                       │
            │                                                                       │
            │                                                                       │
            │                                                                       │
            │                                                                       │
  ┌─────────▼───────────┐           ┌────────────────────────┐                      │
  │    Bowtie2          │           │      Fastp             │                      │
  │ Index DNA sequence  │           │ Trimm and check reads  │                      │
  └──────────────────┬──┘           │ qualittiy              ├─────────────────┐    │
                     │              └────────────────────────┘                 │    │
                     │                                                         │    │
                     │                                                         │    │
                     │                                                         │    │
                     │                                                         │    │
                     │                                                         │    │
             ┌───────▼──────────────┐                                          │    │
             │     Bowtie2          │                                        ┌─▼────▼──────────┐
             │  Align trimmed reads ├────────────────────────────────────────►                 │
             └───────┬──────────────┘                                        │     MultiQC     │
                     │                                                       │  Quality report │
                     │                                                       │                 │
             ┌───────▼──────────────────────────┐                            └──▲──▲───────────┘
             │    Sambamba                      │                               │  │
             │  Quality filters + sorting +     │                               │  │
             │  duplicated reads identification │                               │  │
             └────────────────┬───────────┬─────┘                               │  │
                              │           │                                     │  │
                              │     ┌─────▼─────────────┐                       │  │
                              │     │     Samtools      ├───────────────────────┘  │
                              │     │  Quality controls │                          │
                              │     └───────────────────┘                          │
                              │                                                    │
                              │                                                    │
                              │        ┌──────────────────┐                        │
                              │        │   Picard         ├────────────────────────┘
                              └────────► Quality controls │
                                       └──────────────────┘

```

### Index and genome sequences with [`fair_genome_indexer`](https://github.com/tdayris/fair_genome_indexer/)

| Step                             | Commands                                                                                                         |
| -------------------------------- | ---------------------------------------------------------------------------------------------------------------- |
| Download DNA Fasta from Ensembl  | [ensembl-sequence](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/wrappers/reference/ensembl-sequence.html) |
| Remove non-canonical chromosomes | [pyfaidx](https://github.com/mdshw5/pyfaidx)                                                                     |
| Index DNA sequence               | [samtools](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/wrappers/samtools/faidx.html)                     |
| Creatse sequence Dictionary      | [picard](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/wrappers/picard/createsequencedictionary.html)      |

### Raw-sequences QC with [`fair_fastqc_multiqc`](https://github.com/tdayris/fair_fastqc_multiqc/)

| Step       | Wrapper                                                                                        |
| ---------- | ---------------------------------------------------------------------------------------------- |
| FastQC     | [fastqc-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/wrappers/fastqc.html)     |
| FastScreen | [fastq-screen](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/fastq_screen.html) |
| MultiQC    | [multiqc-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/wrappers/multiqc.html)   |

### Bowtie2 Mapping

| Step          | Meta-Wrapper                                                                                                              | Wrapper                                                                                           |
| ------------- | ------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------- |
| Bowtie2-build | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/meta-wrappers/bowtie2_sambamba.html) | [bowtie2-build](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/wrappers/bowtie2/build.html) |
| Fastp         |                                                                                                                           | [fastp](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/wrappers/fastp.html)                  |
| Bowtie2-align | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/meta-wrappers/bowtie2_sambamba.html) | [bowtie2-align](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/wrappers/bowtie2/align.html) |
| Sambamba sort | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/meta-wrappers/bowtie2_sambamba.html) | [sambamba-sort](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/wrappers/sambamba/sort.html) |

### Filtering

| Step             | Meta-Wrapper                                                                                                   | Wrapper                                                                                                 |
| ---------------- | ------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------- |
| Sambamba-view    | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/meta-wrappers/bowtie2_sambamba.html) | [sambamba-view](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/wrappers/sambamba/view.html)       |
| Sambamba-markdup | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/meta-wrappers/bowtie2_sambamba.html) | [sambamba-markdup](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/wrappers/sambamba/markdup.html) |
| Sambamba-index   | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/meta-wrappers/bowtie2_sambamba.html) | [sambamba-index](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/wrappers/sambamba/index.html)     |


### QC

| Step     | Wrapper                                                                                                                          |
| -------- | -------------------------------------------------------------------------------------------------------------------------------- |
| Picard   | [picard-collectmultiplemetrics](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/wrappers/picard/collectmultiplemetrics.html) |
| Samtools | [samtools-stats](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/wrappers/samtools/stats.html)                               |
| MultiQC  | [multiqc-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/wrappers/multiqc.html)                                     |