[![Snakemake](https://img.shields.io/badge/snakemake-≥7.29.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/tdayris/fair_genome_indexer/workflows/Tests/badge.svg?branch=main)](https://github.com/tdayris/fair_genome_indexer/actions?query=branch%3Amain+workflow%3ATests)

Snakemake workflow used to align reads over a reference genome.

Do not use. Dev. Teaching purpose.

## Usage

The usage of this workflow is described in the [Snakemake workflow catalog](https://snakemake.github.io/snakemake-workflow-catalog?usage=tdayris/fair_bowtie2_mapping)

## Content

### Index and genome sequences


| Step                       | Pipeline                                               | Wrapper                                                                                                                              |
| -------------------------- | --------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| Download DNA Fasta         | [fair_genome_indexer](https://github.com/tdayris/fair_genome_indexer) | [ensembl-sequence](https://snakemake-wrappers.readthedocs.io/en/v2.13.0/wrappers/reference/ensembl-sequence.html)                    |
| Download GTF annotation    | [fair_genome_indexer](https://github.com/tdayris/fair_genome_indexer) | [ensembl-sequence](https://snakemake-wrappers.readthedocs.io/en/v2.13.0/wrappers/reference/ensembl-sequence.html)                    |
| Samtools fasta index       | [fair_genome_indexer](https://github.com/tdayris/fair_genome_indexer) | [ensembl-annotation](https://snakemake-wrappers.readthedocs.io/en/v2.13.0/wrappers/reference/ensembl-annotation.html)                |
| Picard Sequence Dicitonary | [fair_genome_indexer](https://github.com/tdayris/fair_genome_indexer) | [picard-createsequencedictionary](https://snakemake-wrappers.readthedocs.io/en/v2.13.0/wrappers/picard/createsequencedictionary.html) |

### Bowtie2 Mapping

| Step          | Meta-Wrapper                                                                                                   | Wrapper                                                                                           |
| ------------- | ------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------- |
| Bowtie2-build | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v2.13.0/meta-wrappers/bowtie2_sambamba.html) | [bowtie2-build](https://snakemake-wrappers.readthedocs.io/en/v2.13.0/wrappers/bowtie2/build.html) |
| Bowtie2-align | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v2.13.0/meta-wrappers/bowtie2_sambamba.html) | [bowtie2-align](https://snakemake-wrappers.readthedocs.io/en/v2.13.0/wrappers/bowtie2/align.html) |
| Sambamba sort | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v2.13.0/meta-wrappers/bowtie2_sambamba.html) | [sambamba-sort](https://snakemake-wrappers.readthedocs.io/en/v2.13.0/wrappers/sambamba/sort.html) |

### Filtering

| Step             | Meta-Wrapper                                                                                                   | Wrapper                                                                                                 |
| ---------------- | ------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------- |
| Sambamba-view    | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v2.13.0/meta-wrappers/bowtie2_sambamba.html) | [sambamba-view](https://snakemake-wrappers.readthedocs.io/en/v2.13.0/wrappers/sambamba/view.html)       |
| Sambamba-markdup | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v2.13.0/meta-wrappers/bowtie2_sambamba.html) | [sambamba-markdup](https://snakemake-wrappers.readthedocs.io/en/v2.13.0/wrappers/sambamba/markdup.html) |
| Sambamba-index   | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v2.13.0/meta-wrappers/bowtie2_sambamba.html) | [sambamba-index](https://snakemake-wrappers.readthedocs.io/en/v2.13.0/wrappers/sambamba/index.html)     |