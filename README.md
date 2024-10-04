[![Snakemake](https://img.shields.io/badge/snakemake-≥7.29.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/tdayris/fair_genome_indexer/workflows/Tests/badge.svg?branch=main)](https://github.com/tdayris/fair_genome_indexer/actions?query=branch%3Amain+workflow%3ATests)

Snakemake workflow used to align ungapped reads to the genome with Bowtie2.

## Usage

The usage of this workflow is described in the [Snakemake workflow catalog](https://snakemake.github.io/snakemake-workflow-catalog?usage=tdayris/fair_bowtie2_mapping) it is also available [locally](https://github.com/tdayris/fair_bowtie2_mapping/blob/main/workflow/report/usage.rst) on a single page.
 
## Results

A complete description of the results can be found here in [workflow reports](https://github.com/tdayris/fair_bowtie2_mapping/blob/main/workflow/report/results.rst).

## Material and Methods

The tools used in this pipeline are described [here](https://github.com/tdayris/fair_bowtie2_mapping/blob/main/workflow/report/material_methods.rst) textually. Web-links are available below:


### Index and genome sequences with [`fair_genome_indexer`](https://github.com/tdayris/fair_genome_indexer/)

See [`fair_genome_indexer`](https://github.com/tdayris/fair_genome_indexer/) documentation about DNA sequence preparation

### Raw-sequences QC with [`fair_fastqc_multiqc`](https://github.com/tdayris/fair_fastqc_multiqc/)

See  [`fair_fastqc_multiqc`](https://github.com/tdayris/fair_fastqc_multiqc/) documentation about ranw sequences quality controls

### Bowtie2 Mapping

| Step          | Meta-Wrapper                                                                                                              | Wrapper                                                                                           |
| ------------- | ------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------- |
| Bowtie2-build | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v4.6.0/meta-wrappers/bowtie2_sambamba.html) | [bowtie2-build](https://snakemake-wrappers.readthedocs.io/en/v4.6.0/wrappers/bowtie2/build.html) |
| Fastp         |                                                                                                                           | [fastp](https://snakemake-wrappers.readthedocs.io/en/v4.6.0/wrappers/fastp.html)                  |
| Bowtie2-align | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v4.6.0/meta-wrappers/bowtie2_sambamba.html) | [bowtie2-align](https://snakemake-wrappers.readthedocs.io/en/v4.6.0/wrappers/bowtie2/align.html) |
| Sambamba sort | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v4.6.0/meta-wrappers/bowtie2_sambamba.html) | [sambamba-sort](https://snakemake-wrappers.readthedocs.io/en/v4.6.0/wrappers/sambamba/sort.html) |

```
┌───────────────────────────┐   ┌─────────────────────────┐
│Genome indexation (Bowtie2)│   │Sequence cleaning (fastp)│
└──────────────┬────────────┘   └────────┬────────────────┘
               │                         │                 
               │                         │                 
┌──────────────▼───────────────┐         │                 
│Short read alignment (Bowtie2)◄─────────┘                 
└──────────────┬───────────────┘                           
               │                                           
               │                                           
┌──────────────▼───────────────────────────┐               
│Sort and compress aligned reads (sambamba)│               
└──────────────────────────────────────────┘               
```


### Filtering

| Step             | Meta-Wrapper                                                                                                             | Wrapper                                                                                                       |
| ---------------- | ------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------- |
| Sambamba-view    | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v4.6.0/meta-wrappers/bowtie2_sambamba.html) | [sambamba-view](https://snakemake-wrappers.readthedocs.io/en/v4.6.0/wrappers/sambamba/view.html)              |
| Sambamba-markdup | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v4.6.0/meta-wrappers/bowtie2_sambamba.html) | [sambamba-markdup](https://snakemake-wrappers.readthedocs.io/en/v4.6.0/wrappers/sambamba/markdup.html)        |
| Sambamba-index   | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v4.6.0/meta-wrappers/bowtie2_sambamba.html) | [sambamba-index](https://snakemake-wrappers.readthedocs.io/en/v4.6.0/wrappers/sambamba/index.html)            |
| Deeptools        |                                                                                                                          | [deeptools-alignment-sieve](https://snakemake-wrappers.readthedocs.io/en/v4.6.0/bio/deeptools/alignmentsieve) |

```
 ┌─────────────────────────┐                                                
 │Aligned reads (see above)│                                                
 └───────────┬─────────────┘                                                
             │                                                              
             │                                                              
             │                                                              
             │                                                              
 ┌───────────▼─────────────────┐                                            
 │Filter low quality (sambamba)│                                            
 └───────────┬─────────────────┘                                            
             │                                                              
             │                                                              
             │                                                              
 ┌───────────▼──────────────┐           ┌──────────────────────────────────┐
 │Mark duplicates (sambamba)├───────────►Index aligned sequences (sambamba)│
 └───────────┬──────────────┘           └─────────┬────────────────────────┘
             │                                    │                         
             ├────────────────────────────────────┘                         
             │                                                              
┌────────────▼─────────────────┐        ┌──────────────────────────────────┐
│ Sieve alignments (deeptools) ├────────►Index aligned sequences (sambamba)│
│     OPTIONAL                 │        │     OPTIONAL                     │
└──────────────────────────────┘        └──────────────────────────────────┘
```

### QC

| Step     | Wrapper                                                                                                                           |
| -------- | --------------------------------------------------------------------------------------------------------------------------------- |
| Picard   | [picard-collectmultiplemetrics](https://snakemake-wrappers.readthedocs.io/en/v4.6.0/wrappers/picard/collectmultiplemetrics.html) |
| Samtools | [samtools-stats](https://snakemake-wrappers.readthedocs.io/en/v4.6.0/wrappers/samtools/stats.html)                               |
| MultiQC  | [multiqc-wrapper](https://snakemake-wrappers.readthedocs.io/en/v4.6.0/wrappers/multiqc.html)                                     |

```
┌──────────────────────┐        ┌─────────────────────┐                ┌─────────────────────────┐
│ Cleaned reads (fastp)│    ┌───┤Duplicates (sambamba)◄────────────────┤Aligned reads (see above)│
└─────────────────────┬┘    │   └─────────────────────┘                └────┬────────────────────┘
                      │     │                                               │                     
                      ├─────┘                                               │                     
                      │         ┌──────────────────────────┐                │                     
                      ├─────────┤Alignment metrics (picard)◄────────────────┤                     
                      │         └──────────────────────────┘                │                     
                      │                                                     │                     
                      │                                                     │                     
                      │         ┌────────────────────────────┐              │                     
                      ├─────────┤Alignment metrics (samtools)◄──────────────┤                     
                      │         └────────────────────────────┘              │                     
                      │                                                     │                     
┌────────────────┐    │                                                     │                     
│ Quality report │    │         ┌─────────────────────────┐                 │                     
│   (multiqc)    ◄────┼─────────┤Alignment metrics (rseqc)◄─────────────────┤                     
└────────────────┘    │         └─────────────────────────┘                 │                     
                      │                                                     │                     
                      │                                                     │                     
                      │         ┌───────────────────────────┐               │                     
                      ├─────────┤Library metrics (ngsderive)◄───────────────┤                     
                      │         └───────────────────────────┘               │                     
                      │                                                     │                     
                      │                                                     │                     
                      │         ┌─────────────────────────┐                 │                     
                      └─────────┤Coverage metrics (goleft)◄─────────────────┘                     
                                └─────────────────────────┘                                       
```
