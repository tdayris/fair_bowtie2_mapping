# 4.3.0

## Features:

* Mosdepth included in QC

# 4.2.0

## Features:

* Samtools idxstats included

## Fixes:

* Deeptools alignment sieve producing non position-sorted bams 

# 4.1.3

## Features:

* bt2 index explicit

# 4.1.2

## Features:

* Reduce fastp threads

# 4.1.1

## Features:

* Allow local snakemake wrappers
* Allow local modules

# 4.1.0

## Features:

* Fastp output compressed

# 4.0.0

## Features:

* Deeptools alignment sieve (on `params/make_sieve` set to `true`)
* snakemake_wrappers up to 3.13.7

# 3.5.2

## Features:

* fair_genome_indexer up to 3.8.1
* Memory / Time reservations

# 3.5.1

## Features:

* fair_genome_indexer up to 3.8.0
* fair_fastqc_multiqc up to 2.3.5

# 3.5.0

## Features:

* fair_genome_indexer up to 3.7.0
* fair_fastqc_multiqc up to 2.3.4

## Documentation:

* Fix links

# 3.4.1

## Features:

* Better reservations

# 3.4.0

## Features:

* Fix benchmark IO for rseqc bamstat
* Use of wrappers for ngsderive
* tmp/logs/benchmark hold rulenames
* fair_fastqc_multiqc and fair_genome_indexer updates
* Pipeline can now localy import sub-modules


# 3.3.3

## Features:

* fair_fastqc_multiqc update to 2.2.7

# 3.3.2

## Features:

* Use human readable functions to replace raw lookups
* Documentation update

# 3.3.1

## Features:

* fair_genome_indexer update to 3.4.0
* fair_fastqc_multiqc update to 2.2.3
* snakemake-wrappers update to 3.7.0

# 3.3.0

## Featues: 

* All keys in configuration are now optional
* Snakemake wrappers up to version 3.5.0

## Documentation:

* Pipeline display update

# 3.2.0

# Features:

* Snakemake-wrappers update to 3.5.0
* MultiQC configurations
* RSeQC: bamstat, read_distribution, read_gc, infer_experiment, and read_duplication
* Goleft: indexcov
* ngsderive: strandedness, instrument, readlen, and encoding


# 3.1.1

## Fix:

* MultiQC mapping report now holds only the samples belonging to a common `{species}.{build}.{release}`
* Readme points to the right snakemake version

# 3.1.0

Update to Snakemake v8+

## Features:

* Usage updated
* DAG as ascii art
* tempfiles, logs and benchmarks paths reorganized: 
    * `tmp/fair_fastqc_multiqc/{rule_name}/{wildcards}.{extension}`
    * `log/fair_fastqc_multiqc/{rule_name}/{wildcards}.log`
    * `benchmark/fair_fastqc_multiqc/{rule_name}/{wildcards}.tsv`
* use of `lookup` and `collect` rather than hand made functions
* All keys in configuration must be present in configuration file

# 3.0.1

## Features:

* Genome schema validation update

## Fixes:

* Report testing removed as long as TBD issue is opened
* Do not download un-necessary reference files anymore

# 3.0.0

Breaking change: Non canonical chromosomes removed by default

## Features:

* fair_genome_indexer update to [3.0.0](https://github.com/tdayris/fair_genome_indexer/releases/tag/3.0.0)
* snakemake-wrappers updated to [v3.3.3](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/index.html)
* snakemake github action updated to [1.25.1](https://github.com/snakemake/snakemake-github-action/releases/tag/v1.25.1)
* Documentation update

## Fix:

* Wrong parameters passed to MultiQC when importing pipeline with no parameters

## Documentation:

* Pipeline description updated
* Usage generalized
* Gustave Roussy users have a dedicated usage section


# 2.3.0

## Features:

* Include fastqc/multiqc as a module
* Rename meta-wrapper rules to indetify them
* Define generic genomes
* Add configuration keys
* Add MultiQC report label

## Fix:

* Mat&Met typos
* Multiqc report missing tools
* test Makefile update for Snakemake 8+

# 2.2.8

## Features:

* Reduce time and memory reservations
* Snakemake-wrappers up to 3.3.3
* Snakemake 8+ compatibility

## Fix:

* csv.Sniffer getting too much data

# 2.2.7

## Features:

* Snakemake-wrappers up to 3.2.0

## Fix:

* linter-error

# 2.2.6

## Features:

* FastQC in addition to pre-trimming fastp

## Fix:

* Fix wrong reference linking in picard

# 2.2.5

## Features:

* Update fair_genome_indexer to version 2.3.0
* Samtools stats added to report

## Fix:

* Report main page too large
* Fix error removing Bowtie2 QC from MultiQC report
* Fix error making fastp reports non reachable in MultiQC

# 2.2.4

## Features:

* reduce size of configuration provided to fair_genome_indexer
* Include Bowtie2 metrics to MultiQC report

## Fix:

* Fix picard metrics extra parameter missing
* `genomes` not found while importing workflow as module
* Bowtie2-build now takes fasta, if fasta is provided and no bowtie2 index is provided
* Single-end error when mixed single-ended/pair-ended samples were mixed in design

# 2.2.3

## Features:

* Conditionally load fair_genome_indexer
* Let user provide annotation files

# 2.2.2

## Features:

* Documentation update
* Rename loaded rules from fair_genome_indexer rule to indetify them

# 2.2.1

## Features:

* Protect bam output
* Speed-up with wildcards constraints
* Speed-up Picard with BAI inclusion in input files
* Changed output directory order

## Fix:

* snakemake wrappers version fixed


# 2.2.0

## Features

* fastp reports included in snakemake report
* fair_genome_indexer pipeline version updated to 2.2.0
* snakemake wrappers version updated to v3.0.0
* Description of expected results with directoy architecture

# 2.1.1

## Features

* fair_genome_indexer pipeline update

# 2.1.0

## Features:

* Add report texts

## Documentation

* Fastp citation

# 2.0.0

## Features:

* Fastp-trimming (reports included in MultiQC)

# 1.0.0

## Features:

* Download-index reference genome if missing
* Build bowtie2-index if missing
* Align reads
* Perform post-mapping QC with picard and aggregate reports with MultiQC
