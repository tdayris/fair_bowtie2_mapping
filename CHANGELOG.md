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