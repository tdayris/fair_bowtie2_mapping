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