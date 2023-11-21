module fair_genome_indexer:
    snakefile: github("tdayris/fair_genome_indexer", path="workflow/Snakefile", tag="2.0.1")

use rule * from fair_genome_indexer