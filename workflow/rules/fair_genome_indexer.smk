if config.get("fair_genome_indexer_snakefile"):

    module fair_genome_indexer:
        snakefile:
            config.get("fair_genome_indexer_snakefile")
        config:
            config

else:

    module fair_genome_indexer:
        snakefile:
            github(
                "tdayris/fair_genome_indexer", path="workflow/Snakefile", tag="3.8.0"
            )
        config:
            config


use rule * from fair_genome_indexer
