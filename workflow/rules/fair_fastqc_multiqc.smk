if config.get("fair_fastqc_multiqc_snakefile") is not None:

    module fair_fastqc_multiqc:
        snakefile:
            config.get("fair_fastqc_multiqc_snakefile")
        config:
            {**config, "load_fair_genome_indexer": False}

else:

    module fair_fastqc_multiqc:
        snakefile:
            github(
                "tdayris/fair_fastqc_multiqc",
                path="workflow/Snakefile",
                tag="2.3.3",
            )
        config:
            {**config, "load_fair_genome_indexer": False}


use rule * from fair_fastqc_multiqc
