rule picard_create_multiple_metrics:
    input:
        bam="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
        bai="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai",
        ref="reference/{species}.{build}.{release}.{datatype}.fasta",
    output:
        temp(
            multiext(
                "tmp/picard/{species}.{build}.{release}.{datatype}/stats/{sample}",
                ".alignment_summary_metrics",
                ".insert_size_metrics",
                ".insert_size_histogram.pdf",
                ".base_distribution_by_cycle_metrics",
                ".base_distribution_by_cycle.pdf",
                ".gc_bias.detail_metrics",
                ".gc_bias.summary_metrics",
                ".gc_bias.pdf",
            )
        ),
    log:
        "logs/picard/collectmultiplemetrics/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/picard/collectmultiplemetrics/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=config.get("params", {}).get("picard", {}).get("metrics", ""),
    wrapper:
        f"{snakemake_wrappers_version}/bio/picard/collectmultiplemetrics"
