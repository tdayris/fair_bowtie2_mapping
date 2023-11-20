rule picard_create_multiple_metrics_raw:
    input:
        "results/bowtie2/{species}.{build}.{release}.{datatype}/{sample}_raw.bam",
    output:
        temp(multiext(
            "picard/{species}.{build}.{release}.{datatype}/stats/{sample}_raw",
            ".alignment_summary_metrics",
            ".insert_size_metrics",
            ".insert_size_histogram.pdf",
            ".base_distribution_by_cycle_metrics",
            ".base_distribution_by_cycle.pdf",
            ".gc_bias.detail_metrics",
            ".gc_bias.summary_metrics",
            ".gc_bias.pdf",
        )),
    log:
        "logs/picard/collectmultiplemetrics/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/picard/collectmultiplemetrics/{species}.{build}.{release}.{datatype}/{sample}.tsv",
    params:
        extra=config.get("params", {}).get("picard", {}).get("", "--METRIC_ACCUMULATION_LEVEL SAMPLE"),
    wrapper:
        f"{snakemake_wrappers_version}/bio/picard/collectmultiplemetrics"

use rule picard_create_multiple_metrics_raw as picard_create_multiple_metrics with:
    input:
        "results/sambamba/markdup/{species}.{build}.{release}.{datatype}/{sample}.bam",
    output:
        temp(multiext(
            "picard/{species}.{build}.{release}.{datatype}/stats/{sample}",
            ".alignment_summary_metrics",
            ".insert_size_metrics",
            ".insert_size_histogram.pdf",
            ".base_distribution_by_cycle_metrics",
            ".base_distribution_by_cycle.pdf",
            ".gc_bias.detail_metrics",
            ".gc_bias.summary_metrics",
            ".gc_bias.pdf",
        )),