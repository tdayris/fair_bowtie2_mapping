rule fair_bowtie2_mapping_picard_create_multiple_metrics:
    input:
        bam="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
        bai="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai",
        ref=lambda wildcards: select_fasta(wildcards),
    output:
        temp(
            multiext(
                "tmp/fair_bowtie2_mapping_picard_create_multiple_metrics/{species}.{build}.{release}.{datatype}/stats/{sample}",
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
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (3 * 1024) * attempt,
        runtime=lambda wildcards, attempt: int(60 * 0.6) * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping_picard_create_multiple_metrics/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_picard_create_multiple_metrics/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_bowtie2_mapping_picard_collectmultiplemetrics",
            default="",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/picard/collectmultiplemetrics"
