rule fair_bowtie2_mapping_samtools_stats:
    input:
        bam="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
        bai="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai",
    output:
        temp(
            "tmp/fair_bowtie2_mapping/samtools_stats/{species}.{build}.{release}.{datatype}/{sample}.txt"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (2 * 1024) * attempt,
        runtime=lambda wildcards, attempt: int(60 * 0.6) * attempt,
        tmpdir="tmp",
    log:
        "logs/fair_bowtie2_mapping/samtools_stats/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping/samtools_stats/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup(dpath="params/samtools/stats", within=config),
    wrapper:
        "v3.3.3/bio/samtools/stats"
