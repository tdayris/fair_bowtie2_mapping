rule samtools_stats:
    input:
        bam="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
        bai="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai",
    output:
        temp("tmp/samtools/{species}.{build}.{release}.{datatype}/{sample}.txt"),
    threads: 1
    resources:
        # Reserve 2Gb per attempt
        mem_mb=lambda wildcards, attempt: (2 * 1024) * attempt,
        # Reserve 35min per attempt
        runtime=lambda wildcards, attempt: int(60 * 0.6) * attempt,
        tmpdir="tmp",
    log:
        "logs/samtools/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/samtools/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=config.get("params", {}).get("samtools", ""),
    wrapper:
        "v3.3.3/bio/samtools/stats"
