rule samtools_stats:
    input:
        bam="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
        bai="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai",
    output:
        temp("tmp/samtools/{species}.{build}.{release}.{datatype}/{sample}.txt"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (9 * 1024) * attempt,
        runtime=lambda wildcards, attempt: int(60 * 0.6) * attempt,
        tmpdir="tmp",
    log:
        "logs/samtools/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/samtools/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra="",
    wrapper:
        "v3.2.0/bio/samtools/stats"
