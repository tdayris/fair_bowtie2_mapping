rule samtools_stats:
    input:
        bam="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
        bai="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai",
    output:
        temp("tmp/samtools/{species}.{build}.{release}.{datatype}/{sample}.txt"),
    log:
        "logs/samtools/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/samtools/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra="",
    wrapper:
        "v3.2.0/bio/samtools/stats"