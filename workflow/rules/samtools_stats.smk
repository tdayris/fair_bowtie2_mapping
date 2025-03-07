"""
## Memory
Requires a job with at most 386.02  Mb,
 on average 285.2 ± 149.11 Mb, 
on Gustave Roussy's HPC Flamingo, on a 93.0  Mb dataset.
## Time
A job took 0:07:12 to proceed,
on average 0:02:00 ± 0:02:16
"""


rule fair_bowtie2_mapping_samtools_stats:
    input:
        bam="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
        bai="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai",
    output:
        temp(
            "tmp/fair_bowtie2_mapping_samtools_stats/{species}.{build}.{release}.{datatype}/{sample}.txt"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 300 + (100 * attempt),
        runtime=lambda wildcards, attempt: 10 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping_samtools_stats/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_samtools_stats/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_bowtie2_mapping_samtools_stats",
            default="",
        ),
    wrapper:
        "v5.8.3/bio/samtools/stats"


"""
No data
"""


rule fair_bowtie2_mapping_samtools_idxstats:
    input:
        bam="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
        idx="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai",
    output:
        temp(
            "tmp/fair_bowtie2_mapping_samtools_idxstats/{species}.{build}.{release}.{datatype}/{sample}.idxstats"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 300 + (100 * attempt),
        runtime=lambda wildcards, attempt: 10 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping_samtools_idxstats/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_samtools_idxstats/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_bowtie2_mapping_samtools_idxstats",
            default="",
        ),
    wrapper:
        "v5.8.3/bio/samtools/idxstats"
