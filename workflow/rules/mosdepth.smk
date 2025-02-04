rule fair_bowtie2_mapping_mosdepth:
    input:
        bam="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
        bai="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai",
        fasta=lambda wildcards: select_fasta(wildcards),
        fai=lambda wildcards: select_fai(wildcards),
        bed=lambda wildcards: get_intervals(wildcards),
    output:
        temp(
            "tmp/fair_bowtie2_mapping_mosdepth/{species}.{build}.{release}.{datatype}/{sample}.mosdepth.global.dist.txt"
        ),
        temp(
            "tmp/fair_bowtie2_mapping_mosdepth/{species}.{build}.{release}.{datatype}/{sample}.mosdepth.region.dist.txt"
        ),
        temp(
            "tmp/fair_bowtie2_mapping_mosdepth/{species}.{build}.{release}.{datatype}/{sample}.regions.bed.gz"
        ),
        summary=temp(
            "tmp/fair_bowtie2_mapping_mosdepth/{species}.{build}.{release}.{datatype}/{sample}.mosdepth.summary.txt"
        ),
    threads: 10
    resources:
        mem_mb=lambda wildcards, attempt: (2_000 * attempt) + 3_000,
        runtime=lambda wildcards, attempt: 15 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping_mosdepth/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_mosdepth/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_bowtie2_mapping_mosdepth",
            default="--no-per-base --use-median",
        ),
    wrapper:
        "v5.6.0/bio/mosdepth"
