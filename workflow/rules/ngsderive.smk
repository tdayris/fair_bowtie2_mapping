
rule ngsderive_grep_out_gtf_comments:
    input:
        lambda wildcards: get_gtf(wildcards),
    output:
        pipe(
            "tmp/fair_bowtie2_mapping/ngsderive/{species}.{build}.{release}.no_comment.gtf"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 25,
    log:
        "logs/fair_bowtie2_mapping/ngsderive/grep_out/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping/ngsderive/grep_out/{species}.{build}.{release}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_bowtie2_mapping/ngsderive/grep_out_gtf_comments",
            default="-P '^#'",
        ),
    conda:
        "../envs/bash.yaml"
    shell:
        "grep {params.extra} {input} > {output} 2> {log}"


rule ngsderive_sort_gtf:
    input:
        "tmp/fair_bowtie2_mapping/ngsderive/{species}.{build}.{release}.no_comment.gtf",
    output:
        temp(
            "tmp/fair_bowtie2_mapping/ngsderive/{species}.{build}.{release}.sorted.gtf"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 25,
    log:
        "logs/fair_bowtie2_mapping/ngsderive/sort_gtf/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping/ngsderive/sort_gtf/{species}.{build}.{release}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_bowtie2_mapping/ngsderive/sort_gtf",
            default="-k1,1 -k3,3",
        ),
    conda:
        "../envs/bash.yaml"
    shell:
        "sort {params.extra} {input} > {output} 2> {log}"


rule ngsderive_zip_sorted_gtf:
    input:
        "tmp/fair_bowtie2_mapping/ngsderive/{species}.{build}.{release}.sorted.gtf",
    output:
        temp(
            "tmp/fair_bowtie2_mapping/ngsderive/{species}.{build}.{release}.sorted.gtf.gz"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 25,
    log:
        "logs/fair_bowtie2_mapping/ngsderive/gzip/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping/ngsderive/gzip/{species}.{build}.{release}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_bowtie2_mapping/ngsderive/compress_sorted_gtf",
            default="-c",
        ),
    conda:
        "../envs/bash.yaml"
    shell:
        "pbgzip {params.extra} {input} > {output} 2> {log}"


rule tabix_gzipped_gtf:
    input:
        "tmp/fair_bowtie2_mapping/ngsderive/{species}.{build}.{release}.sorted.gtf.gz",
    output:
        temp(
            "tmp/fair_bowtie2_mapping/ngsderive/{species}.{build}.{release}.sorted.gtf.gz.tbi"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 25,
    log:
        "logs/fair_bowtie2_mapping/ngsderive/tabix/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping/ngsderive/tabix/{species}.{build}.{release}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_bowtie2_mapping_ngsderive_tabix_gtf",
            default="-p gff",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/tabix/index"


rule fair_bowtie2_mapping_ngsderive_strandedness:
    input:
        ngs="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
        ngs_bai="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai",
        gene_model="tmp/fair_bowtie2_mapping_ngsderive_{species}.{build}.{release}.sorted.gtf.gz",
        gene_model_tbi="tmp/fair_bowtie2_mapping_ngsderive_{species}.{build}.{release}.sorted.gtf.gz.tbi",
    output:
        tsv=temp(
            "tmp/fair_bowtie2_mapping_ngsderive_strandedness/{species}.{build}/{release}.{datatype}/{sample}.strandedness.tsv"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 25,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping_ngsderive_strandedness/{species}.{build}.{release}.{datatype}/{sample}.strandedness.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_ngsderive_strandedness/{species}.{build}.{release}.{datatype}/{sample}.strandedness.tsv"
    params:
        command="strandedness",
        extra=lookup_config(
            dpath="params/fair_bowtie2_mapping_ngsderive_strandedness",
            default="",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/ngsderive"


"""
Reported on Flamingo on a 6gb dataset (hg38)
* mem 832mb ± 50mb
* time 4min ± 2min
"""


rule fair_bowtie2_mapping_ngsderive_encoding:
    input:
        ngs="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
        ngs_bai="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai",
    output:
        tsv=temp(
            "tmp/fair_bowtie2_mapping_ngsderive_encoding/{species}.{build}/{release}.{datatype}/{sample}.encoding.tsv"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (attempt * 100) + 850,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping_ngsderive_encoding/{species}.{build}.{release}.{datatype}/{sample}.encoding.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_ngsderive_encoding/{species}.{build}.{release}.{datatype}/{sample}.encoding.tsv"
    params:
        command="encoding",
        extra=lookup_config(
            dpath="params/fair_bowtie2_mapping_ngsderive_encoding",
            default="",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/ngsderive"


"""
Reported on Flamingo on a 6gb dataset (hg38)
* time 1:20 ± 1min
* mem 700mb ± 200mb
"""


rule fair_bowtie2_mapping_ngsderive_instrument:
    input:
        ngs="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
        ngs_bai="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai",
    output:
        tsv=temp(
            "tmp/fair_bowtie2_mapping_ngsderive_instrument/{species}.{build}/{release}.{datatype}/{sample}.instrument.tsv"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (attempt * 100) + 900,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping_ngsderive_instrument/{species}.{build}.{release}.{datatype}/{sample}.instrument.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_ngsderive_instrument/{species}.{build}.{release}.{datatype}/{sample}.instrument.tsv"
    params:
        command="instrument",
        extra=lookup_config(
            dpath="params/fair_bowtie2_mapping_ngsderive_instrument",
            default="",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/ngsderive"


"""
Reported on Flamingo on a 6gb dataset
* mem 840mb ± 50mb
* time 2min ± 1min
"""


rule fair_bowtie2_mapping_ngsderive_readlen:
    input:
        ngs="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
        ngs_bai="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai",
    output:
        tsv=temp(
            "tmp/fair_bowtie2_mapping_ngsderive_readlen/{species}.{build}/{release}.{datatype}/{sample}.readlen.tsv"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (attempt * 100) + 875,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping_ngsderive_readlen/{species}.{build}.{release}.{datatype}/{sample}.readlen.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_ngsderive_readlen/{species}.{build}.{release}.{datatype}/{sample}.readlen.tsv"
    params:
        command="readlen",
        extra=lookup_config(
            dpath="params/fair_bowtie2_mapping_ngsderive_readlen",
            default="",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/ngsderive"
