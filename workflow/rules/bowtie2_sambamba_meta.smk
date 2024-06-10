module bowtie2_sambamba_metawrapper:
    meta_wrapper:
        f"{snakemake_wrappers_prefix}/meta/bio/bowtie2_sambamba"
    config:
        config


use rule bowtie2_build from bowtie2_sambamba_metawrapper as fair_bowtie2_mapping_bowtie2_build with:
    input:
        ref=lambda wildcards: select_fasta(wildcards),
    output:
        multiext(
            "reference/bowtie2_index/{species}.{build}.{release}/{species}.{build}.{release}.{datatype}",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: (15 * 1024) * attempt,
        runtime=lambda wildcards, attempt: int(60 * 1.5) * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mappingi_bowtie2_build/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_bowtie2_build/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_bowtie2_mapping_bowtie2_build",
            default="",
        ),


use rule bowtie2_alignment from bowtie2_sambamba_metawrapper as fair_bowtie2_mapping_bowtie2_alignment with:
    input:
        sample=branch(
            lookup(
                query="sample_id == '{sample}' & species == '{species}' & build == '{build}' & release == '{release}' & downstream_file == downstream_file",
                within=samples,
            ),
            then=expand(
                "tmp/fair_bowtie2_mapping_fastp_trimming_pair_ended/{sample}.{stream}.fastq",
                sample="{sample}",
                stream=stream_tuple,
            ),
            otherwise=expand(
                "tmp/fair_bowtie2_mapping_fastp_trimming_single_ended/{sample}.fastq",
                sample="{sample}",
            ),
        ),
        idx=lambda wildcards: branch(
            condition=str(wildcards.datatype).lower(),
            cases={
                "cdna": get_cdna_bowtie2_index(wildcards),
                "dna": get_dna_bowtie2_index(wildcards),
                "transcripts": get_transcripts_bowtie2_index(wildcards),
            },
        ),
    output:
        temp(
            "tmp/fair_bowtie2_mapping_bowtie2_alignment/{species}.{build}.{release}.{datatype}/{sample}.bam"
        ),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: (15 * 1024) * attempt,
        runtime=lambda wildcards, attempt: int(60 * 1.5) * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping_bowtie2_alignment/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_bowtie2_alignment/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_bowtie2_mapping_bowtie2_align",
            default=" --rg-id {sample} --rg 'SM:{sample} LB:{sample} PU:{species}.{build}.{release}.{datatype}.{sample} PL:ILLUMINA'",
        ),


use rule sambamba_sort from bowtie2_sambamba_metawrapper as fair_bowtie2_mapping_sambamba_sort with:
    input:
        "tmp/fair_bowtie2_mapping_bowtie2_alignment/{species}.{build}.{release}.{datatype}/{sample}.bam",
    output:
        temp(
            "tmp/fair_bowtie2_mapping_sambamba_sort/{species}.{build}.{release}.{datatype}/{sample}.bam"
        ),
    threads: 6
    resources:
        mem_mb=lambda wildcards, attempt: (6 * 1024) * attempt,
        runtime=lambda wildcards, attempt: int(60 * 0.75) * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping_sambamba_sort/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_sambamba_sort/{species}.{build}.{release}.{datatype}/{sample}.tsv"


use rule sambamba_view from bowtie2_sambamba_metawrapper as fair_bowtie2_mapping_sambamba_view with:
    input:
        "tmp/fair_bowtie2_mapping_sambamba_sort/{species}.{build}.{release}.{datatype}/{sample}.bam",
    output:
        temp(
            "tmp/fair_bowtie2_mapping_sambamba_view/{species}.{build}.{release}.{datatype}/{sample}.bam"
        ),
    threads: 6
    resources:
        mem_mb=lambda wildcards, attempt: (2 * 1024) * attempt,
        runtime=lambda wildcards, attempt: int(60 * 0.75) * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping_sambamba_view/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_sambamba_view/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_bowtie2_mapping_sambamba_view",
            default="--format 'bam' --filter 'mapping_quality >= 30 and not (unmapped or mate_is_unmapped)' ",
        ),


use rule sambamba_markdup from bowtie2_sambamba_metawrapper as fair_bowtie2_mapping_sambamba_markdup with:
    input:
        "tmp/fair_bowtie2_mapping_sambamba_view/{species}.{build}.{release}.{datatype}/{sample}.bam",
    output:
        protected("results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam"),
    threads: 6
    resources:
        mem_mb=lambda wildcards, attempt: (6 * 1024) * attempt,
        runtime=lambda wildcards, attempt: int(60 * 0.75) * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping_sambamba_markdup/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_sambamba_markdup/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_bowtie2_mapping_sambamba_markdup",
            default="--remove-duplicates --overflow-list-size=500000",
        ),


use rule sambamba_index from bowtie2_sambamba_metawrapper as fair_bowtie2_mapping_sambamba_index with:
    input:
        "results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
    output:
        protected(
            "results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (2 * 1024) * attempt,
        runtime=lambda wildcards, attempt: int(60 * 0.5) * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping_sambamba_index/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_sambamba_index/{species}.{build}.{release}.{datatype}/{sample}.tsv"
