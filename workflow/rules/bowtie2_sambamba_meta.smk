module bowtie2_sambamba_metawrapper:
    meta_wrapper:
        f"{snakemake_wrappers_prefix}/meta/bio/bowtie2_sambamba"
    config:
        config


"""
Reported on Flamingo, on a 6Gb (HG38) dataset with 20 threads:
* time 1h27±15min
* mem 5746±20mb
"""


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
        mem_mb=lambda wildcards, attempt: 6_000 * attempt,
        runtime=lambda wildcards, attempt: int(60 * 1.75) * attempt,
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


"""
Reported on Flamingo on a 6gb hg38 dataset with 6 threads
* mem 3296±250mb
* time 5±12min
"""


use rule sambamba_sort from bowtie2_sambamba_metawrapper as fair_bowtie2_mapping_sambamba_sort with:
    input:
        "tmp/fair_bowtie2_mapping_bowtie2_alignment/{species}.{build}.{release}.{datatype}/{sample}.bam",
    output:
        temp(
            "tmp/fair_bowtie2_mapping_sambamba_sort/{species}.{build}.{release}.{datatype}/{sample}.bam"
        ),
    threads: 6
    resources:
        mem_mb=lambda wildcards, attempt: 3_300 + (700 * attempt),
        runtime=lambda wildcards, attempt: 20 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping_sambamba_sort/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_sambamba_sort/{species}.{build}.{release}.{datatype}/{sample}.tsv"


"""
Reported on Flamingo on a 6gb hg38 dataset with 6 threads
* mem 789±100mb
* time 1±2min
"""


use rule sambamba_view from bowtie2_sambamba_metawrapper as fair_bowtie2_mapping_sambamba_view with:
    input:
        "tmp/fair_bowtie2_mapping_sambamba_sort/{species}.{build}.{release}.{datatype}/{sample}.bam",
    output:
        temp(
            "tmp/fair_bowtie2_mapping_sambamba_view/{species}.{build}.{release}.{datatype}/{sample}.bam"
        ),
    threads: 6
    resources:
        mem_mb=lambda wildcards, attempt: 800 + (200 * attempt),
        runtime=lambda wildcards, attempt: 10 * attempt,
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


"""
Reported on Flamingo on a 6gb hg38 dataset with 6 threads
* time 2±2min
* mem 2,5Gb±1Gb
"""


use rule sambamba_markdup from bowtie2_sambamba_metawrapper as fair_bowtie2_mapping_sambamba_markdup with:
    input:
        "tmp/fair_bowtie2_mapping_sambamba_view/{species}.{build}.{release}.{datatype}/{sample}.bam",
    output:
        protected("results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam"),
    threads: 6
    resources:
        mem_mb=lambda wildcards, attempt: 3_000 + (600 * attempt),
        runtime=lambda wildcards, attempt: 20 * attempt,
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


"""
Reported on Flamingo on a 6gb hg38 dataset with 6 threads
* time 40±20sec
* mem 460±20mb
"""


use rule sambamba_index from bowtie2_sambamba_metawrapper as fair_bowtie2_mapping_sambamba_index with:
    input:
        "results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
    output:
        protected(
            "results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 500 + (100 * attempt),
        runtime=lambda wildcards, attempt: 5 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping_sambamba_index/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_sambamba_index/{species}.{build}.{release}.{datatype}/{sample}.tsv"
