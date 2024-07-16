module bowtie2_sambamba_metawrapper:
    meta_wrapper:
        f"{snakemake_wrappers_prefix}/meta/bio/bowtie2_sambamba"
    config:
        config


"""
## Memory
Requires a job with at most 5753.98  Mb,
 on average 4325.62 ± 2640.25 Mb, 
on Gustave Roussy's HPC Flamingo, on a 93.0  Mb dataset.
## Time
A job took 2:02:58 to proceed,
on average 0:33:06 ± 0:39:30
"""


use rule bowtie2_alignment from bowtie2_sambamba_metawrapper as fair_bowtie2_mapping_bowtie2_alignment with:
    input:
        sample=branch(
            lookup(
                query="sample_id == '{sample}' & species == '{species}' & build == '{build}' & release == '{release}' & downstream_file == downstream_file",
                within=samples,
            ),
            then=expand(
                "tmp/fair_bowtie2_mapping_fastp_trimming_pair_ended/{sample}.{stream}.fastq.gz",
                sample="{sample}",
                stream=stream_tuple,
            ),
            otherwise=expand(
                "tmp/fair_bowtie2_mapping_fastp_trimming_single_ended/{sample}.fastq.gz",
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
        mem_mb=lambda wildcards, attempt: 5_000 + 2_000 * attempt,
        runtime=lambda wildcards, attempt: 120 * attempt,
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
## Memory
Requires a job with at most 3363.84  Mb,
 on average 2259.72 ± 1341.36 Mb, 
on Gustave Roussy's HPC Flamingo, on a 93.0  Mb dataset.
## Time
A job took 0:20:26 to proceed,
on average 0:06:01 ± 0:06:23
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
        runtime=lambda wildcards, attempt: 30 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping_sambamba_sort/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_sambamba_sort/{species}.{build}.{release}.{datatype}/{sample}.tsv"


"""
## Memory
Requires a job with at most 803.36  Mb,
 on average 613.88 ± 351.0 Mb, 
on Gustave Roussy's HPC Flamingo, on a 93.0  Mb dataset.
## Time
A job took 0:10:44 to proceed,
on average 0:03:11 ± 0:03:18
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
        mem_mb=lambda wildcards, attempt: 700 + (200 * attempt),
        runtime=lambda wildcards, attempt: 15 * attempt,
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
## Memory
Requires a job with at most 3861.41  Mb,
 on average 2525.59 ± 1518.73 Mb, 
on Gustave Roussy's HPC Flamingo, on a 93.0  Mb dataset.
## Time
A job took 0:12:12 to proceed,
on average 0:03:20 ± 0:03:53
"""


use rule sambamba_markdup from bowtie2_sambamba_metawrapper as fair_bowtie2_mapping_sambamba_markdup with:
    input:
        "tmp/fair_bowtie2_mapping_sambamba_view/{species}.{build}.{release}.{datatype}/{sample}.bam",
    output:
        protected("results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam"),
    threads: 6
    resources:
        mem_mb=lambda wildcards, attempt: 3_000 + (1_000 * attempt),
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
## Memory
Requires a job with at most 468.46  Mb,
 on average 348.49 ± 185.48 Mb, 
on Gustave Roussy's HPC Flamingo, on a 93.0  Mb dataset.
## Time
A job took 0:03:41 to proceed,
on average 0:01:22 ± 0:01:09
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
