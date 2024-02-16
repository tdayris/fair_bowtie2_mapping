module bowtie2_sambamba_metawrapper:
    meta_wrapper:
        "v3.3.6/meta/bio/bowtie2_sambamba"
    config:
        config


use rule bowtie2_build from bowtie2_sambamba_metawrapper as fair_bowtie2_mapping_bowtie2_build with:
    input:
        ref=getattr(
            lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}'",
                within=genomes,
            ),
            "dna_fasta",
            "reference/sequences/{species}.{build}.{release}.{datatype}.fasta",
        ),
    output:
        multiext(
            "reference/bowtie2_index/{species}.{build}.{release}.{datatype}",
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
        tmpdir="tmp",
    log:
        "logs/fair_bowtie2_mapping/bowtie2_build/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping/bowtie2_build/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=lookup(dpath="params/bowtie2/build", within=config),


use rule bowtie2_alignment from bowtie2_sambamba_metawrapper as fair_bowtie2_mapping_bowtie2_alignment with:
    input:
        sample=branch(
            lookup(
                query="sample_id == '{sample}' & species == '{species}' & build == '{build}' & release == '{release}' & downstream_file == downstream_file",
                within=samples,
            ),
            then=expand(
                "tmp/fair_bowtie2_mapping/fastp_trimming_pair_ended/{sample}.{stream}.fastq",
                sample="{sample}",
                stream=stream_list,
            ),
            otherwise=expand(
                "tmp/fair_bowtie2_mapping/fastp_trimming_single_ended/{sample}.fastq",
                sample="{sample}",
            ),
        ),
        idx=getattr(
            lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}'",
                within=genomes,
            ),
            "bowtie2_index",
            multiext(
                "reference/bowtie2_index/{species}.{build}.{release}.{datatype}",
                ".1.bt2",
                ".2.bt2",
                ".3.bt2",
                ".4.bt2",
                ".rev.1.bt2",
                ".rev.2.bt2",
            ),
        ),
    output:
        temp(
            "tmp/fair_bowtie2_mapping/bowtie2_alignment/{species}.{build}.{release}.{datatype}/{sample}.bam"
        ),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: (15 * 1024) * attempt,
        runtime=lambda wildcards, attempt: int(60 * 1.5) * attempt,
        tmpdir="tmp",
    log:
        "logs/fair_bowtie2_mapping/bowtie2_alignment/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping/bowtie2_alignment/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=config.get("params", {})
        .get("bowtie2", {})
        .get(
            "align",
            " --rg-id {sample} --rg 'SM:{sample} LB:{sample} PU:{species}.{build}.{release}.{datatype}.{sample} PL:ILLUMINA'",
        ),


use rule sambamba_sort from bowtie2_sambamba_metawrapper as fair_bowtie2_mapping_sambamba_sort with:
    input:
        "tmp/fair_bowtie2_mapping/bowtie2_alignment/{species}.{build}.{release}.{datatype}/{sample}.bam",
    output:
        temp(
            "tmp/fair_bowtie2_mapping/sambamba_sort/{species}.{build}.{release}.{datatype}/{sample}.bam"
        ),
    threads: 6
    resources:
        mem_mb=lambda wildcards, attempt: (6 * 1024) * attempt,
        runtime=lambda wildcards, attempt: int(60 * 0.75) * attempt,
        tmpdir="tmp",
    log:
        "logs/fair_bowtie2_mapping/sambamba_sort/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping/sambamba_sort/{species}.{build}.{release}.{datatype}/{sample}.tsv"


use rule sambamba_view from bowtie2_sambamba_metawrapper as fair_bowtie2_mapping_sambamba_view with:
    input:
        "tmp/fair_bowtie2_mapping/sambamba_sort/{species}.{build}.{release}.{datatype}/{sample}.bam",
    output:
        temp(
            "tmp/fair_bowtie2_mapping/sambamba_view/{species}.{build}.{release}.{datatype}/{sample}.bam"
        ),
    threads: 6
    resources:
        mem_mb=lambda wildcards, attempt: (2 * 1024) * attempt,
        runtime=lambda wildcards, attempt: int(60 * 0.75) * attempt,
        tmpdir="tmp",
    log:
        "logs/fair_bowtie2_mapping/sambamba_view/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping/sambamba_view/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup(dpath="params/sambamba/view", within=config),


use rule sambamba_markdup from bowtie2_sambamba_metawrapper as fair_bowtie2_mapping_sambamba_markdup with:
    input:
        "tmp/fair_bowtie2_mapping/sambamba_view/{species}.{build}.{release}.{datatype}/{sample}.bam",
    output:
        protected("results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam"),
    threads: 6
    resources:
        mem_mb=lambda wildcards, attempt: (6 * 1024) * attempt,
        runtime=lambda wildcards, attempt: int(60 * 0.75) * attempt,
        tmpdir="tmp",
    log:
        "logs/fair_bowtie2_mapping/sambamba_markdup/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping/sambamba_markdup/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup(dpath="params/sambamba/markdup", within=config),


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
        tmpdir="tmp",
    log:
        "logs/fair_bowtie2_mapping/sambamba_index/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping/sambamba_index/{species}.{build}.{release}.{datatype}/{sample}.tsv"
