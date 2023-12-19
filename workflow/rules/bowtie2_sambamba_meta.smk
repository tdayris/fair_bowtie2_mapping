module bowtie2_sambamba:
    meta_wrapper:
        "v3.2.0/meta/bio/bowtie2_sambamba"
    config:
        config


use rule bowtie2_build from bowtie2_sambamba with:
    input:
        unpack(get_bowtie2_build_input),
    output:
        multiext(
            "reference/{species}.{build}.{release}.{datatype}",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: (48 * 1024) * attempt,
        runtime=lambda wildcards, attempt: int(60 * 1.5) * attempt,
        tmpdir="tmp",
    log:
        "logs/bowtie2/build/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/bowtie2/build/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=config.get("params", {}).get("bowtie2", {}).get("build", ""),


use rule bowtie2_alignment from bowtie2_sambamba with:
    input:
        unpack(get_bowtie2_alignment_input),
    output:
        temp("tmp/bowtie2/{species}.{build}.{release}.{datatype}/{sample}_raw.bam"),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: (48 * 1024) * attempt,
        runtime=lambda wildcards, attempt: int(60 * 1.5) * attempt,
        tmpdir="tmp",
    log:
        "logs/bowtie2/align/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/bowtie2/align/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=config.get("params", {})
        .get("bowtie2", {})
        .get(
            "align",
            " --rg-id {sample} --rg 'SM:{sample} LB:{sample} PU:{species}.{build}.{release}.{datatype}.{sample} PL:ILLUMINA'",
        ),


use rule sambamba_sort from bowtie2_sambamba with:
    input:
        "tmp/bowtie2/{species}.{build}.{release}.{datatype}/{sample}_raw.bam",
    output:
        temp("tmp/sambamba/sort/{species}.{build}.{release}.{datatype}/{sample}.bam"),
    threads: 6
    resources:
        mem_mb=lambda wildcards, attempt: (15 * 1024) * attempt,
        runtime=lambda wildcards, attempt: int(60 * 0.75) * attempt,
        tmpdir="tmp",
    log:
        "logs/sambamba/sort/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/sambamba/sort/{species}.{build}.{release}.{datatype}/{sample}.tsv"


use rule sambamba_view from bowtie2_sambamba with:
    input:
        "tmp/sambamba/sort/{species}.{build}.{release}.{datatype}/{sample}.bam",
    output:
        temp("tmp/sambamba/view/{species}.{build}.{release}.{datatype}/{sample}.bam"),
    threads: 6
    resources:
        mem_mb=lambda wildcards, attempt: (6 * 1024) * attempt,
        runtime=lambda wildcards, attempt: int(60 * 0.75) * attempt,
        tmpdir="tmp",
    log:
        "logs/sambamba/view/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/sambamba/view/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=config.get("params", {})
        .get("sambamba", {})
        .get(
            "view",
            "--format 'bam' --filter 'mapping_quality >= 30 and not (unmapped or mate_is_unmapped)'",
        ),


use rule sambamba_markdup from bowtie2_sambamba with:
    input:
        "tmp/sambamba/view/{species}.{build}.{release}.{datatype}/{sample}.bam",
    output:
        protected("results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam"),
    threads: 6
    resources:
        mem_mb=lambda wildcards, attempt: (6 * 1024) * attempt,
        runtime=lambda wildcards, attempt: int(60 * 0.75) * attempt,
        tmpdir="tmp",
    log:
        "logs/sambamba/markdup/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/sambamba/markdup/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=config.get("params", {})
        .get("sambamba", {})
        .get("markdup", "--remove-duplicates"),


use rule sambamba_index from bowtie2_sambamba with:
    input:
        "results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
    output:
        protected(
            "results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (4 * 1024) * attempt,
        runtime=lambda wildcards, attempt: int(60 * 0.5) * attempt,
        tmpdir="tmp",
    log:
        "logs/sambamba/index/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/sambamba/index/{species}.{build}.{release}.{datatype}/{sample}.tsv"
