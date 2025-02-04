rule fair_bowtie2_mapping_deeptools_alignment_sieve:
    input:
        aln="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
        aln_idx="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai",
        blacklist=lambda wildcards: lookup_genomes(
            wildcards, key="blacklist", default=[]
        ),
    output:
        temp(
            "tmp/fair_bowtie2_mapping_deeptools_alignment_sieve/{species}.{build}.{release}.{datatype}/{sample}.bam"
        ),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1_000 + 6_000,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping_deeptools_alignment_sieve/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_deeptools_alignment_sieve/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_bowtie2_mapping_deeptools_alignment_sieve", default=""
        ),
    wrapper:
        "v5.6.0/bio/deeptools/alignmentsieve"


use rule fair_bowtie2_mapping_sambamba_sort as fair_bowtie2_mapping_sambamba_sort_sieve with:
    input:
        "tmp/fair_bowtie2_mapping_deeptools_alignment_sieve/{species}.{build}.{release}.{datatype}/{sample}.bam",
    output:
        temp(
            "tmp/fair_bowtie2_mapping_sambamba_sort_sieve/{species}.{build}.{release}.{datatype}/{sample}.bam"
        ),
    log:
        "logs/fair_bowtie2_mapping_sambamba_sort_sieve/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_sambamba_sort_sieve/{species}.{build}.{release}.{datatype}/{sample}.tsv"


use rule fair_bowtie2_mapping_sambamba_index as fair_bowtie2_mapping_sambamba_sieve_index with:
    input:
        "tmp/fair_bowtie2_mapping_sambamba_sort_sieve/{species}.{build}.{release}.{datatype}/{sample}.bam",
    output:
        temp(
            "tmp/fair_bowtie2_mapping_sambamba_sort_sieve/{species}.{build}.{release}.{datatype}/{sample}.bam.bai"
        ),
    log:
        "logs/fair_bowtie2_mapping_sambamba_sieve_index/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_sambamba_sieve_index/{species}.{build}.{release}.{datatype}/{sample}.tsv"
