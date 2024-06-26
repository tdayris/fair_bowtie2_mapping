rule fair_bowtie2_mapping_rseqc_infer_experiment:
    input:
        aln="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
        alnbai="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai",
        refgene="tmp/fair_bowtie2_mapping_ucsc_genepred_to_bed/{species}.{build}.{release}.bed",
    output:
        temp(
            "tmp/fair_bowtie2_mapping_rseqc_infer_experiment/{species}.{build}.{release}.{datatype}/{sample}.infer_experiment.txt"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping_rseqc_infer_experiment/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_rseqc_infer_experiment/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_bowtie2_mapping_rseqc_infer_experiment",
            default="",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/rseqc/infer_experiment"


rule fair_bowtie2_mapping_rseqc_bamstat:
    input:
        aln="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
        alnbai="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai",
    output:
        temp(
            "tmp/fair_bowtie2_mapping_rseqc_bamstat/{species}.{build}.{release}.{datatype}/{sample}.bamstat.txt"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping_rseqc_bamstat/{species}.{build}.{release}/{sample}.{datatype}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_rseqc_bamstat/{species}.{build}.{release}/{sample}.{datatype}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_bowtie2_mapping_rseqc_bamstat",
            default="",
        ),
    conda:
        "../envs/rseqc.yaml"
    script:
        "../scripts/rseqc_bamstat.py"


rule fair_bowtie2_mapping_rseqc_read_gc:
    input:
        aln="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
        alnbai="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai",
    output:
        xls=temp(
            "tmp/fair_bowtie2_mapping_rseqc_read_gc/{species}.{build}.{release}.{datatype}/{sample}.GC.xls"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping_rseqc_read_gc/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_rseqc_read_gc/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_bowtie2_mapping_rseqc_read_gc",
            default="",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/rseqc/read_gc"


rule fair_bowtie2_mapping_rseqc_read_distribution:
    input:
        aln="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
        alnbai="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai",
        refgene="tmp/fair_bowtie2_mapping_ucsc_genepred_to_bed/{species}.{build}.{release}.bed",
    output:
        temp(
            "tmp/fair_bowtie2_mapping_rseqc_read_distribution/{species}.{build}.{release}.{datatype}/{sample}.txt"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping_rseqc_read_distribution/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_rseqc_read_distribution/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_bowtie2_mapping_rseqc_read_distribution",
            default="",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/rseqc/read_distribution"


rule fair_bowtie2_mapping_rseqc_inner_distance:
    input:
        aln="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
        alnbai="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai",
        refgene="tmp/fair_bowtie2_mapping_ucsc_genepred_to_bed/{species}.{build}.{release}.bed",
    output:
        txt=temp(
            "tmp/fair_bowtie2_mapping_rseqc_inner_distance/{species}.{build}.{release}.{datatype}/{sample}.inner_distance.txt"
        ),
        freq=temp(
            "tmp/fair_bowtie2_mapping_rseqc_inner_distance/{species}.{build}.{release}.{datatype}/{sample}.inner_distance_freq.txt"
        ),
        plot_r=temp(
            "tmp/fair_bowtie2_mapping_rseqc_inner_distance/{species}.{build}.{release}.{datatype}/{sample}.inner_distance_plot.r"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping_rseqc_inner_distance/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_rseqc_inner_distance/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_bowtie2_mapping_rseqc_inner_distance",
            default="",
        ),
    conda:
        "../envs/rseqc.yaml"
    script:
        "../scripts/rseqc_inner_distance.py"
