"""
Reported on Flamingo on a 6gb dataset (hg38)
* mem 300mb ± 400mb
* time 1min ± 15s
"""


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
        mem_mb=lambda wildcards, attempt: (attempt * 100) + 600,
        runtime=lambda wildcards, attempt: attempt * 10,
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


"""
Reported on Flamingo on a 6gb dataset (hg38)
* time 1:50 ± 1min
* mem 580mb ± 50mb
"""


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
        mem_mb=lambda wildcards, attempt: (attempt * 100) + 600,
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


"""
Reported on Flamingo on a 6gb dataset (hg38)
* time 1:20 ± 2min
* mem 800mb ± 300mb
"""


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
        mem_mb=lambda wildcards, attempt: (attempt * 300) + 1_100,
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


"""
Reported on Flamingo on a 6gb dataset
* mem 1477mb ± 20mb
* time 3min ± 2min
"""


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
        mem_mb=lambda wildcards, attempt: (attempt * 500) + 1_200,
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


"""
Reported in Flamingo on a 6gb dataset (hg38)
* mem 1.5Gb ± 300mb
* time 2min ± 1min
"""


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
        mem_mb=lambda wildcards, attempt: (attempt * 300) + 1_700,
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
