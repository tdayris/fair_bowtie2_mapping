rule fair_bowtie2_mapping_rseqc_infer_experiment:
    input:
        aln="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
        alnbai="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai",
        refgene="tmp/fair_bowtie2_mapping/ucsc_genepred_to_bed/{species}.{build}.{release}.bed",
    output:
        temp(
            "tmp/fair_bowtie2_mapping/rseqc_infer_experiment/{species}.{build}.{release}.{datatype}/{sample}.infer_experiment.txt"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping/rseqc_infer_experiment/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping/rseqc_infer_experiment/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup(
            dpath="params/fair_bowtie2_mapping/rseqc/infer_experiment", within=config
        ),
    wrapper:
        "v3.5.0/bio/rseqc/infer_experiment"


rule fair_bowtie2_mapping_rseqc_bamstat:
    input:
        aln="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
        alnbai="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai",
    output:
        temp(
            "tmp/fair_bowtie2_mapping/rseqc_bamstat/{species}.{build}.{release}.{datatype}/{sample}.bamstat.txt"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping/rseqc_bamstat/{species}.{build}.{release}/{sample}.{datatype}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping/rseqc_bamstat/{species}.{build}.{release}/{sample}.{datatype}.tsv"
    params:
        extra=lookup(dpath="params/fair_bowtie2_mapping/rseqc/bamstat", within=config),
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
            "tmp/fair_bowtie2_mapping/rseqc_read_gc/{species}.{build}.{release}.{datatype}/{sample}.GC.xls"
        ),
        # plot_r=temp(
        #     "tmp/fair_bowtie2_mapping/rseqc_read_gc/{species}.{build}.{release}.{datatype}/{sample}.GC_plot.r"
        # ),
        # plot_pdf="results/{species}.{build}.{release}.{datatype}/RSeQC/{sample}.GC_plot.pdf",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping/rseqc_bamstat/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping/rseqc_bamstat/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup(dpath="params/fair_bowtie2_mapping/rseqc/read_gc", within=config),
    wrapper:
        "v3.5.0/bio/rseqc/read_gc"


rule fair_bowtie2_mapping_rseqc_read_distribution:
    input:
        aln="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
        alnbai="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai",
        refgene="tmp/fair_bowtie2_mapping/ucsc_genepred_to_bed/{species}.{build}.{release}.bed",
    output:
        temp(
            "tmp/fair_bowtie2_mapping/rseqc_read_distribution/{species}.{build}.{release}.{datatype}/{sample}.txt"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping/rseqc_read_distribution/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping/rseqc_read_distribution/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup(
            dpath="params/fair_bowtie2_mapping/rseqc/read_distribution", within=config
        ),
    wrapper:
        "v3.5.0/bio/rseqc/read_distribution"


rule fair_bowtie2_mapping_rseqc_inner_distance:
    input:
        aln="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
        alnbai="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai",
        refgene="tmp/fair_bowtie2_mapping/ucsc_genepred_to_bed/{species}.{build}.{release}.bed",
    output:
        txt=temp(
            "tmp/fair_bowtie2_mapping/rseqc_inner_distance/{species}.{build}.{release}.{datatype}/{sample}.inner_distance.txt"
        ),
        freq=temp(
            "tmp/fair_bowtie2_mapping/rseqc_inner_distance/{species}.{build}.{release}.{datatype}/{sample}.inner_distance_freq.txt"
        ),
        plot_r=temp(
            "tmp/fair_bowtie2_mapping/rseqc_inner_distance/{species}.{build}.{release}.{datatype}/{sample}.inner_distance_plot.r"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping/rseqc_inner_distance/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping/rseqc_inner_distance/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup(
            dpath="params/fair_bowtie2_mapping/rseqc/inner_distance", within=config
        ),
    conda:
        "../envs/rseqc.yaml"
    script:
        "../scripts/rseqc_inner_distance.py"
