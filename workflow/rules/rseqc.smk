"""
## Memory
Requires a job with at most 611.01  Mb,
 on average 459.18 ± 252.04 Mb, 
on Gustave Roussy's HPC Flamingo, on a 93.0  Mb dataset.
## Time
A job took 0:01:33 to proceed,
on average 0:00:22 ± 0:00:31
"""


rule fair_bowtie2_mapping_rseqc_infer_experiment:
    input:
        aln="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
        alnbai="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai",
        refgene=lambda wildcards: get_genepred_bed(wildcards),
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
        "v5.3.0/bio/rseqc/infer_experiment"


"""
## Memory
Requires a job with at most 583.43  Mb,
 on average 436.03 ± 240.22 Mb, 
on Gustave Roussy's HPC Flamingo, on a 93.0  Mb dataset.
## Time
A job took 0:09:02 to proceed,
on average 0:02:28 ± 0:02:52
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
        mem_mb=lambda wildcards, attempt: (attempt * 100) + 500,
        runtime=lambda wildcards, attempt: attempt * 15,
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
## Memory
Requires a job with at most 1986.88  Mb,
 on average 704.57 ± 564.98 Mb, 
on Gustave Roussy's HPC Flamingo, on a 93.0  Mb dataset.
## Time
A job took 0:11:13 to proceed,
on average 0:02:45 ± 0:03:35
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
        mem_mb=lambda wildcards, attempt: (attempt * 300) + 2_000,
        runtime=lambda wildcards, attempt: attempt * 20,
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
        "v5.3.0/bio/rseqc/read_gc"


"""
## Memory
Requires a job with at most 1564.66  Mb,
 on average 1094.74 ± 639.86 Mb, 
on Gustave Roussy's HPC Flamingo, on a 93.0  Mb dataset.
## Time
A job took 0:14:04 to proceed,
on average 0:04:06 ± 0:04:19
"""


rule fair_bowtie2_mapping_rseqc_read_distribution:
    input:
        aln="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
        alnbai="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai",
        refgene=lambda wildcards: get_genepred_bed(wildcards),
    output:
        temp(
            "tmp/fair_bowtie2_mapping_rseqc_read_distribution/{species}.{build}.{release}.{datatype}/{sample}.txt"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (attempt * 500) + 1_200,
        runtime=lambda wildcards, attempt: attempt * 20,
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
        "v5.3.0/bio/rseqc/read_distribution"


"""
## Memory
Requires a job with at most 1623.0  Mb,
 on average 1035.33 ± 661.93 Mb, 
on Gustave Roussy's HPC Flamingo, on a 93.0  Mb dataset.
## Time
A job took 0:02:53 to proceed,
on average 0:01:10 ± 0:00:49
"""


rule fair_bowtie2_mapping_rseqc_inner_distance:
    input:
        aln="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
        alnbai="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai",
        refgene=lambda wildcards: get_genepred_bed(wildcards),
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
        mem_mb=lambda wildcards, attempt: (attempt * 300) + 1_500,
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
