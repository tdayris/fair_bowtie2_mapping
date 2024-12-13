"""
## Memory
Requires a job with at most 1581.02  Mb,
 on average 1197.14 ± 710.57 Mb, 
on Gustave Roussy's HPC Flamingo, on a 93.0  Mb dataset.
## Time
A job took 0:01:33 to proceed,
on average 0:00:16 ± 0:00:31
"""


rule fair_bowtie2_mapping_goleft_indexcov:
    input:
        aln="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
        alnbai="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai",
        fasta=lambda wildcards: select_fasta(wildcards),
        fai=lambda wildcards: select_fai(wildcards),
    output:
        ped=temp(
            "tmp/fair_bowtie2_mapping_goleft_indexcov/{species}.{release}.{build}.{datatype}/{sample}-indexcov.ped"
        ),
        roc=temp(
            "tmp/fair_bowtie2_mapping_goleft_indexcov/{species}.{release}.{build}.{datatype}/{sample}-indexcov.roc"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (attempt * 100) + 1_500,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping_goleft_indexcov/{sample}.{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_goleft_indexcov/{sample}.{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_bowtie2_mapping_goleft_indexcov",
            default="",
        ),
    wrapper:
        "v5.5.0/bio/goleft/indexcov"
