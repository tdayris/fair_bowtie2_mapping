"""
Reported on Flamingo on a 6gb dataset (hg38)
* mem 1.5Go ± 1mb
* time 20s ± 20s
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
        mem_mb=lambda wildcards, attempt: (attempt * 100) + 1_590,
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
        f"{snakemake_wrappers_prefix}/bio/goleft/indexcov"
