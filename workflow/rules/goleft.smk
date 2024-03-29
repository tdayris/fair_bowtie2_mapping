rule fair_bowtie2_mapping_goleft_indexcov:
    input:
        aln="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
        alnbai="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai",
        fasta=lambda wildcards: select_fasta(wildcards),
        fai=lambda wildcards: select_fai(wildcards),
    output:
        ped=temp(
            "tmp/fair_bowtie2_mapping/goleft/indexcov/{species}.{release}.{build}.{datatype}/{sample}-indexcov.ped"
        ),
        roc=temp(
            "tmp/fair_bowtie2_mapping/goleft/indexcov/{species}.{release}.{build}.{datatype}/{sample}-indexcov.roc"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 30,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping/goleft/indexcov/{sample}.{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping/goleft/indexcov/{sample}.{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_bowtie2_mapping/goleft/indexcov",
            default="",
        ),
    conda:
        "../envs/goleft.yaml"
    script:
        "../scripts/goleft_indexcov.py"
