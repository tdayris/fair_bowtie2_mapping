rule fair_bowtie2_mapping_goleft_indexcov:
    input:
        aln="results/{species}.{build}.{release}.dna/Mapping/{sample}.bam",
        alnbai="results/{species}.{build}.{release}.dna/Mapping/{sample}.bam.bai",
        fasta=getattr(
            lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}'",
                within=genomes,
            ),
            "dna_fasta",
            "reference/sequences/{species}.{build}.{release}.dna.fasta",
        ),
        fai=getattr(
            lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}'",
                within=genomes,
            ),
            "dna_fai",
            "reference/sequences/{species}.{build}.{release}.dna.fasta.fai",
        ),
    output:
        ped=temp(
            "tmp/fair_bowtie2_mapping/goleft/indexcov/{species}.{release}.{build}/{sample}-indexcov.ped"
        ),
        roc=temp(
            "tmp/fair_bowtie2_mapping/goleft/indexcov/{species}.{release}.{build}/{sample}-indexcov.roc"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 30,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping/goleft/indexcov/{sample}.{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping/goleft/indexcov/{sample}.{species}.{build}.{release}.tsv"
    params:
        extra=dlookup(
            dpath="params/fair_bowtie2_mapping/goleft/indexcov",
            within=config,
            default="",
        ),
    # wrapper:
    #     f"{snakemake_wrappers_prefix}/bio/goleft/indexcov"
    conda:
        "../envs/goleft.yaml"
    script:
        "../scripts/goleft_indexcov.py"
