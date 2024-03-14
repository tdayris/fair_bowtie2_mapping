rule fair_bowtie2_mapping_ucsc_genepred_to_bed:
    input:
        "reference/annotation/{species}.{build}.{release}.genePred",
    output:
        temp(
            "tmp/fair_bowtie2_mapping/ucsc_genepred_to_bed/{species}.{build}.{release}.bed"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping/ucsc_genepred_to_bed/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping/ucsc_genepred_to_bed/{species}.{build}.{release}.tsv"
    params:
        extra=lookup(
            dpath="params/fair_bowtie2_mapping/ucsc/genepred2bed", within=config
        ),
    wrapper:
        "v3.5.0/bio/ucsc/genePredToBed"
