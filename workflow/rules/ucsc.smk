"""
Reported on Flamingo on a HG38
* time 5s ± 2s
* mem 450mb ± 50mb
"""


rule fair_bowtie2_mapping_ucsc_genepred_to_bed:
    input:
        "reference/annotation/{species}.{build}.{release}/{species}.{build}.{release}.genePred",
    output:
        temp(
            "tmp/fair_bowtie2_mapping_ucsc_genepred_to_bed/{species}.{build}.{release}.bed"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (attempt * 100) + 500,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping_ucsc_genepred_to_bed/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_ucsc_genepred_to_bed/{species}.{build}.{release}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_bowtie2_mapping_ucsc_genepred2bed",
            default="-tab",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/ucsc/genePredToBed"
