rule fair_bowtie2_mapping_mtnucratiocalculator:
    input:
        "results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
    output:
        txt=report(
            "results/{species}.{build}.{release}.{datatype}/QC/Mitochondrion/{sample}.ratio.txt",
            caption="../report/mtnucratiocalculator.rst",
            category="Quality Controls",
            subcategory="Mapping",
            labels={
                "report": "txt",
                "step": "Mapping",
                "organism": "{species}.{build}.{release}",
            },
        ),
        json=temp(
            "tmp/fair_bowtie2_mapping_mtnucratiocalculator/{species}.{build}.{release}.{datatype}/{sample}.mtnuc.json"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1_000,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping_mtnucratiocalculator/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_mtnucratiocalculator/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        chrom=lookup_config(
            dpath="params/fair_bowtie2_mapping_mtnucratiocalculator", default="MT"
        ),
    wrapper:
        "v5.3.0/bio/mtnucratio"
