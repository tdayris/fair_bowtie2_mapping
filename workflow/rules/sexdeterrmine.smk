rule fair_bowtie2_mapping_samtools_depth:
    input:
        bam=lambda wildcards: get_all_bams_per_genotype(
            wildcards,
            samples,
        ),
        bai=lambda wildcards: get_all_bams_per_genotype(
            wildcards,
            samples,
            index=True,
        ),
        bed=lambda wildcards: get_capturekit_bed(
            wildcards,
            config,
        ),
    output:
        temp(
            "tmp/fair_bowtie2_mapping_samtools_depth/{species}.{build}.{release}.{datatype}.depth.tsv"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1_000,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping_samtools_depth/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_samtools_depth/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_bowtie2_mapping_samtools_depth",
            default="--min-MQ 30 -l 10",
        ),
    wrapper:
        "v5.3.0/bio/samtools/depth"


rule fair_bowtie2_mapping_csvkit_add_header:
    input:
        table="tmp/fair_bowtie2_mapping_samtools_depth/{species}.{build}.{release}.{datatype}.depth.tsv",
    output:
        temp(
            "tmp/fair_bowtie2_mapping_csvkit_add_header/{species}.{build}.{release}.{datatype}.depth.tsv"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1_000,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping_csvkit_add_header/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_csvkit_add_header/{species}.{build}.{release}.{datatype}.benchmark"
    params:
        subcommand="add-header",
        extra=lambda wildcards: "-t -n Chr,Pos,"
        + ",".join(
            get_all_bams_per_genotype(
                wildcards,
                samples,
                samples_only=True,
            ),
        ),
    wrapper:
        "v5.3.0/utils/csvtk"


rule fair_bowtie2_mapping_sexdeterrmine:
    input:
        depth="tmp/fair_bowtie2_mapping_csvkit_add_header/{species}.{build}.{release}.{datatype}.depth.tsv",
    output:
        tsv=report(
            "results/{species}.{build}.{release}.{datatype}/QC/Sex.DetERRmine.tsv",
            caption="../report/sexdeterrmine.rst",
            category="Quality Controls",
            subcategory="General",
            labels={
                "report": "tsv",
                "step": "Mapping",
                "organism": "{species}.{build}.{release}.{datatype}",
            },
        ),
        json=temp(
            "tmp/fair_bowtie2_mapping_sexdeterrmine/{species}.{build}.{release}.{datatype}/sexdeterrmine.json",
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1_000,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping_sexdeterrmine/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_sexdeterrmine/{species}.{build}.{release}.{datatype}.benchmark"
    wrapper:
        "v5.3.0/bio/sexdeterrmine"
