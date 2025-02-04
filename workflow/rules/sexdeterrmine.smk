rule fair_bowtie2_mapping_samtools_depth:
    input:
        bams=lambda wildcards: get_all_bams_per_genotype(
            wildcards,
            samples,
        ),
        bais=lambda wildcards: get_all_bams_per_genotype(
            wildcards,
            samples,
            index=True,
        ),
        bed=lambda wildcards: get_intervals(
            wildcards,
            genomes,
        ),
    output:
        temp(
            "tmp/fair_bowtie2_mapping_samtools_depth/{species}.{build}.{release}.{datatype}.depth.tsv"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 500,
        runtime=lambda wildcards, attempt: attempt * 60 * 2,
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
        "v5.6.0/bio/samtools/depth"


rule fair_bowtie2_mapping_sample_list:
    output:
        "tmp/fair_bowtie2_mapping_sample_list/{species}.{build}.{release}.{datatype}.txt",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 200 + 300,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping_sample_list/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_sample_list/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=lambda wildcards: "\n".join(
            get_all_bams_per_genotype(
                wildcards,
                samples,
                samples_only=True,
            )
        ),
    conda:
        "../envs/bash.yaml"
    shell:
        "echo '{params.extra}' > {output} 2> {log}"


rule fair_bowtie2_mapping_sexdeterrmine:
    input:
        depth="tmp/fair_bowtie2_mapping_samtools_depth/{species}.{build}.{release}.{datatype}.depth.tsv",
        samples="tmp/fair_bowtie2_mapping_sample_list/{species}.{build}.{release}.{datatype}.txt",
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
    params:
        extra=lambda wildcards, input: f" --SampleList {input.samples} ",
    wrapper:
        "v5.6.0/bio/sexdeterrmine"
