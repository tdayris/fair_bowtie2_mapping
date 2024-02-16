rule fair_bowtie2_mapping_fastp_trimming_pair_ended:
    input:
        sample=expand(
            "tmp/fair_fastqc_multiqc/link_or_concat_pair_ended_input/{sample}.{stream}.fastq.gz",
            stream=stream_list,
            allow_missing=True,
        ),
    output:
        trimmed=temp(
            expand(
                "tmp/fair_bowtie2_mapping/fastp_trimming_pair_ended/{sample}.{stream}.fastq",
                stream=stream_list,
                allow_missing=True,
            )
        ),
        html=report(
            "results/QC/report_pe/{sample}.html",
            caption="../report/fastp.rst",
            category="Quality Controls",
            subcategory="Trimming",
            labels={
                "report": "html",
                "sample": "{sample}",
                "library": "pair_ended",
            },
        ),
        json=temp(
            "tmp/fair_bowtie2_mapping/fastp_trimming_pair_ended/{sample}.fastp.json"
        ),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: (4 * 1024) * attempt,
        runtime=lambda wildcards, attempt: int(60 * 0.5) * attempt,
        tmpdir="tmp",
    log:
        "logs/fair_bowtie2_mapping/fastp_trimming_pair_ended/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping/fastp_trimming_pair_ended/{sample}.tsv"
    params:
        adapters=lookup(dpath="params/fastp/adapters", within=config),
        extra=lookup(dpath="params/fastp/extra", within=config),
    wrapper:
        "v3.3.6/bio/fastp"


use rule fair_bowtie2_mapping_fastp_trimming_pair_ended as fair_bowtie2_mapping_fastp_trimming_single_ended with:
    input:
        sample=[
            "tmp/fair_fastqc_multiqc/link_or_concat_single_ended_input/{sample}.fastq.gz"
        ],
    output:
        trimmed=temp(
            "tmp/fair_bowtie2_mapping/fastp_trimming_single_ended/{sample}.fastq"
        ),
        html=report(
            "results/QC/report_se/{sample}.html",
            caption="../report/fastp.rst",
            category="Quality Controls",
            subcategory="Trimming",
            labels={
                "report": "html",
                "sample": "{sample}",
                "library": "single_ended",
            },
        ),
        json=temp(
            "tmp/fair_bowtie2_mapping/fastp_trimming_single_ended/{sample}.fastp.json"
        ),
    log:
        "logs/fair_bowtie2_mapping/fastp_trimming_single_ended/{sample}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping/fastp_trimming_single_ended/{sample}.tsv"
