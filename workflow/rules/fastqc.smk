rule fastqc_pair_ended:
    input:
        unpack(get_fastqc_input),
    output:
        html=report(
            "results/QC/report_pe/{sample}.{stream}.html",
            caption="../report/fastqc.rst",
            category="Quality Controls",
            subcategory="Raw",
            labels={
                "report": "html",
                "sample": "{sample}",
                "library": "pair_ended",
            },
        ),
        zip="results/QC/report_pe/{sample}.{stream}_fastqc.zip",
    threads: 1
    resources:
        # Reserve 2Gb per attempt
        mem_mb=lambda wildcards, attempt: (2 * 1024) * attempt,
        # Reserve 30min per attempt
        runtime=lambda wildcards, attempt: int(60 * 0.5) * attempt,
        tmpdir="tmp",
    log:
        "logs/fastqc/{sample}.{stream}.log",
    benchmark:
        "benchmark/fastqc/{sample}.{stream}.tsv"
    params:
        extra=config.get("params", {}).get("fastqc", ""),
    wrapper:
        "v3.3.3/bio/fastqc"


use rule fastqc_pair_ended as fastqc_single_ended with:
    output:
        html=report(
            "results/QC/report_pe/{sample}.html",
            caption="../report/fastqc.rst",
            category="Quality Controls",
            subcategory="Raw",
            labels={
                "report": "html",
                "sample": "{sample}",
                "library": "single_ended",
            },
        ),
        zip="results/QC/report_pe/{sample}_fastqc.zip",
    log:
        "logs/fastqc/{sample}.log",
    benchmark:
        "benchmark/fastqc/{sample}.tsv"
