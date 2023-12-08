rule fastp_trimming_pair_ended:
    input:
        unpack(get_fastp_trimming_input),
    output:
        trimmed=temp(
            expand(
                "tmp/fastp/trimmed/{sample}.{stream}.fastq",
                stream=["1", "2"],
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
        json=temp("tmp/fastp/report_pe/{sample}.fastp.json"),
    log:
        "logs/fastp/{sample}.log",
    benchmark:
        "benchmark/fastp/{sample}.tsv"
    params:
        adapters=config.get("params", {}).get("fastp", {}).get("adapters", ""),
        extra=config.get("params", {}).get("fastp", {}).get("extra", ""),
    wrapper:
        f"{snakemake_wrappers_version}/bio/fastp"


use rule fastp_trimming_pair_ended as fastp_trimming_single_ended with:
    output:
        trimmed=temp("tmp/fastp/trimmed/{sample}.fastq"),
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
        json=temp("tmp/fastp/report_se/{sample}.fastp.json"),
