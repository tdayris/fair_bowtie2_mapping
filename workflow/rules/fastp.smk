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
        html=temp("tmp/fastp/report_pe/{sample}.html"),
        json=temp("tmp/fastp/report_pe/{sample}.json"),
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
        html=temp("tmp/fastp/report_se/{sample}.html"),
        json=temp("tmp/fastp/report_se/{sample}.json"),
