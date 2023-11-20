rule multiqc_report:
    input:
        unpack(get_multiqc_report_input),
    output:
        "results/QC/MultiQC.html",
        "results/QC/MultiQC_data.zip",
    params:
        extra="--zip-data-dir",
    log:
        "logs/multiqc.log",
    benchmark:
        "benchmark/multiqc.tsv",
    wrapper:
        f"{snakemake_wrappers_version}/bio/multiqc"