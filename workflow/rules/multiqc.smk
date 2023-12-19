rule multiqc_report:
    input:
        unpack(get_multiqc_report_input),
    output:
        report(
            "results/QC/MultiQC.html",
            caption="../report/multiqc.rst",
            category="Quality Controls",
            subcategory="General",
            labels={
                "report": "html",
            },
        ),
        "results/QC/MultiQC_data.zip",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        runtime=lambda wildcards, attempt: int(60 * 0.75) * attempt,
        tmpdir="tmp",
    params:
        extra="--zip-data-dir",
        use_input_files_only=True,
    log:
        "logs/multiqc.log",
    benchmark:
        "benchmark/multiqc.tsv"
    wrapper:
        "v3.2.0/bio/multiqc"
