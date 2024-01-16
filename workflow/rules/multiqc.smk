rule mapping_multiqc_report:
    input:
        unpack(get_multiqc_report_input),
    output:
        report(
            "results/QC/MultiQC_Mapping.html",
            caption="../report/multiqc.rst",
            category="Quality Controls",
            subcategory="General",
            labels={
                "report": "html",
                "step": "Mapping",
            },
        ),
        "results/QC/MultiQC_Mapping_data.zip",
    threads: 1
    resources:
        # Reserve 2Gb per attempt
        mem_mb=lambda wildcards, attempt: (2 * 1024) * attempt,
        # Reserve 30min per attempt
        runtime=lambda wildcards, attempt: int(60 * 0.5) * attempt,
        tmpdir="tmp",
    params:
        extra="--zip-data-dir",
        use_input_files_only=True,
    log:
        "logs/multiqc.log",
    benchmark:
        "benchmark/multiqc.tsv"
    wrapper:
        "v3.3.3/bio/multiqc"
