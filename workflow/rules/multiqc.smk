rule fair_bowtie2_mapping_multiqc_config:
    input:
        "tmp/fair_fastqc_multiqc/bigr_logo.png",
    output:
        temp("tmp/fair_bowtie2_mapping/multiqc_config.yaml"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512,
        runtime=lambda wildcards, attempt: attempt * 5,
        tmpdir=tmp,
    localrule: True
    log:
        "logs/fair_bowtie2_mapping/multiqc_config.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping/multiqc_config.tsv"
    params:
        extra=lambda wildcards, input: {
            "title": "Mapping quality control report",
            "subtitle": "Produced on raw fastq recieved from sequencer",
            "intro_text": (
                "This pipeline building this report has "
                "no information about sequencing protocol, "
                "wet-lab experimental design, nor sample organisms."
            ),
            "report_comment": (
                "This report was generated using: "
                "https://github.com/tdayris/fair_bowtie2_mapping"
            ),
            "show_analysis_paths": False,
            "show_analysis_time": False,
            "custom_logo": input[0],
            "custom_logo_url": "https://bioinfo_gustaveroussy.gitlab.io/bigr/webpage/",
            "custom_logo_title": "Bioinformatics Platform @ Gustave Roussy",
            "report_header_info": [
                {"Contact E-mail": "bigr@gustaveroussy.fr"},
                {"Application type": "Short-gapped reads"},
                {"Project Type": "Mapping"},
            ],
            "software_versions": {
                "Quality controls": {
                    "fastqc": "1.12.1",
                    "fastq_screen": "0.15.3",
                    "bowtie2": "1.3.1",
                    "multiqc": "1.20.0",
                },
                "Mapping": {
                    "bowtie2": "2.5.3",
                    "sambamba": "1.0",
                    "samtools": "1.19.2",
                    "picard": "3.1.1",
                    "rseqc": "5.0.3",
                    "fastp": "0.23.4",
                    "ngsderive": "3.3.2",
                    "goleft": "0.2.4",
                },
                "Pipeline": {
                    "snakemake": "8.5.3",
                    "fair_fastqc_multiqc": "2.1.2",
                    "fair_genome_indexer": "3.2.2",
                },
            },
            "disable_version_detection": True,
            "run_modules": [
                "fastqc",
                "fastq_screen",
                "fastp",
                "bowtie2",
                "samtools",
                "picard",
                "rseqc",
                "ngsderive",
                "goleft_indexcov",
            ],
            "report_section_order": {
                "fastq_screen": {"order": 1000},
                "ngsderive": {"order": 950},
                "fastqc": {"order": 900},
                "fastp": {"order": 890},
                "bowtie2": {"order": 880},
                "picard": {"order": 870},
                "samtools": {"order": 860},
                "rseqc": {"order": 850},
                "goleft_indexcov": {"order": 840},
                "software_versions": {"order": -1000},
            },
        },
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/fair_bowtie2_mapping_multiqc_config.py"


rule fair_bowtie2_mapping_multiqc_report:
    input:
        unpack(get_multiqc_report_input),
    output:
        report(
            "results/{species}.{build}.{release}.{datatype}/QC/MultiQC_Mapping.html",
            caption="../report/multiqc.rst",
            category="Quality Controls",
            subcategory="General",
            labels={
                "report": "html",
                "step": "Mapping",
                "organism": "{species}.{build}.{release}.{datatype}",
            },
        ),
        "results/{species}.{build}.{release}.{datatype}/QC/MultiQC_Mapping_data.zip",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (2 * 1024) * attempt,
        runtime=lambda wildcards, attempt: int(60 * 0.5) * attempt,
        tmpdir=tmp,
    params:
        extra=lookup(dpath="params/fair_bowtie2_mapping/multiqc", within=config),
        use_input_files_only=True,
    log:
        "logs/fair_bowtie2_mapping/multiqc_report/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping/multiqc_report/{species}.{build}.{release}.{datatype}.tsv"
    wrapper:
        "v3.5.0/bio/multiqc"
