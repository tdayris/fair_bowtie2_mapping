# -*- coding: utf-8 -*-

"""Snakemake wrapper for MultiQC configuration file"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2024, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import yaml

from typing import Any

default_config: dict[str, Any] = {
    "title": "Mapping quality control report",
    "subtitle": "Produced on raw fastq recieved from sequencer",
    "intro_text": (
        "This pipeline building this report has "
        "no information about sequencing protocol. "
    ),
    "report_comment": (
        "This report was generated using: "
        "https://github.com/tdayris/fair_bowtie2_mapping"
    ),
    "show_analysis_paths": False,
    "show_analysis_time": False,
    "custom_logo": snakemake.input[0],
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
            "bowtie2": "2.5.3",
            "bowtie1": "1.3.1",
            "multiqc": "1.21.0",
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
            "snakemake": "8.13.0",
            "snakemake-wrappers-utils": "0.6.2",
            "fair_fastqc_multiqc": "2.2.8",
            "fair_genome_indexer": "3.6.0",
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
}

config: dict[str, Any] | None = snakemake.params.get("extra", None)
if config is None:
    config = default_config.copy()


with open(str(snakemake.output[0]), "w") as out_yaml_stream:
    out_yaml_stream.write(yaml.dump(config, default_flow_style=False))
