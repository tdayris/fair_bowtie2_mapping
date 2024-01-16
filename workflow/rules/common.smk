import csv
import pandas
import snakemake
import snakemake.utils

from collections import defaultdict
from pathlib import Path
from typing import Any

snakemake.utils.min_version("7.29.0")

# containerized: "docker://snakemake/snakemake:v7.32.4"
# containerized: "docker://mambaorg/micromamba:git-8440cec-jammy-cuda-12.2.0"
# containerized: "docker://condaforge/mambaforge:23.3.1-1"


# Load and check configuration file
configfile: "config/config.yaml"


snakemake.utils.validate(config, "../schemas/config.schema.yaml")

# Load and check samples properties table
sample_table_path: str = config.get("samples", "config/samples.csv")
with open(sample_table_path, "r") as sample_table_stream:
    dialect: csv.Dialect = csv.Sniffer().sniff(sample_table_stream.readline())
    sample_table_stream.seek(0)

samples: pandas.DataFrame = pandas.read_csv(
    filepath_or_buffer=sample_table_path,
    sep=dialect.delimiter,
    header=0,
    index_col=None,
    comment="#",
    dtype=str,
)
samples = samples.where(samples.notnull(), None)
snakemake.utils.validate(samples, "../schemas/samples.schema.yaml")

# This is here for compatibility with
genome_table_path: str = config.get("genomes")
if genome_table_path:
    with open(genome_table_path, "r") as genome_table_stream:
        dialect: csv.Dialect = csv.Sniffer().sniff(genome_table_stream.readline())
        genome_table_stream.seek(0)

    genomes: pandas.DataFrame = pandas.read_csv(
        filepath_or_buffer=genome_table_path,
        sep=dialect.delimiter,
        header=0,
        index_col=None,
        comment="#",
        dtype=str,
    )
    genomes = genomes.where(genomes.notnull(), None)
else:
    genomes: pandas.DataFrame = samples[
        ["species", "build", "release"]
    ].drop_duplicates(keep="first", ignore_index=True)
    genomes.to_csv("genomes.csv", sep=",", index=False, header=True)
    config["genomes"] = "genomes.csv"

snakemake.utils.validate(genomes, "../schemas/genomes.schema.yaml")

snakemake_wrappers_version: str = "v3.0.0"


report: "../report/workflows.rst"


release_list: list[str] = list(set(genomes.release.tolist()))
build_list: list[str] = list(set(genomes.build.tolist()))
species_list: list[str] = list(set(genomes.species.tolist()))
stream_list: list[str] = ["1", "2"]


wildcard_constraints:
    sample=r"|".join(samples.sample_id),
    release=r"|".join(release_list),
    build=r"|".join(build_list),
    species=r"|".join(species_list),
    stream=r"|".join(stream_list),


def get_reference_genome_data(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame
) -> dict[str, str | None]:
    """
    Return genome information for a given set of {species, build, release} wildcards

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    genomes   (pandas.DataFrame)      : Describe genomes and reference file(s)

    Return (dict[str, str | None]):
    Genome information
    """
    result: str | None = genomes.loc[
        (genomes["species"] == str(wildcards.species))
        & (genomes["build"] == str(wildcards.build))
        & (genomes["release"] == str(wildcards.release))
    ]
    if len(result) > 0:
        return next(iter(result.to_dict(orient="index").values()))
    return defaultdict(lambda: None)


def get_sample_information(
    wildcards: snakemake.io.Wildcards, samples: pandas.DataFrame
) -> dict[str, str | None]:
    """
    Return sample information for a given {sample} wildcards

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe samples and their input files

    Return (dict[str, str | None]):
    Sample information
    """
    result: str | None = samples.loc[(samples["sample_id"] == str(wildcards.sample))]
    if len(result) > 0:
        return next(iter(result.to_dict(orient="index").values()))
    return defaultdict(lambda: None)


def get_fastp_trimming_input(
    wildcards: snakemake.io.Wildcards,
    samples: pandas.DataFrame = samples,
    config: dict[str, Any] = config,
) -> dict[str, list[str]]:
    """
    Return expected input files for Bowtie2 mapping, according to user-input,
    and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome
    config    (dict[str, Any])        : Configuration file

    Return (dict[str, list[str]]):
    Dictionnary of all input files as required by Fastp's snakemake-wrapper
    """
    sample_data: dict[str, str | None] = get_sample_information(wildcards, samples)
    downstream_file: str | None = sample_data.get("downstream_file")
    if downstream_file or not pandas.isna(downstream_file):
        return {
            "sample": [sample_data["upstream_file"], downstream_file],
        }

    return {
        "sample": [sample_data["upstream_file"]],
    }


def get_fastqc_input(
    wildcards: snakemake.io.Wildcards, samples: pandas.DataFrame = samples
) -> str:
    """
    Return expected input files for FastQC, according to user-input,
    and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe samples and their genome

    Return (str):
    Path to a fastq file, as required by FastQC's snakemake-wrapper
    """
    sample_data: dict[str, str | None] = get_sample_information(wildcards, samples)
    downstream_file: str | None = sample_data.get("downstream_file")
    if "stream" in wildcards.keys():
        if wildcards.stream == "1":
            return {"fastq": sample_data["upstream_file"]}
        elif wildcards.stream == "2" and downstream_file:
            return {"fastq": sample_data["downstream_file"]}
        raise ValueError("Could not guess which fastq we're talking about")
    return {"fastq": sample_data["upstream_file"]}


def get_bowtie2_build_input(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
):
    """
    Return expected input files for Bowtie2 indexing, according to user-input,
    and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    genomes   (pandas.DataFrame)      : Describe genomes and reference file(s)

    Return (dict[str, str]):
    Dictionnary of all input files as required by Bowtie2's snakemake-wrapper
    """
    species: str = str(wildcards["species"])
    build: str = str(wildcards["build"])
    release: str = str(wildcards["release"])
    datatype: str = "dna"

    idx: str | None = get_reference_genome_data(wildcards, genomes).get("fasta")
    if idx:
        return {"ref": idx}
    else:
        return {"ref": f"reference/{species}.{build}.{release}.{datatype}.fasta"}


def get_bowtie2_alignment_input(
    wildcards: snakemake.io.Wildcards,
    samples: pandas.DataFrame = samples,
    config: dict[str, Any] = config,
    genomes: pandas.DataFrame = genomes,
) -> dict[str, list[str] | str]:
    """
    Return expected input files for Bowtie2 mapping, according to user-input,
    and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome
    config    (dict[str, Any])        : Configuration file

    Return ( dict[str, list[str] | str]):
    Dictionnary of all input files as required by Bowtie2's snakemake-wrapper
    """
    sample_data: dict[str, str | None] = get_sample_information(wildcards, samples)
    if not sample_data:
        raise KeyError(
            f"Could not find sample {str(wildcards.sample)} in the sample.csv file."
        )

    species: str = str(sample_data["species"])
    build: str = str(sample_data["build"])
    release: str = str(sample_data["release"])
    datatype: str = "dna"

    idx: str | None = get_reference_genome_data(wildcards, genomes).get("bowtie2_index")
    if idx:
        idx = [str(file) for file in Path(idx) if str(file).endswith(".bt2")]
    else:
        idx = multiext(
            f"reference/{species}.{build}.{release}.{datatype}",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        )

    results: dict[str, list[str]] = {
        "idx": idx,
        "sample": [],
    }
    downstream_file: str | None = sample_data.get("downstream_file")
    if downstream_file:
        results["sample"] = expand(
            "tmp/fastp/trimmed/{sample}.{stream}.fastq",
            stream=["1", "2"],
            sample=[str(wildcards.sample)],
        )
    else:
        results["sample"] = ["tmp/fastp/trimmed/{sample}.fastq"]

    return results


def get_picard_create_multiple_metrics_input(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> dict[str, str]:
    """
    Return expected input files for Picard CreateMultipleMetrics, according to user-input,
    and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome
    genomes   (pandas.DataFrame)      : Describe genome

    Return (dict[str, str]):
    Dictionnary of all input files as required by Picard's snakemake-wrapper
    """
    species: str = str(wildcards.species)
    build: str = str(wildcards.build)
    release: str = str(wildcards.release)
    datatype: str = "dna"
    sample: str = str(wildcards.sample)
    idx: str = get_reference_genome_data(wildcards, genomes).get(
        "fasta", f"reference/{species}.{build}.{release}.{datatype}.fasta"
    )
    return {
        "bam": f"results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
        "bai": f"results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai",
        "ref": idx,
    }


def get_multiqc_report_input(
    wildcards: snakemake.io.Wildcards, samples: pandas.DataFrame = samples
) -> dict[str, list[str]]:
    """
    Return expected input files for MultiQC report, according to user-input,
    and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome

    Return (dict[str, list[str]]):
    Dictionnary of all input files as required by MultiQC's snakemake-wrapper
    """
    results: dict[str, list[str]] = {
        "picard_qc": [],
        "fastp": [],
        "fastqc": [],
        "bowtie2": [],
        "samtools": [],
    }
    datatype: str = "dna"
    sample_iterator = zip(
        samples.sample_id,
        samples.species,
        samples.build,
        samples.release,
    )
    for sample, species, build, release in sample_iterator:
        results["picard_qc"] += multiext(
            f"tmp/picard/{species}.{build}.{release}.{datatype}/stats/{sample}",
            ".alignment_summary_metrics",
            ".insert_size_metrics",
            ".insert_size_histogram.pdf",
            ".base_distribution_by_cycle_metrics",
            ".base_distribution_by_cycle.pdf",
            ".gc_bias.detail_metrics",
            ".gc_bias.summary_metrics",
            ".gc_bias.pdf",
        )
        sample_data: dict[str, str | None] = get_sample_information(
            snakemake.io.Wildcards(fromdict={"sample": sample}), samples
        )
        if sample_data.get("downstream_file"):
            results["fastp"].append(f"tmp/fastp/report_pe/{sample}.fastp.json")
            results["fastqc"].append(f"results/QC/report_pe/{sample}.1_fastqc.zip")
            results["fastqc"].append(f"results/QC/report_pe/{sample}.2_fastqc.zip")
        else:
            results["fastp"].append(f"tmp/fastp/report_se/{sample}.fastp.json")
            results["fastqc"].append(f"results/QC/report_pe/{sample}_fastqc.zip")

        results["bowtie2"].append(
            f"logs/bowtie2/align/{species}.{build}.{release}.{datatype}/{sample}.log"
        )

        results["samtools"].append(
            f"tmp/samtools/{species}.{build}.{release}.{datatype}/{sample}.txt"
        )

    return results


def get_fair_bowtie2_mapping_target(
    wildcards: snakemake.io.Wildcards,
    samples: pandas.DataFrame = samples,
    config: dict[str, Any] = config,
) -> dict[str, list[str]]:
    """
    Return the expected list of output files at the end of the pipeline

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome
    config    (dict[str, Any])        : Configuration file

    Return (dict[str, List(str)]):
    Dictionnary of expected output files
    """
    results: dict[str, list[str]] = {
        "multiqc": [
            "results/QC/MultiQC_FastQC.html",
            "results/QC/MultiQC_FastQC_data.zip",
            "results/QC/MultiQC_Mapping.html",
            "results/QC/MultiQC_Mapping_data.zip",
        ],
        "bams": [],
        "bais": [],
    }
    sample_iterator = zip(
        samples.sample_id,
        samples.species,
        samples.build,
        samples.release,
    )
    for sample, species, build, release in sample_iterator:
        results["bams"].append(
            f"results/{species}.{build}.{release}.dna/Mapping/{sample}.bam"
        )
        results["bais"].append(
            f"results/{species}.{build}.{release}.dna/Mapping/{sample}.bam.bai"
        )

    return results
