import csv
import os
import pandas
import snakemake
import snakemake.utils

from collections import defaultdict
from pathlib import Path
from typing import Any, NamedTuple

snakemake.utils.min_version("8.1.0")


container: "docker://snakemake/snakemake:v8.5.3"


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


report: "../report/workflows.rst"


snakemake_wrappers_prefix: str = "v3.5.2"
release_list: list[str] = list(set(genomes.release.tolist()))
build_list: list[str] = list(set(genomes.build.tolist()))
species_list: list[str] = list(set(genomes.species.tolist()))
datatype_list: list[str] = ["dna", "cdna", "transcripts"]
stream_list: list[str] = ["1", "2"]
tmp: str = f"{os.getcwd()}/tmp"


wildcard_constraints:
    sample=r"|".join(samples.sample_id),
    release=r"|".join(release_list),
    build=r"|".join(build_list),
    species=r"|".join(species_list),
    datatype=r"|".join(datatype_list),
    stream=r"|".join(stream_list),


def dlookup(
    dpath: str | None = None,
    query: str | None = None,
    cols: list[str] | None = None,
    within=None,
    default: str | dict[str, Any] | None = None,
) -> str:
    """
    Return lookup() results or defaults

    dpath   (str | Callable | None): Passed to dpath library
    query   (str | Callable | None): Passed to DataFrame.query()
    cols    (list[str] | None):      The columns to operate on
    within  (object):                The dataframe or mappable object
    default (str):                   The default value to return
    """
    value = None
    try:
        value = lookup(dpath=dpath, query=query, cols=cols, within=within)
    except LookupError:
        value = default
    except WorkflowError:
        value = default
    except KeyError:
        value = default
    except AttributeError:
        value = default

    return value


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
    species: str = str(wildcards.species)
    build: str = str(wildcards.build)
    release: str = str(wildcards.release)
    datatype: str = str(wildcards.datatype)
    results: dict[str, list[str]] = {
        "config": "tmp/fair_bowtie2_mapping/multiqc_config.yaml",
        "logo": "tmp/fair_fastqc_multiqc/bigr_logo.png",
        "fastp_pair_ended": collect(
            "tmp/fair_bowtie2_mapping/fastp_trimming_pair_ended/{sample.sample_id}.fastp.json",
            sample=lookup(
                query=f"downstream_file == downstream_file & species == '{species}' & build == '{build}' & release == '{release}'",
                within=samples,
            ),
        ),
        "fastp_single_ended": collect(
            "tmp/fair_bowtie2_mapping/fastp_trimming_single_ended/{sample.sample_id}.fastp.json",
            sample=lookup(
                query=f"downstream_file != downstream_file & species == '{species}' & build == '{build}' & release == '{release}'",
                within=samples,
            ),
        ),
        "fastqc_pair_ended": collect(
            "results/QC/report_pe/{sample.sample_id}.{stream}_fastqc.zip",
            sample=lookup(
                query=f"downstream_file == downstream_file & species == '{species}' & build == '{build}' & release == '{release}'",
                within=samples,
            ),
            stream=stream_list,
        ),
        "fastqc_single_ended": collect(
            "results/QC/report_pe/{sample.sample_id}_fastqc.zip",
            sample=lookup(
                query=f"downstream_file != downstream_file & species == '{species}' & build == '{build}' & release == '{release}'",
                within=samples,
            ),
        ),
        "bowtie2": [],
        "samtools": [],
        "picard_qc": [],
        "ngsderive_readlen": [],
        "ngsderive_instrument": [],
        "ngsderive_encoding": [],
        "ngsderive_strandedness": [],
        "ngsderive_endedness": [],
        "goleft": [],
        "rseqc_infer_experiment": collect(
            "tmp/fair_bowtie2_mapping/rseqc_infer_experiment/{sample.species}.{sample.build}.{sample.release}.dna/{sample.sample_id}.infer_experiment.txt",
            sample=lookup(
                query=f"species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        "rseqc_bamstat": collect(
            "tmp/fair_bowtie2_mapping/rseqc_bamstat/{sample.species}.{sample.build}.{sample.release}.dna/{sample.sample_id}.bamstat.txt",
            sample=lookup(
                query=f"species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        "rseqc_read_gc": collect(
            "tmp/fair_bowtie2_mapping/rseqc_read_gc/{sample.species}.{sample.build}.{sample.release}.dna/{sample.sample_id}.GC.xls",
            sample=lookup(
                query=f"species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        "rseqc_read_distribution": collect(
            "tmp/fair_bowtie2_mapping/rseqc_read_distribution/{sample.species}.{sample.build}.{sample.release}.dna/{sample.sample_id}.txt",
            sample=lookup(
                query=f"species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        "rseqc_inner_distance": collect(
            "tmp/fair_bowtie2_mapping/rseqc_inner_distance/{sample.species}.{sample.build}.{sample.release}.dna/{sample.sample_id}.inner_distance_freq.txt",
            sample=lookup(
                query=f"species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
    }
    sample_iterator = zip(
        samples.sample_id,
        samples.species,
        samples.build,
        samples.release,
    )
    for sample, species, build, release in sample_iterator:
        results["bowtie2"].append(
            f"logs/fair_bowtie2_mapping/bowtie2_alignment/{species}.{build}.{release}.{datatype}/{sample}.log"
        )
        results["samtools"].append(
            f"tmp/fair_bowtie2_mapping/samtools_stats/{species}.{build}.{release}.{datatype}/{sample}.txt"
        )
        results["picard_qc"] += multiext(
            f"tmp/fair_bowtie2_mapping/picard_create_multiple_metrics/{species}.{build}.{release}.{datatype}/stats/{sample}",
            ".alignment_summary_metrics",
            ".insert_size_metrics",
            ".insert_size_histogram.pdf",
            ".base_distribution_by_cycle_metrics",
            ".base_distribution_by_cycle.pdf",
            ".gc_bias.detail_metrics",
            ".gc_bias.summary_metrics",
            ".gc_bias.pdf",
        )

        results["ngsderive_readlen"].append(
            f"tmp/fair_bowtie2_mapping/ngsderive/readlen/{species}.{build}/{release}.{datatype}/{sample}.readlen.tsv"
        )

        results["ngsderive_instrument"].append(
            f"tmp/fair_bowtie2_mapping/ngsderive/instrument/{species}.{build}/{release}.{datatype}/{sample}.instrument.tsv"
        )

        results["ngsderive_encoding"].append(
            f"tmp/fair_bowtie2_mapping/ngsderive/encoding/{species}.{build}/{release}.{datatype}/{sample}.encoding.tsv"
        )

        # results["ngsderive_strandedness"].append(
        #     f"tmp/fair_bowtie2_mapping/ngsderive/strandedness/{species}.{build}/{release}.{datatype}/{sample}.strandedness.tsv"
        # )

        results["goleft"].append(
            f"tmp/fair_bowtie2_mapping/goleft/indexcov/{species}.{release}.{build}/{sample}-indexcov.ped"
        )
        results["goleft"].append(
            f"tmp/fair_bowtie2_mapping/goleft/indexcov/{species}.{release}.{build}/{sample}-indexcov.roc"
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
        ],
        "multiqc_mapping": [],
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
        results["multiqc_mapping"].append(
            f"results/{species}.{build}.{release}.dna/QC/MultiQC_Mapping.html"
        )
        results["bams"].append(
            f"results/{species}.{build}.{release}.dna/Mapping/{sample}.bam"
        )
        results["bais"].append(
            f"results/{species}.{build}.{release}.dna/Mapping/{sample}.bam.bai"
        )

    results["multiqc_mapping"] = list(set(results["multiqc_mapping"]))
    return results
