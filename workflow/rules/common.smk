import csv
import os
import pandas
import snakemake
import snakemake.utils

from typing import Any, Callable, NamedTuple

snakemake_min_version: str = "8.13.0"
snakemake.utils.min_version(snakemake_min_version)

snakemake_docker_image: str = "docker://snakemake/snakemake:v8.13.0"


container: snakemake_docker_image


# Load and check configuration file
default_config_file: str = "config/config.yaml"


configfile: default_config_file


snakemake.utils.validate(config, "../schemas/config.schema.yaml")


# Load and check samples properties table
def load_table(path: str) -> pandas.DataFrame:
    """
    Load a table in memory, automatically inferring column separators

    Parameters:
    path (str): Path to the table to be loaded

    Return
    (pandas.DataFrame): The loaded table
    """
    with open(path, "r") as table_stream:
        dialect: csv.Dialect = csv.Sniffer().sniff(table_stream.readline())
        table_stream.seek(0)

    # Load table
    table: pandas.DataFrame = pandas.read_csv(
        path,
        sep=dialect.delimiter,
        header=0,
        index_col=None,
        comment="#",
        dtype=str,
    )

    # Remove empty lines
    table = table.where(table.notnull(), None)

    return table


def load_genomes(
    path: str | None = None, samples: pandas.DataFrame | None = None
) -> pandas.DataFrame:
    """
    Load genome file, build it if genome file is missing and samples is not None.

    Parameters:
    path    (str)               : Path to genome file
    samples (pandas.DataFrame)  : Loaded samples
    """
    if path is not None:
        genomes: pandas.DataFrame = load_table(path)

        if samples is not None:
            genomes = used_genomes(genomes, samples)
        return genomes

    elif samples is not None:
        return samples[["species", "build", "release"]].drop_duplicates(
            ignore_index=True
        )

    raise ValueError(
        "Provide either a path to a genome file, or a loaded samples table"
    )


def used_genomes(
    genomes: pandas.DataFrame, samples: pandas.DataFrame | None = None
) -> tuple[str]:
    """
    Reduce the number of genomes to download to the strict minimum
    """
    if samples is None:
        return genomes

    return genomes.loc[
        genomes.species.isin(samples.species.tolist())
        & genomes.build.isin(samples.build.tolist())
        & genomes.release.isin(samples.release.tolist())
    ]


# Load and check samples properties tables
try:
    if (samples is None) or samples.empty():
        sample_table_path: str = config.get("samples", "config/samples.csv")
        samples: pandas.DataFrame = load_table(sample_table_path)
except NameError:
    sample_table_path: str = config.get("samples", "config/samples.csv")
    samples: pandas.DataFrame = load_table(sample_table_path)

snakemake.utils.validate(samples, "../schemas/samples.schema.yaml")


# Load and check genomes properties table
genomes_table_path: str = config.get("genomes", "config/genomes.csv")
try:
    if (genomes is None) or genomes.empty:
        genomes: pandas.DataFrame = load_genomes(genomes_table_path, samples)
except NameError:
    genomes: pandas.DataFrame = load_genomes(genomes_table_path, samples)

snakemake.utils.validate(genomes, "../schemas/genomes.schema.yaml")


report: "../report/workflows.rst"


snakemake_wrappers_prefix: str = "v3.12.0"
release_tuple: tuple[str] = tuple(set(genomes.release.tolist()))
build_tuple: tuple[str] = tuple(set(genomes.build.tolist()))
species_tuple: tuple[str] = tuple(set(genomes.species.tolist()))
datatype_tuple: tuple[str] = ("dna", "cdna", "all", "transcripts")
gxf_tuple: tuple[str] = ("gtf", "gff3")
id2name_tuple: tuple[str] = ("t2g", "id_to_gene")
tmp: str = f"{os.getcwd()}/tmp"
samples_id_tuple: tuple[str] = tuple(samples.sample_id)
stream_tuple: tuple[str] = ("1", "2")


wildcard_constraints:
    sample=r"|".join(samples_id_tuple),
    release=r"|".join(release_tuple),
    build=r"|".join(build_tuple),
    species=r"|".join(species_tuple),
    datatype=r"|".join(datatype_tuple),
    stream=r"|".join(stream_tuple),


def lookup_config(
    dpath: str, default: str | None = None, config: dict[str, Any] = config
) -> str:
    """
    Run lookup function with default parameters in order to search a key in configuration and return a default value
    """
    value: str | None = default

    try:
        value = lookup(dpath=dpath, within=config)
    except LookupError:
        value = default
    except WorkflowError:
        value = default

    return value


def lookup_genomes(
    wildcards: snakemake.io.Wildcards,
    key: str,
    default: str | list[str] | None = None,
    genomes: pandas.DataFrame = genomes,
) -> str:
    """
    Run lookup function with default parameters in order to search user-provided sequence/annotation files
    """
    query: str = "species == '{wildcards.species}' & build == '{wildcards.build}' & release == '{wildcards.release}'".format(
        wildcards=wildcards
    )

    query_result: str | float = getattr(
        lookup(query=query, within=genomes), key, default
    )
    if (query_result != query_result) or (query_result is None):
        # Then the result of the query is nan
        return default
    return query_result


def get_dna_fasta(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final DNA fasta sequences
    """
    default: str = "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.dna.fasta".format(
        wildcards=wildcards
    )
    return lookup_genomes(wildcards, key="dna_fasta", default=default, genomes=genomes)


def get_cdna_fasta(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final cDNA fasta sequences
    """
    default: str = "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.cdna.fasta".format(
        wildcards=wildcards
    )
    return lookup_genomes(wildcards, key="cdna_fasta", default=default, genomes=genomes)


def get_transcripts_fasta(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final cDNA transcripts fasta sequences
    """
    default: str = "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.transcripts.fasta".format(
        wildcards=wildcards
    )
    return lookup_genomes(
        wildcards, key="transcripts_fasta", default=default, genomes=genomes
    )


def select_fasta(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Evaluates the {datatype} wildcard, and return the right fasta file
    """
    return branch(
        condition=str(wildcards.datatype).lower(),
        cases={
            "dna": get_dna_fasta(wildcards),
            "cdna": get_cdna_fasta(wildcards),
            "transcripts": get_transcripts_fasta(wildcards),
        },
    )


def get_dna_fai(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final DNA fasta sequences index
    """
    default: str = "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.dna.fasta.fai".format(
        wildcards=wildcards
    )
    return lookup_genomes(wildcards, key="dna_fai", default=default, genomes=genomes)


def get_cdna_fai(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final cDNA fasta sequences index
    """
    default: str = "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.cdna.fasta.fai".format(
        wildcards=wildcards
    )
    return lookup_genomes(wildcards, key="cdna_fai", default=default, genomes=genomes)


def get_transcripts_fai(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final cDNA transcripts fasta sequences index
    """
    default: str = "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.transcripts.fasta.fai".format(
        wildcards=wildcards
    )
    return lookup_genomes(
        wildcards, key="transcripts_fai", default=default, genomes=genomes
    )


def select_fai(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Evaluates the {datatype} wildcard, and return the right fasta index file
    """
    return branch(
        condition=str(wildcards.datatype).lower(),
        cases={
            "dna": get_dna_fai(wildcards),
            "cdna": get_cdna_fai(wildcards),
            "transcripts": get_transcripts_fai(wildcards),
        },
    )


def get_gtf(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final genome annotation (GTF formatted)
    """
    default: str = "reference/annotation/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.gtf".format(
        wildcards=wildcards
    )
    return lookup_genomes(wildcards, key="gtf", default=default, genomes=genomes)


def get_gff(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final genome annotation (GFF3 formatted)
    """
    default: str = "reference/annotation/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.gff3".format(
        wildcards=wildcards
    )
    return lookup_genomes(wildcards, key="gff3", default=default, genomes=genomes)


def get_dna_bowtie2_index(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final Bowtie2 index DNA index
    """
    default: list[str] = multiext(
        "reference/bowtie2_index/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.dna".format(
            wildcards=wildcards
        ),
        ".1.bt2",
        ".2.bt2",
        ".3.bt2",
        ".4.bt2",
        ".rev.1.bt2",
        ".rev.2.bt2",
    )
    return lookup_genomes(
        wildcards, key="bowtie2_dna_index", default=default, genomes=genomes
    )


def get_cdna_bowtie2_index(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final Bowtie2 index cDNA index
    """
    default: list[str] = multiext(
        "reference/bowtie2_index/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.cdna".format(
            wildcards=wildcards
        ),
        ".1.bt2",
        ".2.bt2",
        ".3.bt2",
        ".4.bt2",
        ".rev.1.bt2",
        ".rev.2.bt2",
    )
    return lookup_genomes(
        wildcards, key="bowtie2_cdna_index", default=default, genomes=genomes
    )


def get_transcripts_bowtie2_index(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final Bowtie2 index trnascripts index
    """
    default: list[str] = multiext(
        "reference/bowtie2_index/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.transcripts".format(
            wildcards=wildcards
        ),
        ".1.bt2",
        ".2.bt2",
        ".3.bt2",
        ".4.bt2",
        ".rev.1.bt2",
        ".rev.2.bt2",
    )
    return lookup_genomes(
        wildcards, key="bowtie2_transcripts_index", default=default, genomes=genomes
    )


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
    # Include requirements of previous pipelines
    fair_genome_indexer_targets: dict[
        str, list[str] | str
    ] = fair_genome_indexer.get_fair_genome_indexer_target(
        wildcards=wildcards, genomes=genomes, samples=samples
    )

    fair_fastqc_multiqc_targets: dict[str, list[str]] = dict(
        fair_fastqc_multiqc.rules.fair_fastqc_multiqc_target.input
    )

    # Build dict of expected outputs
    usable_genomes: list[NameTyple] | NamedTuple = lookup(
        query="species != '' & build != '' & release != ''",
        within=used_genomes(genomes, samples),
    )
    if not isinstance(usable_genomes, list):
        usable_genomes = [usable_genomes]

    genome_properties: list[str] = [
        ".".join([usable_genome.species, usable_genome.build, usable_genome.release])
        for usable_genome in usable_genomes
    ]
    results: dict[str, list[str] | str] = (
        fair_fastqc_multiqc_targets | fair_genome_indexer_targets
    )
    results["multiqc_mapping"] = (
        f"results/{genome_property}.dna/QC/MultiQC_Mapping.html"
        for genome_property in genome_properties
    )

    results["bams"] = []
    results["bais"] = []

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
