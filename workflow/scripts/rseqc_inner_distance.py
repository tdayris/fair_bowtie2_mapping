# coding: utf-8

__author__ = "Thibault Dayris"
__mail__ = "thibault.dayris@gustaveroussy.fr"
__copyright__ = "Copyright 2024, Thibault Dayris"
__license__ = "MIT"

from snakemake import shell
from tempfile import TemporaryDirectory

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

with TemporaryDirectory() as tempdir:
    shell(
        "inner_distance.py "
        "--input-file {snakemake.input.aln} "
        "--out-prefix {tempdir}/out "
        "--refgene {snakemake.input.refgene} "
        "{log}"
    )

    if "txt" in snakemake.output.keys():
        shell(
            "mv --verbose {tempdir}/out.inner_distance.txt {snakemake.output.txt} {log}"
        )

    if "freq" in snakemake.output.keys():
        shell(
            "mv --verbose {tempdir}/out.inner_distance_freq.txt {snakemake.output.freq} {log}"
        )

    if "plot_r" in snakemake.output.keys():
        shell(
            "mv --verbose {tempdir}/out.inner_distance_plot.r {snakemake.output.plot_r} {log}"
        )
