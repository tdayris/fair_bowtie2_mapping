"""
## Memory
Requires a job with at most 255.86  Mb,
 on average 219.45 ± 96.33 Mb, 
on Gustave Roussy's HPC Flamingo, on a 1.0  Mb dataset.
## Time
A job took 0:00:05 to proceed,
on average 0:00:05 ± 0:00:01
"""


rule fair_bowtie2_mapping_multiqc_config:
    input:
        "tmp/fair_fastqc_multiqc_bigr_logo.png",
    output:
        temp("tmp/fair_bowtie2_mapping_multiqc_config.yaml"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 300 + attempt * 100,
        runtime=lambda wildcards, attempt: attempt * 5,
        tmpdir=tmp,
    localrule: True
    log:
        "logs/fair_bowtie2_mapping_config.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_config.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_bowtie2_mapping_multiqc_config", default=None
        ),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/fair_bowtie2_mapping_multiqc_config.py"


"""
## Memory
Requires a job with at most 6812.27  Mb,
 on average 5839.23 ± 2574.42 Mb, 
on Gustave Roussy's HPC Flamingo, on a 1.0  Mb dataset.
## Time
A job took 0:01:54 to proceed,
on average 0:01:38 ± 0:00:42
"""


rule fair_bowtie2_mapping_multiqc_report:
    input:
        config="tmp/fair_bowtie2_mapping_multiqc_config.yaml",
        picard_qc=collect(
            "tmp/fair_bowtie2_mapping_picard_create_multiple_metrics/{sample.species}.{sample.build}.{sample.release}.dna/stats/{sample.sample_id}{ext}",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            ext=[
                ".alignment_summary_metrics",
                ".insert_size_metrics",
                ".insert_size_histogram.pdf",
                ".base_distribution_by_cycle_metrics",
                ".base_distribution_by_cycle.pdf",
                ".gc_bias.detail_metrics",
                ".gc_bias.summary_metrics",
                ".gc_bias.pdf",
            ],
        ),
        idxstats=collect(
            "tmp/fair_bowtie2_mapping_samtools_idxstats/{sample.species}.{sample.build}.{sample.release}.{datatype}/{sample.sample_id}.idxstats",
            sample=lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}'",
                within=samples,
            ),
            allow_missing=True,
        ),
        fastp_pair_ended=collect(
            "tmp/fair_bowtie2_mapping_fastp_trimming_pair_ended/{sample.sample_id}.fastp.json",
            sample=lookup(
                query="downstream_file == downstream_file & species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        fastp_single_ended=collect(
            "tmp/fair_bowtie2_mapping_fastp_trimming_single_ended/{sample.sample_id}.fastp.json",
            sample=lookup(
                query="downstream_file != downstream_file & species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        fastqc_pair_ended=collect(
            "results/QC/report_pe/{sample.sample_id}.{stream}_fastqc.zip",
            sample=lookup(
                query="downstream_file == downstream_file & species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            stream=stream_tuple,
        ),
        fastqc_single_ended=collect(
            "results/QC/report_pe/{sample.sample_id}_fastqc.zip",
            sample=lookup(
                query="downstream_file != downstream_file & species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        bowtie2=collect(
            "logs/fair_bowtie2_mapping_bowtie2_alignment/{sample.species}.{sample.build}.{sample.release}.{datatype}/{sample.sample_id}.log",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            allow_missing=True,
        ),
        samtools=collect(
            "tmp/fair_bowtie2_mapping_samtools_stats/{sample.species}.{sample.build}.{sample.release}.{datatype}/{sample.sample_id}.txt",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            allow_missing=True,
        ),
        rseqc_infer_experiment=collect(
            "tmp/fair_bowtie2_mapping_rseqc_infer_experiment/{sample.species}.{sample.build}.{sample.release}.{datatype}/{sample.sample_id}.infer_experiment.txt",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            allow_missing=True,
        ),
        rseqc_bamstat=collect(
            "tmp/fair_bowtie2_mapping_rseqc_bamstat/{sample.species}.{sample.build}.{sample.release}.{datatype}/{sample.sample_id}.bamstat.txt",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            allow_missing=True,
        ),
        rseqc_read_gc=collect(
            "tmp/fair_bowtie2_mapping_rseqc_read_gc/{sample.species}.{sample.build}.{sample.release}.{datatype}/{sample.sample_id}.GC.xls",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            allow_missing=True,
        ),
        rseqc_read_distribution=collect(
            "tmp/fair_bowtie2_mapping_rseqc_read_distribution/{sample.species}.{sample.build}.{sample.release}.{datatype}/{sample.sample_id}.txt",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            allow_missing=True,
        ),
        rseqc_inner_distance=collect(
            "tmp/fair_bowtie2_mapping_rseqc_inner_distance/{sample.species}.{sample.build}.{sample.release}.{datatype}/{sample.sample_id}.inner_distance_freq.txt",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            allow_missing=True,
        ),
        goleft_indexcov_ped=collect(
            "tmp/fair_bowtie2_mapping_goleft_indexcov/{sample.species}.{sample.release}.{sample.build}.{datatype}/{sample.sample_id}-indexcov.ped",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            allow_missing=True,
        ),
        goleft_indexcov_roc=collect(
            "tmp/fair_bowtie2_mapping_goleft_indexcov/{sample.species}.{sample.release}.{sample.build}.{datatype}/{sample.sample_id}-indexcov.roc",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            allow_missing=True,
        ),
        ngsderive_readlen=collect(
            "tmp/fair_bowtie2_mapping_ngsderive_readlen/{sample.species}.{sample.build}/{sample.release}.{datatype}/{sample.sample_id}.readlen.tsv",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            allow_missing=True,
        ),
        ngsderive_instrument=collect(
            "tmp/fair_bowtie2_mapping_ngsderive_instrument/{sample.species}.{sample.build}/{sample.release}.{datatype}/{sample.sample_id}.instrument.tsv",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            allow_missing=True,
        ),
        ngsderive_encoding=collect(
            "tmp/fair_bowtie2_mapping_ngsderive_encoding/{sample.species}.{sample.build}/{sample.release}.{datatype}/{sample.sample_id}.encoding.tsv",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            allow_missing=True,
        ),
        mosdepth_global=collect(
            "tmp/fair_bowtie2_mapping_mosdepth/{sample.species}.{sample.build}.{sample.release}.{datatype}/{sample.sample_id}.mosdepth.global.dist.txt",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            allow_missing=True,
        ),
        mosdepth_region=collect(
            "tmp/fair_bowtie2_mapping_mosdepth/{sample.species}.{sample.build}.{sample.release}.{datatype}/{sample.sample_id}.mosdepth.region.dist.txt",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            allow_missing=True,
        ),
        mosdepth_summary=collect(
            "tmp/fair_bowtie2_mapping_mosdepth/{sample.species}.{sample.build}.{sample.release}.{datatype}/{sample.sample_id}.mosdepth.summary.txt",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            allow_missing=True,
        ),
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
        mem_mb=lambda wildcards, attempt: 7_000 + (1_000 * attempt),
        runtime=lambda wildcards, attempt: 30 * attempt,
        tmpdir=tmp,
    params:
        extra=lookup_config(
            dpath="params/fair_bowtie2_mapping_multiqc",
            default="--verbose --no-megaqc-upload --no-ansi --force",
        ),
        use_input_files_only=True,
    log:
        "logs/fair_bowtie2_mapping_multiqc_report/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_multiqc_report/{species}.{build}.{release}.{datatype}.tsv"
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/multiqc"
