SHELL=/usr/bin/bash
.ONESHELL:
.SHELLFLAGS := -euic
.DELETE_ON_ERROR: linter_info.txt format.txt pipeline_pair_ended.txt pipeline_single_ended.txt pipeline_complex.txt

MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

# List of snakefiles to lint
snakefiles := \
	"../workflow/rules/bowtie2_sambamba_meta.smk" \
	"../workflow/rules/multiqc.smk" \
	"../workflow/rules/picard_metrics.smk" \
	"../workflow/rules/fair_genome_indexer_pipeline.smk" \
	"../workflow/rules/common.smk" \
	"../workflow/Snakefile"


.PHONY: all
all: pipeline_pair_ended.txt pipeline_single_ended.txt pipeline_complex.txt


pipeline_pair_ended.txt: format.txt linter.txt
	mamba activate test-snakemake-wrappers && \
	snakemake -s ../workflow/Snakefile \
		--configfile config_pair_ended/config.yaml \
		--cores 1 \
		--use-conda \
		--use-singularity \
		--workflow-profile ../profiles 2>&1 | tee pipeline_pair_ended.txt


pipeline_single_ended.txt: format.txt linter.txt
	mamba activate test-snakemake-wrappers && \
	snakemake -s ../workflow/Snakefile \
		--configfile config_single_ended/config.yaml \
		--cores 1 \
		--use-conda \
		--use-singularity \
		--workflow-profile ../profiles 2>&1 | tee pipeline_pair_ended.txt


pipeline_complex.txt: format.txt linter.txt
	mamba activate test-snakemake-wrappers && \
	snakemake -s ../workflow/Snakefile \
		--configfile config_complex/config.yaml \
		--cores 1 \
		--use-conda \
		--use-singularity \
		--workflow-profile ../profiles 2>&1 | tee pipeline_pair_ended.txt


format.txt:
	mamba activate test-snakemake-wrappers && \
	snakefmt $(snakefiles) 2>&1 | tee format.txt


linter.txt:
	mamba activate test-snakemake-wrappers && \
	snakemake -s ../workflow/Snakefile --lint 2>&1 | tee linter.txt


.PHONY: clean
clean:
	rm -rv linter.txt \
		   format.txt \
		   pipeline_pair_ended.txt \
		   pipeline_single_ended.txt \
		   pipeline_complex.txt \
		   reference results \
		   .snakemake