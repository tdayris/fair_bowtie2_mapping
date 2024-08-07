SHELL=/usr/bin/bash
.ONESHELL:
.SHELLFLAGS := -euic
.DELETE_ON_ERROR: linter_info.txt format.txt pipeline.txt

MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

# List of snakefiles to lint
snakefiles := \
	"../workflow/rules/bowtie2_sambamba_meta.smk" \
	"../workflow/rules/common.smk" \
	"../workflow/rules/deeptools.smk" \
	"../workflow/rules/fair_genome_indexer.smk" \
	"../workflow/rules/fair_fastqc_multiqc.smk" \
	"../workflow/rules/fastp.smk" \
	"../workflow/rules/goleft.smk" \
	"../workflow/rules/multiqc.smk" \
	"../workflow/rules/ngsderive.smk" \
	"../workflow/rules/picard_metrics.smk" \
	"../workflow/rules/rseqc.smk" \
	"../workflow/rules/samtools_stats.smk" \
	"../workflow/rules/mosdepth.smk" \
	"../workflow/Snakefile"


pythonfiles := \
	"../workflow/scripts/fair_bowtie2_mapping_multiqc_config.py" \
	"../workflow/scripts/regenerate_genome.py" \
	"../workflow/scripts/rseqc_bamstat.py" \
	"../workflow/scripts/rseqc_inner_distance.py" \
	"../workflow/scripts/ucsc_genepred_to_bed.py"


.PHONY: all
all: report.txt


dag.png:
	snakemake -s ../workflow/Snakefile \
		--rulegraph --forceall | dot -Tpng > dag.png && \
	cp --verbose dag.png ../dag.png


report.txt: pipeline.txt
	snakemake -s ../workflow/Snakefile --profile local \
		--report report.html 2>&1 | tee report.txt


pipeline.txt: linter.txt
	snakemake -s ../workflow/Snakefile --profile local  \
		--cores 2 \
		--software-deployment-method apptainer conda \
		2>&1 | tee pipeline.txt


format.txt:
	snakefmt $(snakefiles) 2>&1 | tee format.txt


linter.txt: format.txt black.txt
	snakemake -s ../workflow/Snakefile  --profile local \
		--lint 2>&1 | tee linter.txt


black.txt:
	black $(pythonfiles) 2>&1 | tee black.txt



.PHONY: clean
clean:
	rm -rvf linter.txt \
		   format.txt \
		   reference \
		   results \
		   tmp \
		   logs \
		   benchmark \
		   report.html \
		   report.txt \
		   pipeline.txt \
		   genomes.csv \
		   dag.png \
		   .cache \
		   .conda \
		   .condarc \
		   .java \
		   .wget-hsts \
		   .snakemake \
		   black.txt \
		   report.zip \
		   test.png
