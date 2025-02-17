Material and methods
=====================

Genome DNA sequence and annotations were download from Ensembl. 
Pyfaidx_ [#pyfaidxpaper]_ was used to filter non-cannonical 
chromosomes. Agat_ [#agatpaper]_ was used to correct common 
issues found in Ensembl genome annotation files, filter non-
cannonical chromosomes, and remove transcripts with TSL being
equal to NA. Samtools_ [#samtoolspaper]_ and Picard_ [#picardpaper]_ 
were used to index genome sequences.

Raw fastq file quality was assessed with FastQC_ [#fastqcpaper]_.
FastQ files formats were checked using Fastq_utils_ [#fastqutilspaper]_, additional
qualities were gathered using Seqkit_ [#seqkitpaper]_ and Fastqinfo_ [#fastqinfopaper]_.

Raw fastq files were trimmed using Fastp_ [#fastppaper]_ . Cleaned 
reads were aligned over indexed Ensembl genome with Bowtie2_ 
[#bowtie2paper]_. Sambamba_ [#sambambapaper]_ was used to sort, 
filter, mark duplicates, and compress aligned reads. Quality 
controls were done on cleaned, sorted, deduplicated aligned reads 
using Picard_ [#picardpaper]_ and Samtools_ [#samtoolspaper]_.
Additonal quality assessments are done with RSeQC_ [#rseqcpaper]_,
NGSderive_ [#ngsderivepaper]_, GOleft_ [#goleftpaper]_, and
Mosdepth_ [#mosdepthpaper]_ to determine mapping quality. Additional
sample quality checks were done with MTNucRatioCalculator_ [#mtnucratiocalculator]_,
and SexDetERRmine_ [#sexdeterrmine]_ to verify library and sample
expectations.
Quality repord produced during both trimming and mapping steps 
have been aggregated with MultiQC_ [#multiqcpaper]_. 

On user demand, alignment sieve are produced using Deeptools_ [#deeptoolspaper]_.

The whole pipeline_ [#fairbowtiemapping]_ was powered by Snakemake_ [#snakemakepaper]_. 
This pipeline is freely available on Github_, details about 
installation usage, and resutls can be found on the 
`Snakemake workflow`_ page. It relies on `fair-genome-indexer`_ [#fairgenomeindexer]_ and
`fair-fastqc-multiqc`_ [#fairfastqcmultiqc]_ pipelines.

.. [#pyfaidxpaper] Shirley, Matthew D., et al. Efficient" pythonic" access to FASTA files using pyfaidx. No. e1196. PeerJ PrePrints, 2015.
.. [#agatpaper] Dainat J. AGAT: Another Gff Analysis Toolkit to handle annotations in any GTF/GFF format.  (Version v0.7.0). Zenodo. https://www.doi.org/10.5281/zenodo.3552717
.. [#samtoolspaper] Li, Heng, et al. "The sequence alignment/map format and SAMtools." bioinformatics 25.16 (2009): 2078-2079.
.. [#picardpaper] McKenna, Aaron, et al. "The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data." Genome research 20.9 (2010): 1297-1303.
.. [#fastqcpaper] Andrews, S. Fastqc. "A quality control tool for high throughput sequence data. Augen, J.(2004). Bioinformatics in the post-genomic era: Genome, transcriptome, proteome, and information-based medicine." (2010).
.. [#fastqutilspaper] Nuno Fonseca, & Jonathan Manning. (2023). nunofonseca/fastq_utils: 0.25.2 (0.25.2). Zenodo. https://doi.org/10.5281/zenodo.7755574
.. [#seqkitpaper] Shen, Wei, et al. "SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation." PloS one 11.10 (2016): e0163962.
.. [#fastqinfopaper] Kiu R, fastq-info: compute estimated sequencing depth (coverage) of prokaryotic genomes
.. [#fastppaper] Chen, Shifu, et al. "fastp: an ultra-fast all-in-one FASTQ preprocessor." Bioinformatics 34.17 (2018): i884-i890.
.. [#bowtie2paper] Langmead, Ben, and Steven L. Salzberg. "Fast gapped-read alignment with Bowtie 2." Nature methods 9.4 (2012): 357-359.
.. [#sambambapaper] Tarasov, Artem, et al. "Sambamba: fast processing of NGS alignment formats." Bioinformatics 31.12 (2015): 2032-2034.
.. [#rseqcpaper] Wang, Liguo, Shengqin Wang, and Wei Li. "RSeQC: quality control of RNA-seq experiments." Bioinformatics 28.16 (2012): 2184-2185.
.. [#ngsderivepaper] McLeod, Clay, et al. "St. Jude Cloud: a pediatric cancer genomic data-sharing ecosystem." Cancer discovery 11.5 (2021): 1082-1099.
.. [#goleftpaper] Pedersen, Brent S., et al. "Indexcov: fast coverage quality control for whole-genome sequencing." Gigascience 6.11 (2017): gix090.
.. [#mosdepthpaper] Pedersen, Brent S., and Aaron R. Quinlan. "Mosdepth: quick coverage calculation for genomes and exomes." Bioinformatics 34.5 (2018): 867-868.
.. [#mtnucratiocalculator] Alexander Peltzerand Alexander Regueiro. Apeltzer/mtnucratiocalculator: 0.7.1. 0.7.1, Zenodo, 10 Sept. 2024, doi:10.5281/zenodo.13739840.
.. [#sexdeterrmine] Lamnidis, T.C. et al., 2018. Ancient Fennoscandian genomes reveal origin and spread of Siberian ancestry in Europe. Nature communications, 9(1), p.5018. Available at: http://dx.doi.org/10.1038/s41467-018-07483-5.
.. [#multiqcpaper] Ewels, Philip, et al. "MultiQC: summarize analysis results for multiple tools and samples in a single report." Bioinformatics 32.19 (2016): 3047-3048.
.. [#deeptoolspaper] Ramírez, Fidel, et al. "deepTools: a flexible platform for exploring deep-sequencing data." Nucleic acids research 42.W1 (2014): W187-W191.
.. [#fairbowtiemapping] Dayris, T. (2024). fair-bowtie2-mapping (Version 4.4.1) [Computer software]. https://github.com/tdayris/fair_bowtie2_mapping
.. [#snakemakepaper] Köster, Johannes, and Sven Rahmann. "Snakemake—a scalable bioinformatics workflow engine." Bioinformatics 28.19 (2012): 2520-2522.
.. [#fairgenomeindexer] Dayris, T. (2024). fair-genome-indexer (Version 3.9.3) [Computer software]. https://github.com/tdayris/fair_genome_indexer
.. [#fairfastqcmultiqc] Dayris, T. (2024). fair-fastqc-multiqc (Version 2.4.2) [Computer software]. https://github.com/tdayris/fair_fastqc_multiqc


.. _Sambamba: https://snakemake-wrappers.readthedocs.io/en/v5.6.0/wrappers/sambamba.html
.. _Bowtie2: https://snakemake-wrappers.readthedocs.io/en/v5.6.0/wrappers/bowtie2.html
.. _Fastp: https://snakemake-wrappers.readthedocs.io/en/v5.6.0/wrappers/fastp.html
.. _Mosdepth: https://snakemake-wrappers.readthedocs.io/en/v5.6.0/wrappers/mosdepth.html
.. _Picard: https://snakemake-wrappers.readthedocs.io/en/v5.6.0/wrappers/picard/collectmultiplemetrics.html
.. _MultiQC: https://snakemake-wrappers.readthedocs.io/en/v5.6.0/wrappers/multiqc.html
.. _Snakemake: https://snakemake.readthedocs.io
.. _Github: https://github.com/tdayris/fair_bowtie2_mapping
.. _`Snakemake workflow`: https://snakemake.github.io/snakemake-workflow-catalog?usage=tdayris/fair_bowtie2_mapping
.. _Agat: https://agat.readthedocs.io/en/latest/index.html
.. _Samtools: https://snakemake-wrappers.readthedocs.io/en/v5.6.0/wrappers/samtools/faidx.html
.. _FastQC: https://snakemake-wrappers.readthedocs.io/en/v5.6.0/wrappers/fastqc.html
.. _Fastq_utils: https://github.com/nunofonseca/fastq_utils
.. _Seqkit: https://bioinf.shenwei.me/seqkit/
.. _Fastqinfo: https://github.com/raymondkiu/fastq-info
.. _Pyfaidx: https://github.com/mdshw5/pyfaidx
.. _GOleft: https://github.com/brentp/goleft
.. _NGSderive: https://stjudecloud.github.io/ngsderive/
.. _RSeQC: https://rseqc.sourceforge.net/
.. _Deeptools: https://snakemake-wrappers.readthedocs.io/en/v5.6.0/wrappers/deeptools/alignmentsieve.html
.. _SexDetERRmine: ???
.. _MTNucRatioCalculator: ???
.. _`fair-genome-indexer`: https://github.com/tdayris/fair_genome_indexer
.. _`fair-fastqc-multiqc`: https://github.com/tdayris/fair_fastqc_multiqc
.. _`fair-bowtie2-mapping`: https://github.com/tdayris/fair_bowtie2_mapping

:Authors:
    Thibault Dayris


:Version: 4.4.5 of 2025-02-17
