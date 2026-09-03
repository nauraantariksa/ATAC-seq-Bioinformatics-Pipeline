This GitHub page is a collection of scripts and general code used for the analysis of ATAC-seq datasets

The pipeline of ATAC analyses is identical to that of ChIP-seq, except during **Filtering&Sorting and Deduplication**, where ChrM reads are removed. Additionally, peak calling with MACS2 uses different parameters. As such, this page will only elaborate on those differences. For:
1) Trimming
2) QC
3) Alignment
please refer to the ChIP-Seq Bioinformatics Pipeline

Make sure these dependencies are installed in your system:

anaconda3
python (version 3.9)
pairtools (https://github.com/open2c/pairtools)
samtools
bedtools

In anaconda3, make separate environments for the following Python suites used for the data analyses:

deeptools (in a deeptools_env, perform $conda install -c conda-forge -c bioconda deeptools)
macs2 (in a macs2_env, perform $conda install -c bioconda macs2)

For other analyses, install or clone Github repositories:

R
Homer
SEACR
