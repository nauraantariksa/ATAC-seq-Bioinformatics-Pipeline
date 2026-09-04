# Table of Contents
- [Data analyses](#data-analyses)
  - [List of tools/software](#list-of-toolssoftware)
  - [Filtering & Sorting](#filtering--sorting)
  - [Peak Calling with MACS2](#peak-calling-with-macs2)
  - [Differential chromatin accessibility analyses (CSB+ vs CSB-)](#differential-chromatin-accessibility-analyses-csb-vs-csb-)

# Data analyses
This GitHub page is a collection of scripts and general code used for the analysis of ATAC datasets. Differential analyses were given with a comparison between CS1AN CSB+ and CSB- cell lines as an example.

## List of tools/software
The pipeline of ATAC analyses is identical to that of ChIP-seq, except during Filtering&Sorting and Deduplication, where ChrM reads are removed. Additionally, peak calling with MACS2 uses different parameters. As such, this page will only elaborate on those differences. For trimming, QC, and alignment, please refer to the ChIP-Seq Bioinformatics Pipeline in https://github.com/nauraantariksa/ChIP-Bioinformatics-Pipeline

Make sure these dependencies are installed in your system:
- anaconda3
- python (version 3.9)
- fastp v0.23.4
- FastQC v0.11.9
- Bowtie2 v2.5.1
- SAMtools v1.22.1
- Picard v2.25.1
- BEDTools v2.31.1

In anaconda3, make separate environments for the following Python suites used for the data analyses:
- deepTools v3.5.6
- MACS2 v2.2.9.1

For other analyses, install or clone Github repositories:
- R
- HOMER v5.1

For visualization:
- IGV v2.19.8

## Filtering & Sorting
ATAC data requires the removal of mitochondrial reads
```bash
#PBS -lwalltime=72:00:00
#PBS -lselect=1:ncpus=128:mem=900gb

module load SAMtools/1.18-GCC-12.3.0
module load BEDTools/2.31.0-GCC-12.3.0
module load tools/prod
module load picard/2.25.1-Java-11

sam_dir=/PATH/TO/SAM/FILES
tmp_dir=$EPHEMERAL/tmp_filtering
mkdir -p $tmp_dir

for file in $sam_dir/*.sam
do
name=`basename $file .sam`
samtools view -F 2308 -b -q 10 $sam_dir/$name.sam > $tmp_dir/$name.bam
samtools sort -o $tmp_dir/$name.sort.bam -O 'bam' -T $tmp_dir/temp_$name $tmp_dir/$name.bam
java -jar $EBROOTPICARD/picard.jar MarkDuplicates REMOVE_DUPLICATES=TRUE I=$tmp_dir/$name.sort.bam O=$tmp_dir/$name.dedup.sort.bam M=$tmp_dir/$name.metrics.txt TMP_DIR=$tmp_dir
bedtools intersect -v -abam $tmp_dir/$name.dedup.sort.bam -b '/rds/general/project/diantonio/live/tools/blacklisted_regions/hg38-blacklist.v2.bed' > $tmp_dir/$name.bl.dedup.sort.bam
samtools index $tmp_dir/$name.bl.dedup.sort.bam
rm $tmp_dir/$name.bam $tmp_dir/$name.sort.bam $tmp_dir/$name.dedup.sort.bam
done

cp $tmp_dir/*.bl.dedup.sort.bam* $tmp_dir/*.metrics.txt "$PBS_O_WORKDIR"
```

## Peak Calling with MACS2
ATAC peaks were identified using MACS2 using the default settings with the paired-end mode on.
```bash
eval "$(/rds/general/user/nfa23/home/anaconda3/bin/conda shell.bash hook)"
conda activate macs2_env
macs2 callpeak -f BAMPE --keep-dup all -t ATAC_BIO1_NoMito.bam -n ATAC_BIO1
```

## Differential Chromatin Accessibility Analyses (CSB+ vs CSB-)
Creating a consensus peak requires cloning Robert Hansel-Hertsch's (2016) GitHub at https://github.com/sblab-bioinformatics/dna-secondary-struct-chrom-lands

```bash
## Consensus peaks in Cell Line 1
mergePeaks.sh CELL_1_BIO1.narrowPeak CELL_1_BIO2.narrowPeak CELL_1_BIO3.narrowPeak CELL_1_BIO4.narrowPeak
| awk '$5 > 2' > CELL_1.merge.narrowPeak # awk '$5 > 2' refers to peaks appearing at least 3/4 replicates

## Consensus peaks in Cell Line 2
mergePeaks.sh CELL_2_BIO1.narrowPeak CELL_2_BIO2.narrowPeak CELL_2_BIO3.narrowPeak CELL_2_BIO4.narrowPeak
| awk '$5 > 2' > CELL_2.merge.narrowPeak # awk '$5 > 2' refers to peaks appearing at least 3/4 replicates

## Union set of Cell Line 1 & 2
mergePeaks.sh CELL_1.merge.narrowPeak CELL_2.merge.narrowPeak \
| sortBedAsBam.py -i - -b <BAM> > CELL_1-2.union.merge.bed

## Read counts in peaks
for bam in *.bam
do
coverageBed -g genome.txt -counts -a CELL_1-2.union.merge.bed -b $bam > ${name}.union.bed
done

echo 'chrom start end peakFiles nfiles nreads sample_id' | tr ' ' '\t' > CELL_1-2.union.merge.counts.bed
tableCat.py -i *.union.bed -r '.union.bed' >> CELL_1-2.union.merge.counts.bed
```

Differential analyses were performed locally on R, using code adapted from Robert Hansel-Hertsch (2016)
```bash
library(data.table)
library(edgeR)
library(reshape2)
cnt<- fread('/Users/naura/Downloads/CS1AN_WTvCSB_ATAC.union.merge.counts.bed')
cnt <- cnt[chrom != 'chrY']
cnt[, locus := paste(chrom, start, end, sep= '_')]
cntct<- dcast.data.table(data= cnt, locus ~ sample_id, value.var= 'nreads') 
y<- data.frame(cntct[, 2:ncol(cntct), with= FALSE])
row.names(y)<- cntct$locus

# Lib size *without* chrM from above:
libSize<- c(
  CS1AN_CSB_ATAC_Bio1 = 83485375,
  CS1AN_CSB_ATAC_Bio2 = 95964802,
  CS1AN_CSB_ATAC_Bio3 = 81124923,
  CS1AN_CSB_ATAC_Bio4 = 81934961,
  CS1AN_WT_ATAC_Bio1 = 76845250,
  CS1AN_WT_ATAC_Bio2 = 87851889,
  CS1AN_WT_ATAC_Bio3 = 86445406,
  CS1AN_WT_ATAC_Bio4 = 83731461)

group <- factor(c('CSB', 'CSB', 'CSB', 'CSB', 'WT', 'WT', 'WT', 'WT'))

y<- DGEList(counts=y, group=group)
stopifnot(rownames(y$samples) == names(libSize))
y$samples$lib.size<- libSize
y<- calcNormFactors(y, method= 'none')
y<- estimateDisp(y)
y<- estimateCommonDisp(y)
y<- estimateTagwiseDisp(y)

et<- exactTest(y, pair = c('CSB', 'WT'))

detable<- data.frame(topTags(et, n= Inf)$table)
detable$locus<- rownames(detable)
detable<- data.table(detable)
detable<- merge(detable, unique(cnt[, list(chrom, start, end, locus)]), by= 'locus')

# Exclude sex chromosomes due to clonal Y-loss artifact between CSB/WT lines
#detable <- detable[!chrom %in% c('chrX', 'chrY')]

pal<- colorRampPalette(c("white", "lightblue", "yellow", "red"), space = "Lab")
pdf('maplot.atac-CSBvWT.pdf', w= 12/2.54, h= 12/2.54, pointsize= 10)
par(las= 1, mgp= c(1.75, 0.5, 0), bty= 'l', mar= c(3, 3, 3, 0.5))
smoothScatter(x= detable$logCPM, y= detable$logFC, xlab= 'logCPM', ylab= 'logFC',
              main= "Differential chromatin (ATAC)\n[CSB - WT]", colramp= pal, col= 'blue')
lines(loess.smooth(x= detable$logCPM, y= detable$logFC, span= 0.1), lwd= 2, col= 'grey60')
abline(h= 0, col= 'grey30')
points(x= detable$logCPM, y= detable$logFC, col= ifelse(detable$FDR < 0.05, '#FF000080', 'transparent'), cex= 0.5, pch= '.')
mtext(side= 3, line= -1.2, text= sprintf('FDR < 0.05: %s', nrow(detable[FDR < 0.05 & logFC > 0])), adj= 1)
mtext(side= 1, line= -1.2, text= sprintf('FDR < 0.05: %s', nrow(detable[FDR < 0.05 & logFC < 0])), adj= 1)
grid(col= 'grey50')
dev.off()

write.table(detable, "diff.atac-CSBvWT.txt", row.names= FALSE, col.names= TRUE, sep= '\t', quote= FALSE)
```
