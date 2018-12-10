#####HW4

Homework 4-Ghaidaa Kashgari 
Summarize partitions of a genome assembly
First, download necessary modules.

Prepare environment
$ module load jje/jjeutils jje/kent
Calculate the following for all sequences ≤ 100kb and all sequences > 100kb:
Total number of nucleotides
Total number of Ns
Total number of sequences
Download most current fasta file for D. melanogaster wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/dmel-all-chromosome-r6.24.fasta.gz

Unzip the file
gunzip dmel-all-chromosome-r6.24.fasta.gz

shorter sequences
bioawk -c fastx '{ if(length($seq) < 100000) { print ">"$name; print $seq }}' dmel-all-chromosome-r6.24.fasta > shorter.kb.fasta

faSize shorter.kb.fasta

Answer:

Total number of nucleotides = 6178042
Total number of Ns = 662593
Total number of sequences = 1863
6178042 bases (662593 N's 5515449 real 5515449 upper 0 lower) in 1863 sequences in 1 files
Total size: mean 3316.2 sd 7116.2 min 544 (211000022279089) max 88768 (Unmapped_Scaffold_8_D1580_D1567) median 1567
N count: mean 355.7 sd 1700.6
U count: mean 2960.5 sd 6351.5
L count: mean 0.0 sd 0.0
%0.00 masked total, %0.00 masked real

longer sequences
bioawk -c fastx '{ if(length($seq) > 100000) { print ">"$name; print $seq }}' dmel-all-chromosome-r6.24.fasta > longer.kb.fasta

faSize longer.kb.fasta

Answer:

Total number of nucleotides = 137547960
Total number of Ns = 490385
Total number of sequences = 7
137547960 bases (490385 N's 137057575 real 137057575 upper 0 lower) in 7 sequences in 1 files
Total size: mean 19649708.6 sd 12099037.5 min 1348131 (4) max 32079331 (3R) median 23542271
N count: mean 70055.0 sd 92459.2
U count: mean 19579653.6 sd 12138278.9
L count: mean 0.0 sd 0.0
%0.00 masked total, %0.00 masked real

Plots
Sequence length distribution
Sequence GC% distribution
Cumulative genome size sorted from largest to smallest sequences
###Generating Sequence Length Distribution File

Rename fasta file "assembly"
mv dmel-all-chromosome-r6.24.fasta assembly.txt

Get just headers from all chromosome "assembly" file grep ">" assembly.txt > headers.txt

View the field with lengths
awk '{ print $7 }' headers.txt > lengths.txt

Remove characters and keep numerical length values
tr -d "length=;" < lengths.txt > lengthsnum.txt

Sort length from longest to shortest sort -rn lengthsnum.txt > sortlength.txt

Add Length Header
echo $'Length' | cat - sortlength.txt > sortlengthhead.txt

For shorter sequences,

	library(ggplot2)

	lengthdistribution <-	read.table("shorter_sort.txt", header = TRUE)
	LD = ggplot(data = lengthdistribution)
	LD + geom_histogram(mapping = aes(x=Length), bins 	= 1000) + scale_x_log10() + ggtitle("Sequence 	Length Distibution") + xlab("Sequence Length (bp)")
� ###Generating Sequence Length Distribution Plot

Load RStudio
module load rstudio

Open rstudio
rstudio

Load ggplot
library(ggplot2)

Load data into rstudio
lengthdistribution <- read.table("sortlengthhead.txt", header = TRUE)

Load data into ggplot
LD = ggplot(data = lengthdistribution)

Make labeled histogram out of data, with a log scale and bin sizes of 1000
LD + geom_histogram(mapping = aes(x=Length), bins = 1000) + scale_x_log10() + ggtitle("Sequence Length Distibution") + xlab("Sequence Length (bp)")

###Generating Sequence GC% Distribution File

View the GC content bioawk -c fastx '{ print ">"$name; print gc($seq) }' assembly.txt > GC.txt

gets just the GC% of each sequence grep -v ">" GC.txt > GCval.txt

sorts GC content from highest to lowest
sort -rn GCval.txt > GCsort.txt

shows answer less GCsort.txt

###Generating Sequence GC% Distribution Plot

Open rstudio
rstudio

Load ggplot
library(ggplot2)

Load GCsort data into Rstudio GCDistribution = read.table("GCsort.txt")

Load data into ggplot
GCD = ggplot(data = GCDistribution)

Make labeled histogram with bin sizes of 0.01% GCD + geom_histogram(mapping = aes(x=V1), binwidth = 0.01) + ggtitle("Sequence GC% Distribution") + xlab("Sequence GC Content (%)")

Save plot ggsave("SeqGCDis.png", width = 6, height = 6)

###Generating CDF File

Adds Dmel Assembly ID
bioawk '{print $1 "\t Dmel"}' < sortlength.txt > CDF.txt

Adds Length and Assembly Headers to list of sorted Lengths
echo $'Length\tAssembly' | cat - CDF.txt > CDFHead.txt

create CDF from file with headers Length and Assembly
plotCDF2 CDFHead.txt GenomeCDF.png

Plots for longer sequences
bioawk -c fastx '{ print $name, length($seq) }' longer.kb.fasta > lengths.txt

big <-read.table("lengths.txt", header=TRUE)

big$cut <-cut(x=big$Length, breaks=50)
	
bigplot <-ggplot(data=big)

bigplot + geom_bar(mapping = aes(x=cut)) + labs(title="Sequence Length Distribution", x="Length", y="Count (Number of Contigs)") + theme_bw()+ theme(axis.text.x = element_text(angle = 60, hjust = 1))
Genome assembly
Assemble a genome from MinION reads1.
Go to an open interactive queue and request 32 cores
qrsh -q abio,free128,free88i,free72i,free32i,free64 -pe openmp 32

Download reads
wget https://hpc.oit.uci.edu/~solarese/ee282/iso1_onp_a2_1kb.fastq.gz

Unzip the fasta file
gunzip iso1_onp_a2_1kb.fastq.gz

Use minimap to overlap reads minimap -t 32 -Sw5 -L100 -m0 iso*.fastq iso*.fastq | gzip -1 > iso1_onp.paf.gz

Use miniasm to construct an assembly
miniasm -f iso1_onp_a2_1kb.fastq iso1_onp.paf.gz > reads.gfa

Convert gfa file to fasta file
awk '/^S/{print ">"$2"\n"$3}' reads.gfa | fold -w 60 > unitigs.fa

Assembly Assessment
Calculate N50
bioawk -c fastx ' { li=length($seq); l=li+l; print li; } END { print l; } ' unitigs.fa \ | sort -rn \ | gawk 'NR == 1 { l=$1; } NR > 1 { li=$1; lc=li+lc; if(lc/l >= 0.5) { print li; exit; } } '\ | less -S

####N50=4,494,246 (Reference N50=21,485,538)

Compare your assembly to the contig assembly (not the scaffold assembly!) from Drosophila melanogaster on FlyBase
faSplitByN dmel-all-chromosome-r6.24.fasta.gz unitigs.fa 10

bioawk -c fastx ' {li=length($seq); l=li+l; print li; } END {print l;} ' unitigs.fa | sort -rn | gawk ' NR == 1 {l = $1;} NR >1 {li=$1; lc=li+lc; if(lc/l >= 0.5) {print li; exit;} }' | less -S

21485538

MuMmer plot script
#!/bin/bash
#$ -N MUMmer_hw4
#$ -q abio128,abio,bsg,bsg2
#$ -pe openmp 32-128
#$ -R Y
#$ -m beas

###Loading of binaries via module load or PATH 	reassignment
source /pub/jje/ee282/bin/.qmbashrc
module load gnuplot/4.6.0

###Query and Reference Assignment. State my 	prefix for output filenames
REF="/data/users/bcraver/myrepos/HW4/assembly/	unitigs.fa"
PREFIX="flybase"
SGE_TASK_ID=1
QRY=$(ls /pub/jje/ee282/bcraver/nanopore_assembly/nanopore_assembly/data/processed/	unitigs.fa | head -n $SGE_TASK_ID | tail -n 1)
PREFIX=${PREFIX}_$(basename ${QRY} .fa)

#nucmer -l 100 -c 125 -d 10 -banded -D 5 -prefix 	${PREFIX} ${REF} ${QRY}
mummerplot --fat --layout --filter -p ${PREFIX} $	{PREFIX}.delta \
  -R ${REF} -Q ${QRY} --png
Compare your assembly to both the contig assembly and the scaffold assembly from the Drosophila melanogaster on FlyBase using a contiguity plot.
#!/usr/bin/env bash
module load rstudio/0.99.9.9 
module load perl 
module load jje/jjeutils/0.1a 
module load jje/kent

bioawk -c fastx ' { print length($seq) } ' /pub/jje/ee282/	bcraver/nanopore_assembly/nanopore_assembly/data/	processed/unitigs.fa \
| sort -rn \
| awk ' BEGIN { print "Assembly\tLength\nUnitigs_Ctg\t0" } 	{ print "Unitigs_Ctg\t" $1 } ' \
> unitigs.sorted.text

bioawk -c fastx ' { print length($seq) } ' dmel-all-	chromosome-r6.24.fasta.gz \
| sort -rn \
| awk ' BEGIN { print 	"Assembly\tLength\ndmel_scaffold\t0" } { print 	"dmel_scaffold\t" $1 } ' \
> dmel_scaffold.sorted.text

bioawk -c fastx ' { print length($seq) } ' dmel-all-	chromosome-r6.24.ctg.fa \
| sort -rn \
| awk ' BEGIN { print "Assembly\tLength\ndmel_ctg\t0" } 	{ print "dmel_ctg\t" $1 } ' \
> dmel_ctg.sorted.text

plotCDF2 *.sorted.text /dev/stdout \
| tee hw4_plotcdf2.png \
| display  
Calculate BUSCO score
#!/bin/bash
#
#$ -N hw4_busco
#$ -q free128,abio128,bio,abio
#$ -pe openmp 32
#$ -R Y
module load augustus/3.2.1
module load blast/2.2.31 hmmer/3.1b2 boost/1.54.0
source /pub/jje/ee282/bin/.buscorc

INPUTTYPE="geno"
MYLIBDIR="/data/users/bcraver/busco/lineages/"
MYLIB="diptera_odb9"
OPTIONS="-l ${MYLIBDIR}${MYLIB}"
##OPTIONS="${OPTIONS} -sp 4577"
QRY="/pub/jje/ee282/bcraver/nanopore_assembly/	nanopore_assembly/data/processed/unitigs.fa"
MYEXT=".fa" ###Please change this based on your qry file. 	I.e. .fasta or .fa or .gfa

#my busco run
BUSCO.py -c 128 ${NSLOTS} -i ${QRY} -m ${INPUTTYPE} -o $	(basename ${QRY} ${MYEXT})_${MYLIB}${SPTAG} ${OPTIONS} -f"
C:0.5%[S:0.5%,D:0.0%],F:1.1%,M:98.4%,n:2799
INFO	13 Complete BUSCOs (C)
INFO	13 Complete and single-copy BUSCOs (S)
INFO	0 Complete and duplicated BUSCOs (D)
INFO	32 Fragmented BUSCOs (F)
INFO	2754 Missing BUSCOs (M)
INFO	2799 Total BUSCO groups searched
