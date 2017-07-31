#!/usr/bin/Rscript
#PBS -l nodes=1:ppn=8 -j oe -N TruePaiR
#
# Example run:
# Just produce stats: 
# Rscript TruePaiR.R EXAMPLE.fq
#
# Keep intermediate files (can check individual read pairs):
# Rscript TruePaiR.R EXAMPLE.fq TRUE

ags <- commandArgs(TRUE)
fq<-ags[1]; verbose<-ags[2]
options(warn=-1)

library(ShortRead); library(seqinr)

if(is.na(verbose)) { verbose="FALSE" }

suppressPackageStartupMessages(library(ShortRead)); suppressPackageStartupMessages(library(seqinr)); suppressPackageStartupMessages(library(systemPipeR))
moduleload("bowtie2")

redundCorrect <-
function(ids, u1a10_bank=as.character())
{
#	Returns the appropriate value for the number of new unique reads
#	that has been shown to have a pair (One hit, one pair)
#
#	sam_df - SAM file read in as a data frame
#	u1a10_bank - optional: use if reads have already been associated with pairs
#
#	For map hits, the value (per line in SAM file) will always be:
#	0 - both reads have already been associated with a pair
#	1 - one, but not both reads has yet to be associated with a pair
#	2 - neither read has previously been associated with a pair

	u1_redund <- as.character(ids) %in% as.character(u1a10_bank)
	u1_redund <- gsub(TRUE, 1, u1_redund)
	u1_redund <- gsub(FALSE, 0, u1_redund)
	count <- sum(as.integer(u1_redund))

	return(count)
}

printSumm <- 
function(num_u1, num_a10, perc_u1, perc_a10, and1, and2, or, verbose)
{
	count=0

	# Imperfect Pairs: 1 count per read that has been associated with a U1/A10 pair (1-2 mismatches)
	if(file.info("and1.sam")$size != 0 | file.info("and2.sam")$size != 0) {
		system("cat and1.sam and2.sam > and-all.sam")

		ab <- read.table("and-all.sam", sep="\t", row.names=NULL, fill=TRUE)
		ab <- ab[!duplicated(ab),]
		ab1 <- ab[!duplicated(as.character(ab[,1])),]
		ab3 <- ab[!duplicated(as.character(ab[,3])),]	
	} else { system("cp and1.sam and-all.sam"); ab1 <- data.frame(); ab3 <- data.frame() }
	
	if(file.info("or.sam")$size != 0) {
		or <- read.table("or.sam", sep="\t", fill=TRUE)
		# Read can only appear once from U1 and A10 subsets
		or1 <- or[!duplicated(as.character(or[,1])),]
		or3 <- or[!duplicated(as.character(or[,3])),]

		if(length(ab1) != 0) { 
			or_redund1 <- redundCorrect(or1[,1], ab1[,1]) 
			or_redund3 <- redundCorrect(or3[,3], ab1[,1])
		} else {
			or_redund1 <- redundCorrect(or1[,1])
			or_redund3 <- redundCorrect(or3[,3])
		}
	} else { or1 <- data.frame(); or3 <- data.frame(); ab3 <- data.frame(); or_redund1 <- 0; or_redund3 <- 0; }
	
	if(length(ab1)==0) { ab1_len <- 0 } else { ab1_len <- length(ab1[,1]) }
	if(length(ab3)==0) { ab3_len <- 0 } else { ab3_len <- length(ab3[,1]) }
	if(length(or1)==0) { or1_len <- 0 } else { or1_len <- length(or1[,1]) }
	if(length(or3)==0) { or3_len <- 0 } else { or3_len <- length(or3[,1]) }

	imp_count <- ab1_len + ab3_len + or1_len + or3_len - or_redund1 - or_redund3

	if(file.info("and-all.sam")$size != 0) {
		system("cat and-all.sam or.sam > ALL.sam")
	} else { system("cp or.sam ALL.sam") }
	unlink("and-all.sam")

	imp_pairs <- imp_count / length(reads)
	imp_pairs <- imp_pairs * 100
	imp_pairs <- format(round(imp_pairs, digits=1), nsmall = 1)
	imp_pairs <- paste(imp_pairs, "%", sep="")

	# Perfect Pairs: 1 count per read that has been associated with a U1/A10 pair (no mismatches)
	if(file.info("ALL.sam")$size == 0) {
                return(cat("\n\n",
                        "\nLibrary Size: ", length(reads),
                        "\n\n\nNumber of Reads with:\n",
                        "\nPercent U1: ", num_u1, " (", perc_u1, ")",
                        "\nPercent A10: ", num_a10, " (", perc_a10, ")",
                        "\nPercent Pairs (1-2 mismatch):  0 ( 0% )",
                        "\nPercent Pairs (no mismatch): 0 ( 0% )",
                        "\n\n"))
	} else { 
		perf <- readLines("ALL.sam") 
		perf <- perf[grepl("XM:i:0",perf,fixed=TRUE)]
		write.table(perf, "Perfect.sam", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\n")
		if(file.info("Perfect.sam")$size == 0) {
			perf <- data.frame()
		} else { perf <- read.table("Perfect.sam", sep="\t", fill=TRUE, row.names=NULL) }
	
		if(file.info("Perfect.sam")$size==0) {
			perf_count <- 0; perf_pairs <- 0.0;
		} else if(length(perf[,1]) == 1) {
			if(grepl("Binary", perf)) { 
				perf_count <- 0; perf_pairs <- 0.0 
			} else { 
				perf_count <- 2
				perf_pairs <- perf_count / length(reads)
				perf_pairs <- perf_pairs * 100
				perf_pairs <- format(round(perf_pairs, digits=1), nsmall = 1)
			}
		} else {
			perf1 <- perf[!duplicated(as.character(perf[,1])),]
			perf3 <- perf[!duplicated(as.character(perf[,3])),]

			perf_count <- length(perf1[,1]) + length(perf3[,1])
			perf_pairs <- perf_count / length(reads)
			perf_pairs <- perf_pairs * 100
			perf_pairs <- format(round(perf_pairs, digits=1), nsmall = 1)
		}
		perf_pairs <- paste(perf_pairs, "%", sep="")

		if(verbose != "TRUE") { unlink("u1*.fa"); unlink("a10*.fa"); unlink("u1a10*.fa"); unlink("and1*.sam"); unlink("and2*.sam"); unlink("or*.sam"); unlink("ALL.sam"); unlink("*.bt2"); unlink("Perfect.sam")  }
		return(cat("\n\n",
			"\nLibrary Size: ", length(reads),
			"\n\n\nNumber of Reads with:\n",
			"\nPercent U1: ", num_u1, " (", perc_u1, ")",
			"\nPercent A10: ", num_a10, " (", perc_a10, ")",
			"\nPercent Pairs (1-2 mismatch): ", imp_count, " (", imp_pairs, ")",
			"\nPercent Pairs (no mismatch): ", perf_count, "(", perf_pairs, ")",
			"\n\n"))
	}
}

if(grepl(".fq", fq, fixed=TRUE) | grepl(".fastq", fq, fixed=TRUE)) { reads <- readFastq(fq) 
} else { reads <- readFasta(fq) }
reads <- reads[width(reads)>23]

## Create file of U1 AND A10 sRNA
u1a10 <- reads[substr(sread(reads),1,1)=="T" & substr(sread(reads),10,10)=="A"]
u1a10_f10 <- substr(sread(u1a10), 1, 10)
write.fasta(as.list(u1a10_f10), as.list(as.character(id(u1a10))), "u1a10_f10.fa")

## Create file of ONLY U1 sRNA
u1 <- reads[substr(sread(reads),1,1)=="T" & !(sread(reads) %in% sread(u1a10))]
u1_f10 <- substr(sread(u1), 1, 10)
write.fasta(as.list(u1_f10), as.list(as.character(id(u1))), "u1_f10.fa")

## Create file of ONLY A10 sRNA
a10 <- reads[substr(sread(reads),10,10)=="A" & !(sread(reads) %in% sread(u1a10))]
a10_f10 <- substr(sread(a10), 1, 10)
write.fasta(as.list(a10_f10), as.list(as.character(id(a10))), "a10_f10.fa")

## Find pairs of reads with both U1 AND A10
if(file.info("u1a10_f10.fa")$size != 0) {
	system("bowtie2-build u1a10_f10.fa u1a10_f10.fa")
	system("bowtie2 --no-unal --no-hd --no-sq -f -x u1a10_f10.fa -U u1_f10.fa -S and1.sam")
	system("bowtie2 --no-unal --no-hd --no-sq -f -x u1a10_f10.fa -U a10_f10.fa -S and2.sam")
} else {
	file.create("and1.sam")
	file.create("and2.sam")
}
## Find pairs for reads with only U1 OR A10
if(file.info("u1_f10.fa")$size != 0) {
	system("bowtie2-build u1_f10.fa u1_f10.fa")
	system("bowtie2 --no-unal --no-hd --no-sq -f -x u1_f10.fa -U a10_f10.fa -S or.sam")
} else {
	file.create("or.sam")
}

num_u1 <- length(u1) + length(u1a10)
perc_u1 <- num_u1 * 100
perc_u1 <- format(round(perc_u1 / length(reads), digits=1), nsmall = 1)
perc_u1 <- paste(perc_u1, "%", sep="")

num_a10 <- length(a10) + length(u1a10)
perc_a10 <- num_a10 * 100
perc_a10 <- format(round(perc_a10 / length(reads), digits=1), nsmall = 1)
perc_a10 <- paste(perc_a10, "%", sep="")

res <- printSumm(num_u1, num_a10, perc_u1, perc_a10, "and1.sam", "and2.sam", "or.sam", verbose)
