# This R script For the R script of the DADA2 analysis,
# the workflow of callahan_mcmurdie_rosen_han_johnson_holmes_2016 was used as a guideline.
# The following parameters were set as default for this workflow:
# 1) truncLen of 175 for both reads, this truncate reads after 175 bases. Reads shorter than this were discarded.
# 2) maxN of 0, after truncation reads with "N" were discarded.
# 3) maxEE of 2, after truncation, reads with higher than 2 expected errors were discarded.
# Expected errors were calculated from the nominal definition of the quality score, EE=$sum(10 ^{(-Q/10)})$.
# 4) truncQ of 2, truncates reads at the first instance of a quality score less than or equal to 2.
# 5) phix was set to TRUE, discarding reads that match against a phiX genome.
# 6) Compress was set to TRUE, making the output files gzipped.
# And 7) multithread was set to TRUE, filtering the input files in parallel.
# The samples with a sequence depth of more than twenty were kept, the others were discarded.
# All the samples that pass these steps were used to estimate the error rate based on the quality scores
# associated with all of those reads.
# From this error rate estimation, the program formed a sample inference with the core DADA2 algorithm.
# After this, the reads were merged together to create a sequence table. From this table chimeras were removed.
# The clean sequence table was used to assign taxonomy using the UNITE general FASTA release for Fungi database
# version 02/02/2019 \citep{unite}, with the use of BlastN taking the top result with a minimum e-value of 0.00001.


#remove all stored variables before start
rm(list = ls(all = TRUE))

#import the library dada2 and print the version for control.
#import the shortread library
library(dada2); packageVersion("dada2")
library(ShortRead)

# set the path to all the demultiplexed files and print the files.
path <- "../"
list.files(path)

# store all the forward fq and reverse fq files in the right variable.
fastqFs <- sort(list.files(path, pattern="-R1fw.fastq", full.names = TRUE))
fastqRs <- sort(list.files(path, pattern="-R2rv.fastq", full.names = TRUE))
fastqF2s <- sort(list.files(path, pattern="-R2fw.fastq", full.names = TRUE))
fastqR2s <- sort(list.files(path, pattern="-R1rv.fastq", full.names = TRUE))
# make a list of all the sample names, which is the beginning of the filename.
sample.names <- sapply(strsplit(basename(fastqFs), "-R"),`[`, 1)
sample.names2 <- sapply(strsplit(basename(fastqF2s), "-R"),'[', 1)
#print(sample.names)

# Plots the quality score of some forward and reverse sequences. From here the best trimming length can be determined.
plotQualityProfile(fastqFs[1:2])
plotQualityProfile(fastqRs[1:2])

plotQualityProfile(fastqF2s[1:2])
plotQualityProfile(fastqR2s[1:2])

# make a directory and all files before filtering them
filtF <- file.path(path, "filtered1", paste0(sample.names, "1_F1_filt1.fastq.gz"))
filtR <- file.path(path, "filtered1", paste0(sample.names, "1_R1_filt1.fastq.gz"))
filtF2 <- file.path(path, "filtered1", paste0(sample.names2, "2_F2_filt1.fastq.gz"))
filtR2 <- file.path(path, "filtered1", paste0(sample.names2, "2_R2_filt1.fastq.gz"))

# filter all the files by trimming them by the set length, also delete sequences with a N, a low Estimated Error.
control <- function(fastq) {
outy <- tryCatch(
    {
        for (i in fastq){
            readFastq(i)
}
    },
    error=function(cond) {
            message(paste("file corrupted: ", i))
            message("Here's the original error message:")
            message(cond)
            # Choose a return value in case of error
            return(NA)
        },
        finally={
            message(paste("Processed file: ", i))
            message("no problem found. ")
        }
    )
    return(outy)
}
y <- lapply(fastqRs, control)
w <- lapply(fastqFs, control)
v <- lapply(fastqF2s, control)
z <- lapply(fastqR2s, control)

out <- filterAndTrim(fastqFs, filtF, fastqRs, filtR, truncLen=c(175,175),
                maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                compress=TRUE, multithread=TRUE)

out2 <- filterAndTrim(fastqF2s, filtF2, fastqR2s, filtR2, truncLen=c(175,175),
                maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                compress=TRUE, multithread=TRUE)

# only keep the samples with more the 20 reads. store the sample names of the representive samples for Forward and Reverse
keep <- out[,"reads.out"] > 20
keep2 <- out2[,"reads.out"] >20

filtFk <- file.path(filtF)[keep]
filtFk2 <- file.path(filtF2)[keep2]
filtFs <- filtFk[!is.na(filtFk)]
filtFs2 <- filtFk2[!is.na(filtFk2)]
sample.namesF <- sapply(strsplit(basename(filtFs), "_F1"),'[', 1)
sample.namesF2 <- sapply(strsplit(basename(filtFs2), "_F2"), '[',1)

filtRs <- file.path(filtR)[keep]
filtRs2 <- file.path(filtR2)[keep2]
filtRs <- filtRs[!is.na(filtRs)]
filtRs2 <- filtRs2[!is.na(filtRs2)]
sample.namesR <- sapply(strsplit(basename(filtRs), "_R1"),'[', 1)
sample.namesR2 <- sapply(strsplit(basename(filtRs2), "_R2"),'[',1)

# Calculate the estimated error rate for the forward and reverse reads.
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
errF2 <- learnErrors(filtFs2, multithread=TRUE)
errR2 <- learnErrors(filtRs2, multithread=TRUE)

# polt the errors of the forward and reverse fastq files.
# This shows the probablity that a base actually is an other base,
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
plotErrors(errF2, nominalQ=TRUE)
plotErrors(errR2, nominalQ=TRUE)

# magic part :)
# the sample names in these objects are initially the file names of the samples,
# this sets them to the sample names for the rest of the workflow
derepFs <- derepFastq(filtFs)
derepFs2 <- derepFastq(filtFs2)
names(derepFs) <- sample.namesF
names(derepFs2) <- sample.namesF2
derepRs <- derepFastq(filtRs)
derepRs2 <- derepFastq(filtRs2)
names(derepRs) <- sample.namesR
names(derepRs2)<- sample.namesR2

# Do the dada2 algorithm on the made files with the estimated error rates.
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaFs2 <- dada(derepFs2, err=errF2, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaRs2 <- dada(derepRs2, err=errR2, multithread=TRUE)
#dadaFs[[1]]
#dadaRs[[1]]

# merge the pairs from forward and reverse reads together to get a representive sequence.
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs,)
mergers2 <- mergePairs(dadaFs2, derepFs2, dadaRs2, derepRs2)
# Inspect the merger data.frame from the first sample and the class, length and names as a check up.
# head(mergers[[1]])
# class(mergers)
# length(mergers)
# names(mergers) # the names() function gives us the name of each element of the list

# make a sequence table of the merges.
seqtab <- makeSequenceTable(mergers)
seqtab2 <- makeSequenceTable(mergers2)
#class(seqtab)
#dim(seqtab)
#View(t(seqtab))
# summary.matrix(t(seqtab))
saveRDS(seqtab, "./seqtab.rds")
saveRDS(seqtab2, "./seqtab2.rds")

# Inspect distribution of sequence lengths
#table(nchar(getSequences(seqtab)))

# remove bimera's from the sequence table. save this table.
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
saveRDS(seqtab.nochim, "./ESVtable_NOCHIM.rds")
write.table(seqtab.nochim, file="ESVtable_NOCHIM.txt", row.names=TRUE, col.names=TRUE)
write.table(seqtab.nochim, file="ESVtable_NOCHIM.csv", row.names=TRUE, col.names=TRUE)

seqtab2.nochim2 <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab2.nochim2)
sum(seqtab2.nochim2)/sum(seqtab2)
saveRDS(seqtab2.nochim2, "./ESVtable_NOCHIM2.rds")
write.table(seqtab2.nochim2, file="ESVtable_NOCHIM2.txt", row.names=TRUE, col.names=TRUE)
write.table(seqtab2.nochim2, file="ESVtable_NOCHIM2.csv", row.names=TRUE, col.names=TRUE)

merged_seqtabnochim <- mergeSequenceTables(seqtab.nochim, seqtab2.nochim2, orderBy = "abundance")
saveRDS(merged_seqtabnochim, "./ESVtable_NOCHIM_MER.rds")
write.table(merged_seqtabnochim, file="ESVtable_NOCHIM_MER.txt", row.names=TRUE, col.names=TRUE)
write.table(merged_seqtabnochim, file="ESVtable_NOCHIM_MER.csv", row.names=TRUE, col.names=TRUE)

# assign taxonomy to the no chimera table, form the last known fungi database. save the data in differend formats.
axa <- assignTaxonomy(merged_seqtabnochim, "/cluster/home/veerle/test/sh_general_release_dynamic_02.02.2019.fasta", multithread=TRUE)
saveRDS(axa, "./tax_TableFinal.rds")
write.table(axa, file="ESVtable_final.txt", row.names=TRUE, col.names=TRUE)
write.table(axa, file="ESVtable_final.csv", row.names=TRUE, col.names=TRUE)

