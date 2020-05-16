# This is inspired by https://benjjneb.github.io/dada2/tutorial_1_8.html
# The script required paired end data
# Currently it is tested with mothur miseq example dataset
# Script can be run with Rscript dada2_16S_worflow.R --input inputfolder --output --outputfolder

library(dada2); packageVersion("dada2")
library(argparser, quietly=TRUE)

p <- arg_parser("DADA2 in R")

p = add_argument(p, "--input", help="Enter input folder", type='character')
p = add_argument(p, "--output", help="Enter output folder", type='character')
p = add_argument(p, "--reference", help="Enter reference folder", type='character')
p = add_argument(p, "--plots", help="Enter plot folder path", type='character')

# Parse the command line arguments
argv <- parse_args(p)

path = argv$input
outdir = argv$output
reference = argv$reference
p_outdir = argv$plots

#print(path)
#print(outdir)
#print(reference)

files <- list.files(path)
#fastqs <- files[grepl(".fastq.gz$", files)] # gz
fastq_files <- files[grepl(".fastq", files)]


######################################################################################################
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
######################################################################################################
#fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE)) # for gz
#fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE)) # for gz
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)  #change based on sample names


####################################################################
# Visualizing the quality profiles of the forward and reverse reads
####################################################################
pdf(paste0(p_outdir,"dada2_plotQualityProfile.pdf"), onefile=T)
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
dev.off()

###################################
## Filter and trimming on raw reads
###################################
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))




#############################################################
# Perform filtering and trimming on forward and reverse reads
##############################################################
# place under filtered subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# Modify and add new parameter based on your data.
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                                   maxN=0, maxEE=c(2,2), truncQ=5, rm.phix=TRUE,
                                                 compress=TRUE, multithread=TRUE)


## Examine quality profiles of filtered reads
pdf(paste0(p_outdir, "QualityProfile.filt_plot.pdf"),onefile=T)
plotQualityProfile(filtFs[1:2])
plotQualityProfile(filtRs[1:2])
dev.off()


#####################################################
# Learn the Error Rates for forward and reverse reads
#####################################################

errF <- learnErrors(filtFs, multithread=TRUE)

errR <- learnErrors(filtRs, multithread=TRUE)

## Plot estimated error as sanity check
pdf(paste0(p_outdir, "ErrorsRates_F.pdf"), onefile=T)
plotErrors(errF, nominalQ=TRUE)
dev.off()

pdf(paste0(p_outdir,"ErrorsRates_R.pdf"), onefile=T)
plotErrors(errR, nominalQ=TRUE)
dev.off()


################
# Dereplication
################
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


###################
# Sample Inference
###################
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)


####################
# Merge paired reads
####################
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])


##########################
# Construct sequence table
###########################

seqtab <- makeSequenceTable(mergers)
## Get dimensions
#dim(seqtab)

## Inspect distribution of sequence lengths
#table(nchar(getSequences(seqtab)))


##################
# Remove chimeras
##################
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#################
#Save into Rds
#################
saveRDS(seqtab, paste0(outdir,"/seqtab.Rds"))
saveRDS(seqtab.nochim, paste0(outdir,"/seqtab.nochim.Rds"))


#################################
# rack reads through the pipeline
#################################
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track, paste0(outdir,"/track_reads.txt"),sep="\t",quote=FALSE)

#################
# Assign taxonomy
#################
#IMP : download taxonomy file wget https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz?download=1
taxa <- assignTaxonomy(seqtab.nochim, paste0(reference,"/silva_nr_v132_train_set.fa"), multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

write.table(taxa, file = paste0(outdir,"/taxa.tsv"), quote=FALSE)

