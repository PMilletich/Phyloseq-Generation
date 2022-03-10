library(dada2) #devtools::install_github("benjjneb/dada2")
library(phyloseq) #BiocManager::install("phyloseq")

library(knitr)
library(ShortRead)
library(BiocManager)
library(installr)
library(devtools)
library(microseq)

#Create path to assembled Fastqs from previous step 
path <- "/Volumes/Seagate_Turpin/Forward_Reverse/Combined/Merge_Filter"
#Create list of all fastq files 
fnFs <- sort(list.files(path, pattern=".fastq", full.names = TRUE))
#Create a list of the sample_names by splitting the string 
sample.names <- sapply(strsplit(basename(fnFs), ".assembled"), `[`, 1)
#Create a pathway for filtered fastq files 
filtFs <- file.path(path, "NoTrun_Q10_min400", paste0(sample.names, "_filt.fastq.gz"))
names(filtFs) <- sample.names

#https://rdrr.io/github/benjjneb/dada2/man/filterAndTrim.html
out <- filterAndTrim(fnFs, filtFs, truncLen=0, minLen = 400,
                    maxN=0, maxEE=2, truncQ=10, rm.phix=TRUE,
                    compress=TRUE, multithread=TRUE)
head(out)

plotQualityProfile(filtFs[1:4])

filtFs <- sort(list.files(paste(path, "/NoTrun_Q10_min400", sep = ""), pattern="filt.fastq.gz", full.names = TRUE))
names(filtFs) <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)

errF <- learnErrors(filtFs, multithread=TRUE)
#124215057 total bases in 303192 reads from 7 samples will be used for learning the error rates.
plotErrors(errF, nominalQ=TRUE)

dadaFs_pool <- dada(filtFs, err=errF, multithread=TRUE, pool = "pseudo")
seqtab_pooled <- makeSequenceTable(dadaFs_pool)
seqtab.nochim_pooled <- removeBimeraDenovo(seqtab_pooled, 
                                           method="consensus", multithread=TRUE, verbose=TRUE)
#Identified 7385 bimeras out of 17536 input sequences.
sum(seqtab.nochim_pooled)/sum(seqtab_pooled)
#[1] 0.9937659
taxa_pool <- assignTaxonomy(seqtab.nochim_pooled, "/Volumes/Seagate_Turpin/Trish/silva_nr_v138_train_set.fa", multithread=TRUE)
dim(taxa_pool)
#[1] 10151     6

#Not enough memory to run at once
taxa_1 = taxa_pool[1:2000,]
taxa_2 = taxa_pool[2001:4000,]
taxa_3 = taxa_pool[4001:6000,]
taxa_4 = taxa_pool[6001:9000,]
taxa_5 = taxa_pool[9001:nrow(taxa_pool),]

taxa_1 <- addSpecies(taxa_1, "/Volumes/Seagate_Turpin/Trish/silva_species_assignment_v138.fa")
taxa_2 <- addSpecies(taxa_2, "/Volumes/Seagate_Turpin/Trish/silva_species_assignment_v138.fa")
taxa_3 <- addSpecies(taxa_3, "/Volumes/Seagate_Turpin/Trish/silva_species_assignment_v138.fa")
taxa_4 <- addSpecies(taxa_4, "/Volumes/Seagate_Turpin/Trish/silva_species_assignment_v138.fa")
taxa_5 <- addSpecies(taxa_5, "/Volumes/Seagate_Turpin/Trish/silva_species_assignment_v138.fa")

taxa_total = rbind(taxa_1, taxa_2)
taxa_total = rbind(taxa_total, taxa_3)
taxa_total = rbind(taxa_total, taxa_4)
taxa_total = rbind(taxa_total, taxa_5)

saveRDS(seqtab_pooled, "/Volumes/Seagate_Turpin/Forward_Reverse/Combined/Merge_Filter/seqtab.13.RDS")
saveRDS(seqtab.nochim_pooled, "/Volumes/Seagate_Turpin/Forward_Reverse/Combined/Merge_Filter/seqtab.13.nochim.RDS")
saveRDS(taxa_total, "/Volumes/Seagate_Turpin/Forward_Reverse/Combined/Merge_Filter/taxa.13.RDS")


phyloseq_final = phyloseq(otu_table(seqtab.nochim_pooled, taxa_are_rows = FALSE),
                          tax_table(taxa_total))
phyloseq_final
saveRDS(phyloseq_final, "/Volumes/Seagate_Turpin/Autoimmune_Analysis/phyloseq.13.RDS")
table(sample_sums(phyloseq_final))
