library(dada2) #devtools::install_github("benjjneb/dada2")
library(phyloseq) #BiocManager::install("phyloseq")

#Create a path to where the appropriate merged fastq files 
path <- "./Merged_Files"
#Create list of all fastq files in the path 
fnFs <- sort(list.files(path, pattern=".fastq", full.names = TRUE))
#Split to only include the ID; 10708.assembled.fastq.gz -> 10708
sample.names <- sapply(strsplit(basename(fnFs), ".assembled"), `[`, 1)

#Create a list of the desired output using the IDs; 10708_filt.fastq.gz
filtFs <- file.path(path, "NoTrun_Q5_EE0.5", paste0(sample.names, "_filt.fastq.gz"))
names(filtFs) <- sample.names

#truncLen=0; No trun length since we've already merged 
#MaxN = 0; Remove sequences with any N 
#maxEE = 0.5, Remove sequences with expected error greater than 0.5
#truncQ=5, Truncate reads at the first instance of a quality score less than or equal to 5
#rm.phix = TRUE, discard reads that match against the phiX genome
#compress = TRUE, gzip output file 
#multithread = TRUE, Allow parallelization to increase speed 
out <- filterAndTrim(fnFs, filtFs, truncLen=0, 
                     maxN=0, maxEE=0.5, truncQ=5, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out) #reads.in reads.out

#Visoulixe quality of reads 
plotQualityProfile(filtFs[1:4])

#Create list of filtered fastqs 
filtFs <- sort(list.files(paste(path, "/NoTrun_Q5_EE0.5", sep = ""), pattern="filt.fastq.gz", full.names = TRUE))
#Name the filtered files by the ID
names(filtFs) <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)

#Learn Error rate
errF <- learnErrors(filtFs, multithread=TRUE)
#124215057 total bases in 303192 reads from 7 samples will be used for learning the error rates.
plotErrors(errF, nominalQ=TRUE)

#https://benjjneb.github.io/dada2/pseudo.html#Pseudo-pooling
#pooling information across samples can increase sensitivity to sequence variants 
#that may be present at very low frequencies in multiple samples.
dadaFs_pool <- dada(filtFs, err=errF, multithread=TRUE, pool = "pseudo")
seqtab_pooled <- makeSequenceTable(dadaFs_pool)

#Remove chimeras 
seqtab.nochim_pooled <- removeBimeraDenovo(seqtab_pooled, 
                                           method="consensus", multithread=TRUE, verbose=TRUE)
#Identified 93995 bimeras out of 152083 input sequences.

#Determine % of non-chimera samples 
sum(seqtab.nochim_pooled)/sum(seqtab_pooled)
#[1] 0.8468572

seqtab.nochim_pooled = readRDS("~/Desktop/TEDDY_CeD/seqtab.1.nochim.RDS")

#'McLaren, Michael R., & Callahan, Benjamin J. (2021). Silva 138.1 prokaryotic SSU taxonomic training data formatted for DADA2 [Data set]. Zenodo. https://doi.org/10.5281/zenodo.4587955

#Assign Taxonomy [Excluding species] using silva: https://zenodo.org/record/3731176#.Ys3JdOzMKAk
taxa_pool <- assignTaxonomy(seqtab.nochim_pooled, "~/Downloads/rdp_train_set_18.fa.gz", multithread=TRUE)

#Not enough memory to run at once, so split into partial tables 
taxa_1 = taxa_pool[1:10000,]
taxa_2 = taxa_pool[10001:20000,]
taxa_3 = taxa_pool[20001:30000,]
taxa_4 = taxa_pool[30001:40000,]
taxa_5 = taxa_pool[40001:nrow(taxa_pool),]

#Assign species to the Silva: https://zenodo.org/record/3731176#.Ys3JdOzMKAk
taxa_1 <- addSpecies(taxa_1, "~/Downloads/rdp_species_assignment_18.fa.gz")
taxa_2 <- addSpecies(taxa_2, "~/Downloads/rdp_species_assignment_18.fa.gz")
taxa_3 <- addSpecies(taxa_3, "~/Downloads/rdp_species_assignment_18.fa.gz")
taxa_4 <- addSpecies(taxa_4, "~/Downloads/rdp_species_assignment_18.fa.gz")
taxa_5 <- addSpecies(taxa_5, "~/Downloads/rdp_species_assignment_18.fa.gz")

#Bind all the subsets tofether 
taxa_total = rbind(taxa_1, taxa_2)
taxa_total = rbind(taxa_total, taxa_3)
taxa_total = rbind(taxa_total, taxa_4)
taxa_total = rbind(taxa_total, taxa_5)

#Save the dada2 output
saveRDS(seqtab_pooled, "RDP_seqtab.1.RDS")
saveRDS(seqtab.nochim_pooled, "seqtab.1.nochim.RDS")
saveRDS(taxa_total, "RDP_taxa.1.RDS")

#Create and Save Phyloseq Object 
phyloseq_final = phyloseq(otu_table(seqtab.nochim_pooled, taxa_are_rows = FALSE),
                          tax_table(taxa_total))
phyloseq_final
saveRDS(phyloseq_final, "RDP_phyloseq.1.RDS")
summary(sample_sums(phyloseq_final))
