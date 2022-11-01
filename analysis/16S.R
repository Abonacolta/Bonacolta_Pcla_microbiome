```{bash}
ls *_L001_R1_001.fastq.gz | cut -f 1-2 -d "_" > samples

for sample in $(cat samples)
do

    echo "On sample: $sample"

    cutadapt -a ^ACGCGHNRAACCTTACC...TTGYACWCACYGCCCGT \
    -A ^ACGGGCRGTGWGTRCAA...GGTAAGGTTYNDCGCGT \
    -m 200 -M 300 --discard-untrimmed \
    -o ${sample}_L001_R1_001_trimmed.fastq.gz -p ${sample}_L001_R2_001_trimmed.fastq.gz \
    ${sample}_L001_R1_001.fastq.gz ${sample}_L001_R2_001.fastq.gz \
    >> cutadapt_primer_trimming_stats.txt 2>&1

done
```


library(dada2); packageVersion("dada2")
path <- "/projects/scratch/marinemicrobes/pcla/16S/data" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)
fnFs <- sort(list.files(path, pattern="_L001_R1_001_trimmed.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_L001_R2_001_trimmed.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_L001"), `[`, 1)
plotQualityProfile(fnFs[1:3])
plotQualityProfile(fnRs[1:3])
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=10, truncLen=c(250,210),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE, matchIDs=TRUE)
head(out)
errF <- learnErrors(filtFs, multithread=14)
errR <- learnErrors(filtRs, multithread=14)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
sam.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
names(derepFs) <- sam.names
names(derepRs) <- sam.names
ddFs <- dada(derepFs, err=NULL, selfConsist=TRUE)
ddRs <- dada(derepRs, err=NULL, selfConsist=TRUE)
plotErrors(ddFs)
plotErrors(ddRs)
dadaFs <- dada(derepFs, err=ddFs[[1]]$err_out, pool=TRUE, multithread=14)
dadaRs <- dada(derepRs, err=ddRs[[1]]$err_out, pool=TRUE, multithread=14)
dadaFs[[1]]
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
seqtab.all <- makeSequenceTable(mergers)
seqtab <- removeBimeraDenovo(seqtab.all, minFoldParentOverAbundance=16, method="consensus", multi=TRUE)
dim(seqtab)
table(nchar(getSequences(seqtab))) # Inspect distribution of sequence lengths
seqtab.chim <- getSequences(seqtab.all)[!getSequences(seqtab.all) %in% getSequences(seqtab)]
ref_fasta <- "/projects/scratch/marinemicrobes/silva_nr99_v138_train_set.fa.gz"
taxa <- assignTaxonomy(seqtab, refFasta=ref_fasta, multithread=14)
colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
taxa <- addSpecies(taxa, "/projects/scratch/marinemicrobes/silva_species_assignment_v138.fa.gz")
taxa.print <- taxa
save.image("/projects/scratch/marinemicrobes/pcla/16S/analysis/pcal_16S_peg.RData")
