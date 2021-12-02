library(tidyverse)
library(vcfR)
library(qvalue)
library(OutFLANK)
library(glue)
library(runner)
library(SNPRelate)
library(SeqArray)

##### files and information required
vcf_file <- "D:/Victoria University/Trevally/outflank_test/Allison/test_vcf.vcf.gz"
bim_file <- "D:/Victoria University/Trevally/outflank_test/Allison/binary_fileset_new1.bim"
fam_file <- "D:/Victoria University/Trevally/outflank_test/Allison/binary_fileset.fam"
chr_file <- "D:/Victoria University/Trevally/outflank_test/Allison/chromosome_G_aus_testContig_new.tsv"
out_dir  <- "D:/Victoria University/Trevally/outflank_test/Allison/"
out_ext  <- "TRE_tot_raw"
chr_ext  <- "Contig"

##### set parameters
qval           <- 0.05
LTrim          <- 0.1
RTrim          <- 0.1
Hmin           <- 0.1
sliding_window <- 500000

#extension to genotype files
dir.create(out_dir,recursive = T) #dir for Outflank
outflank_file_ext  <- paste0(out_dir,"/",out_ext)
#load data
bim <- read_tsv(bim_file,col_names = F, 
                col_types = cols(
                  X1 = col_character(),
                  X2 = col_character(),
                  X3 = col_integer(),
                  X4 = col_integer(),
                  X5 = col_character(),
                  X6 = col_character()
                )) 
colnames(bim) <- c("CHR","LOC","PAR","POS","REF","ALT")
bim$CHR <-  as.numeric(str_remove(bim$CHR,chr_ext))

fam <- read_delim(fam_file,col_names = F,delim = " ") ; colnames(fam) <- c("POP","IND","PAT","MAT","SEX","PHE") 
#sample_info <- read_tsv(glue("./data/{dataset}/{dataset}_sample_info.tsv"))
#load chromosome information
Chr_info     <- read_tsv(chr_file)
Chr_info$CHR <-  as.numeric(str_remove(Chr_info$CHR,chr_ext))
Chr_info <- Chr_info %>% mutate(tot= cumsum(Length)-Length) 
pops   <- as.character(unique(fam$POP))
N_pops <- length(pops)

#N_chr <- nrow(Chr_info)

#add genome information to bim file
bim <- left_join(bim,Chr_info) %>%
  # Add a cumulative position of each SNP
  arrange(CHR, POS) %>%
  mutate(BPcum  = tot+POS)

### read VCF file
vcfR <- vcfR::read.vcfR(vcf_file)
#convert vcf to outFLANK format
geno <- extract.gt(vcfR) # Character matrix containing the genotypes
chromosome <- getCHROM(vcfR) # linkage group
position   <- getPOS(vcfR) # Positions in bp
ID         <- paste0(chromosome,"_",position)
G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))
G[geno %in% c("0/0", "0|0")]                <- 0
G[geno %in% c("0/1", "1/0", "1|0", "0|1")]  <- 1
G[geno %in% c("1/1", "1|1")]                <- 2
G[geno %in% c(NA)]                          <- 9
snp <- t(G)
### estimate FST values
if(file.exists(glue::glue("{outflank_file_ext}_outliers.tsv"))){
  
    fst <- read_tsv(glue::glue("{outflank_file_ext}_FST_estimates.tsv"))  
    
} else {
    fst <- MakeDiploidFSTMat(snp, locusNames = ID, popNames = fam$POP)
    fst <- fst %>% mutate(POS       = as.integer(str_extract(LocusName,"[:digit:]*$"))
                         ,FST_ratio = FST/FSTNoCorr
                         ,LOC       = str_replace(LocusName,"_",":"))
    ### write FST estimates ### 
    write_tsv(fst,glue::glue("{outflank_file_ext}_FST_estimates.tsv"))
    ### sumary plot SNPs for LG
    png(glue::glue("{outflank_file_ext}_FST_correction.png"),res = 300, units = "in", width = 12, height = 5)
    par(mfrow=c(1,2))
    plot(fst$FST, fst$FSTNoCorr, 
         xlab="FST", ylab="uncorrected FST",
         pch=20, col="blue") ; abline(0,1)
    hist(fst$FSTNoCorr,xlim=c(0,1)  ,breaks=seq(0,1, by=0.001))
    dev.off()
} 





### identify outliers
fst_trim <-OutFLANK(fst, NumberOfSamples = N_pops, qthreshold = qval, Hmin = Hmin,
                    LeftTrimFraction = LTrim, RightTrimFraction = RTrim)
fst_trim$results <- left_join(fst_trim$results,bim) 
fst_trim$results <- fst_trim$results %>% dplyr::filter(!base::is.na(FST))
#test <- fst_trim$results
### select outliers
candidates     <- which(fst_trim$results$OutlierFlag==TRUE)
FST_candidates <- fst_trim$results[candidates,]
### thin outliers to remove 
runner             <-  runner(x    = FST_candidates
                              ,idx = "BPcum"  
                              ,f   = function(x){x[["LocusName"]][[which.max(x[["FST"]])]] }     
                              ,k   = sliding_window
                              ,at  = seq(from =1, to = max(FST_candidates$BPcum), by = sliding_window)
                              ,na_pad = FALSE)
FST_select      <- FST_candidates %>% dplyr::filter(LocusName %in% runner)

if(nrow(FST_select) != 0 ){
  outliers <- dplyr::select(FST_select,LocusName)
  write_tsv(outliers, glue::glue("{outflank_file_ext}_outliers.tsv"))
}
### create summary plot
png(glue::glue("{outflank_file_ext}_outfliank_output.png"),res = 300, units = "in", width = 16, height = 10)
par(mfrow=c(2,2))
plot(fst_trim$results$He, fst_trim$results$FST, pch=20, col="grey")
points(FST_candidates$He, FST_candidates$FST, pch=21, col="blue")
points(FST_select$He, FST_select$FST, pch=21, col="green")

plot(fst_trim$results$BPcum, fst_trim$results$FST, pch=20, col="grey")
points(FST_candidates$BPcum, FST_candidates$FST, pch=21, col="blue")
points(FST_select$BPcum, FST_select$FST, pch=21, col="green")

OutFLANKResultsPlotter(fst_trim, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)

hist(fst_trim$results$pvaluesRightTail)
dev.off()


POS_df <- FST_select %>% dplyr::select(LocusName) %>% tidyr::separate(LocusName,sep="_",remove=TRUE, into = c("LG","POS"))
N_all  <- nrow(FST_candidates)
N_SNPs <- nrow(POS_df)

#write files
write_tsv(POS_df        ,glue::glue("{outflank_file_ext}_{N_SNPs}_adaptive_sites_q{qval}.tsv"))
write_tsv(FST_candidates,glue::glue("{outflank_file_ext}_{N_all}_candidates_q{qval}.tsv"))
write_tsv(FST_select    ,glue::glue("{outflank_file_ext}_{N_SNPs}_outliers_q{qval}.tsv"))

#filter SNPs to new VCF
fst2 <- left_join(fst,FST_select)
SNP_for_VCF <- which(!is.na(fst2$qvalues))
New_vcfR <- vcfR[SNP_for_VCF,]
#test_vcfR <- New_vcfR@gt
colnames(New_vcfR@gt) <- c("FORMAT", fam$IND)
#export VCF and VCF.gz
write.vcf(New_vcfR, file=glue::glue("{outflank_file_ext}_adaptiveSNPs.vcf.gz"), mask=FALSE)
R.utils::gunzip(glue::glue("{outflank_file_ext}_adaptiveSNPs.vcf.gz"), remove = FALSE )

### convert vcf to gds and perform PCA on adaptive dataset
SeqArray::seqVCF2GDS(vcf.fn  = glue::glue("{outflank_file_ext}_adaptiveSNPs.vcf"),         out.fn    = glue::glue("{outflank_file_ext}_adaptiveSNPs_SeqArray.gds"))
#SNPRelate::snpgdsClose(gds_outflank)
SeqArray::seqGDS2SNP(gdsfile = glue::glue("{outflank_file_ext}_adaptiveSNPs_SeqArray.gds"),out.gdsfn = glue::glue("{outflank_file_ext}_adaptiveSNPs_SNPrelate.gds"))

### load SNPrelate gds file
#SNPRelate::snpgdsClose(gds_outflank)
gds_outflank <- SNPRelate::snpgdsOpen(glue::glue("{outflank_file_ext}_adaptiveSNPs_SNPrelate.gds"),readonly = F)
add.gdsn(gds_outflank,"pop.loc",fam$POP)

pca <- snpgdsPCA(gds_outflank, snp.id= c(1:N_SNPs), num.thread=2, autosome.only = F)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

tab <- data.frame(sample.id = pca$sample.id,
                  pop = factor(fam$POP)[match(pca$sample.id, fam$IND)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
#calculate genetic relationship matrix
GRM <- snpgdsGRM(gds_outflank, sample.id = fam$IND, method = "Eigenstrat", autosome.only = FALSE)
#perform abd plot PCA
png(filename =  glue::glue("{outflank_file_ext}_PCA.png"), width = 8, height = 4, units = "in", res = 300) 
    ggplot(tab,aes(EV1,EV2, color = pop))+
    geom_point()+
    xlab(paste0("eigenvector1 - ",round(pc.percent[1],2),"%"))+
    ylab(paste0("eigenvector2 - ",round(pc.percent[2],2),"%"))
dev.off()


