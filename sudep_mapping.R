setwd("~/Desktop/JohnS/probextant")
library(abind)
library(doBy)
library(DOQTL)

#Collapse into 8 state. Only need to do once!!
setwd("~/Desktop/JohnS/probextant/MRCA")
sam<-scan("MRCAlist.csv", what="list", quote=NULL)
for (i in sam) {
  hap36<-read.csv(i)
  hap8<-hap36[,1:3]
  AA<-hap36$AA+(hap36$AB/2)+(hap36$AC/2)+(hap36$AD/2)+(hap36$AE/2)+(hap36$AF/2)+(hap36$AG/2)+(hap36$AH/2)
  BB<-hap36$BB+(hap36$AB/2)+(hap36$BC/2)+(hap36$BD/2)+(hap36$BE/2)+(hap36$BF/2)+(hap36$BG/2)+(hap36$BH/2)
  CC<-hap36$CC+(hap36$AC/2)+(hap36$BC/2)+(hap36$CD/2)+(hap36$CE/2)+(hap36$CF/2)+(hap36$CG/2)+(hap36$CH/2)
  DD<-hap36$DD+(hap36$AD/2)+(hap36$BD/2)+(hap36$CD/2)+(hap36$DE/2)+(hap36$DF/2)+(hap36$DG/2)+(hap36$DH/2)
  EE<-hap36$EE+(hap36$AE/2)+(hap36$BE/2)+(hap36$CE/2)+(hap36$DE/2)+(hap36$EF/2)+(hap36$EG/2)+(hap36$EH/2)
  FF<-hap36$FF+(hap36$AF/2)+(hap36$BF/2)+(hap36$CF/2)+(hap36$DF/2)+(hap36$EF/2)+(hap36$FG/2)+(hap36$FH/2)
  GG<-hap36$GG+(hap36$AG/2)+(hap36$BG/2)+(hap36$CG/2)+(hap36$DG/2)+(hap36$EG/2)+(hap36$FG/2)+(hap36$GH/2)
  HH<-hap36$HH+(hap36$AH/2)+(hap36$BH/2)+(hap36$CH/2)+(hap36$DH/2)+(hap36$EH/2)+(hap36$FH/2)+(hap36$GH/2)
  hap8<-cbind(hap8, AA, BB, CC, DD, EE, FF, GG, HH)
  write.csv(hap8, file=paste("8", i, sep="_"), row.names=FALSE)
  rm(hap8)
  rm(hap36)
  
}


#Creating the prob array. Only need to do once!
Sam.List<-scan("MRCAlist2.csv", what="list", quote=NULL)

dat1<-read.csv("8_CC001-Uncb38V01.csv")
Snp.Info<-dat1[,1:3]
g.Info<-dat1[,4:11]
gT1<-t(g.Info)

dat2<-read.csv("8_CC002-Uncb38V01.csv")
g.Info<-dat2[,4:11]
gT2<-t(g.Info)

y<-abind(gT1, gT2, along=0)

for( i in 3:length(Sam.List))
{
  a<-Sam.List[[i]]
  g.Info<-read.csv(a)
  g.Info<-g.Info[,4:11]
  gT<-t(g.Info)
  y<-abind(y, gT, along=1)
}

#name the 3 sides of the array
markers<-list(Snp.Info$marker)
haps<-list('A','B','C','D','E','F','G','H')
dimnames(y)[[1]] <- unlist(Sam.List)
dimnames(y)[[2]] <- unlist(haps)
dimnames(y)[[3]] <- unlist(markers)

#check to see if number of markers, haps, and samples are correct
str(y)
#is the haplotype frequency correct
freqsy <- apply(y, c(1,2), mean)
colMeans(freqsy)



library(ggplot2)
library(plyr)
devtools::load_all("~/Desktop/argyle")
devtools::load_all("~/Desktop/mouser/")

###Sanity check. Does the frequency make sense? For b38 it does not at chr5 and chr13!
## load SNP positions
colnames(Snp.Info) <- c("marker", "chr", "pos")

## take mean across samples to get allele freq at each marker
freqs <- colMeans(y, dims = 1)

## turn it into a dataframe for plotting
freqs2 <- reshape2::melt(freqs, value.name = "freq")
colnames(freqs2)[1:2] <- c("strain","marker")

## add SNP positions
freqs2 <- merge(freqs2, Snp.Info, all.x = TRUE)

## make the plot
ggmanhattan(freqs2) +
  geom_line(aes(y = freq, colour = strain, group = strain:chr)) +
  geom_hline(yintercept = 0.125, lty = "dashed", colour = "darkgrey") +
  scale_y_continuous("allele frequency\n") +
  scale_color_CC() +
  facet_grid(strain ~ .) +
  theme_classic() +
  theme(axis.title.x = element_blank())


##traitmapping with CC epi data
sudep <- read.csv("~/Desktop/JohnS/probextant/sudep.csv")

GIGA_snps <- read.csv("~/Desktop/JohnS/probextant/MRCA/map.csv")
rownames(sudep) <- as.character(sudep$samples)

#doqtl
covar = data.frame(sex = as.numeric(sudep$sex == "M"), row.names = as.character(sudep$samples))
sudepqtl = scanone(pheno = sudep, pheno.col = "sudep", probs = y, addcovar = covar, snps = GIGA_snps)
plot(sudepqtl, main = "sudep")

coefplot(sudepqtl, chr=17)

sudep_perm <- scanone.perm(pheno = sudep, pheno.col = "sudep", probs = y, addcovar = covar, snps = GIGA_snps, nperm = 20)
plot(sudep_perm, main = "sudep_perm")



##trying rqtl2
library("qtl2", lib.loc="/Library/Frameworks/R.framework/Versions/3.3/Resources/library")
library("qtl2convert", lib.loc="/Library/Frameworks/R.framework/Versions/3.3/Resources/library")

#convert from the doqtl version to qtl2 version. 
probsy <- probs_doqtl_to_qtl2(probs = y, map = GIGA_snps, chr_column="chromosome", pos_column="position", marker_column="marker")

#reordering the chromosomes, otherwise it goes 1,10,11,..
probsy <- subset(probsy, chr=c(1:19, "X")) 

#seems to always need a covariate.
addcovar <- subset(sudep, select=c("sex"))
rownames(addcovar) <- sudep$sample

#preparing the phenotype. 
pheno <- as.matrix(sudep$sudep)
rownames(pheno) <- sudep$sample

#not always necessary to have kinship. Various types to chose, see documentation if interested in learning more.
kinshipadd <- calc_kinship(probsy, type = c("overall"), use_allele_probs=FALSE)

#this is the qtl mapping function.
out <- scan1(probsy, pheno, kinship = kinshipadd, model="binary")

#Need to change map to rqtl2 friendly version for plotting qtl and allelic effects
map <- map_df_to_list(GIGA_snps, chr_column = "chromosome", pos_column = "position", marker_column = "marker", Xchr = c("X"))

#ploting out QTL
plot(out, map, lodcolumn=1, col="slateblue")

#scan1blup is the way to show allelic effects. Specify what chromosome.
blup2 <- scan1blup(probsy[,"17"], pheno, kinshipadd)
plot_coefCC(blup2, map["17"])


#an alternative way to plot. Sometimes doesn't show all allelic effects is frequency is skewed
# ymx <- maxlod(out)
# coef_c2 <- scan1coef(probsy[,"17"], pheno, kinshipadd)
# plot_coefCC(coef_c2, map["17"], bgcolor="gray95")
# legend("topleft", col=CCcolors, names(CCcolors), ncol=2, lwd=2, bg="gray95")