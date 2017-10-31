# Windowed analysis of TE distrubtions across the genomes of the core Armillaria species

# Aims:
# 1) Identify most meaningful interval to investigate clustering of TEs (test 10k, 50k, 100k)
# 2) Invesitigate whether there are species differences in how TEs are clustered among the genomes w.r.t. coding genes
# 3) Identify cutoffs for what constitues TE-rich regions to test for associations

# Read data
genome_files <-  grep(".genome",list.files("."),value=TRUE)
speciesNames <-  sub("(.*).genome", "\\1",genome_files)

for (species in speciesNames){
  TEcov_10k <- read.table(paste(species, "_TEcov_10k.bed",sep=""),h=F,col.names=c("chr","start","end","TEstart","TEend","seg_lenght","TEcov"))
  mergedCov_10k <- cbind(TEcov_10k,read.table(paste(species, "_genecov_10k.bed",sep=""),h=F,col.names=c("chr","start","end","gene_start","gene_end","seg_lenght","genecov")))
  mergedCov_10k <- subset(mergedCov_10k, select=c("chr","start","end","TEcov","genecov"))
  assign(paste(species,"_mergedCov_10k",sep=""),mergedCov_10k)
}

for (species in speciesNames){
  TEcov_50k <- read.table(paste(species, "_TEcov_50k.bed",sep=""),h=F,col.names=c("chr","start","end","TEstart","TEend","seg_lenght","TEcov"))
  mergedCov_50k <- cbind(TEcov_50k,read.table(paste(species, "_genecov_50k.bed",sep=""),h=F,col.names=c("chr","start","end","gene_start","gene_end","seg_lenght","genecov")))
  mergedCov_50k <- subset(mergedCov_50k, select=c("chr","start","end","TEcov","genecov"))
  assign(paste(species,"_mergedCov_50k",sep=""),mergedCov_50k)
}

for (species in speciesNames){
  TEcov_100k <- read.table(paste(species, "_TEcov_100k.bed",sep=""),h=F,col.names=c("chr","start","end","TEstart","TEend","seg_lenght","TEcov"))
  mergedCov_100k <- cbind(TEcov_100k,read.table(paste(species, "_genecov_100k.bed",sep=""),h=F,col.names=c("chr","start","end","gene_start","gene_end","seg_lenght","genecov")))
  mergedCov_100k <- subset(mergedCov_100k, select=c("chr","start","end","TEcov","genecov"))
  assign(paste(species,"_mergedCov_100k",sep=""),mergedCov_100k)
}

require(hexbin)
pdf(file="10k_coverage_plots.pdf")
plot(hexbin(Acepistipes_mergedCov_10k$genecov,Acepistipes_mergedCov_10k$TEcov),colorcut = seq(0,1,length.out=25), colramp = function(n) cols(25))
plot(hexbin(Agallica_mergedCov_10k$genecov,Agallica_mergedCov_10k$TEcov),colorcut = seq(0,1,length.out=25), colramp = function(n) cols(25))
plot(hexbin(Aostoyae_mergedCov_10k$genecov,Aostoyae_mergedCov_10k$TEcov),colorcut = seq(0,1,length.out=25), colramp = function(n) cols(25))
plot(hexbin(Asolidipes_mergedCov_10k$genecov,Asolidipes_mergedCov_10k$TEcov),colorcut = seq(0,1,length.out=25), colramp = function(n) cols(25))
plot(hexbin(Gnecrorhiza_mergedCov_10k$genecov,Gnecrorhiza_mergedCov_10k$TEcov),colorcut = seq(0,1,length.out=25), colramp = function(n) cols(25))
plot(hexbin(Ctorrendii_mergedCov_10k$genecov,Ctorrendii_mergedCov_10k$TEcov),colorcut = seq(0,1,length.out=25), colramp = function(n) cols(25))
dev.off()

pdf(file="50k_coverage_plots.pdf")
plot(hexbin(Acepistipes_mergedCov_50k$genecov,Acepistipes_mergedCov_50k$TEcov),colorcut = seq(0,1,length.out=25), colramp = function(n) cols(25))
plot(hexbin(Agallica_mergedCov_50k$genecov,Agallica_mergedCov_50k$TEcov),colorcut = seq(0,1,length.out=25), colramp = function(n) cols(25))
plot(hexbin(Aostoyae_mergedCov_50k$genecov,Aostoyae_mergedCov_50k$TEcov),colorcut = seq(0,1,length.out=25), colramp = function(n) cols(25))
plot(hexbin(Asolidipes_mergedCov_50k$genecov,Asolidipes_mergedCov_50k$TEcov),colorcut = seq(0,1,length.out=25), colramp = function(n) cols(25))
plot(hexbin(Gnecrorhiza_mergedCov_50k$genecov,Gnecrorhiza_mergedCov_50k$TEcov),colorcut = seq(0,1,length.out=25), colramp = function(n) cols(25))
plot(hexbin(Ctorrendii_mergedCov_50k$genecov,Ctorrendii_mergedCov_50k$TEcov),colorcut = seq(0,1,length.out=25), colramp = function(n) cols(25))
dev.off()

pdf(file="100k_coverage_plots.pdf")
plot(hexbin(Acepistipes_mergedCov_100k$genecov,Acepistipes_mergedCov_100k$TEcov),colorcut = seq(0,1,length.out=25), colramp = function(n) cols(25))
plot(hexbin(Agallica_mergedCov_100k$genecov,Agallica_mergedCov_100k$TEcov),colorcut = seq(0,1,length.out=25), colramp = function(n) cols(25))
plot(hexbin(Aostoyae_mergedCov_100k$genecov,Aostoyae_mergedCov_100k$TEcov),colorcut = seq(0,1,length.out=25), colramp = function(n) cols(25))
plot(hexbin(Asolidipes_mergedCov_100k$genecov,Asolidipes_mergedCov_100k$TEcov),colorcut = seq(0,1,length.out=25), colramp = function(n) cols(25))
plot(hexbin(Gnecrorhiza_mergedCov_100k$genecov,Gnecrorhiza_mergedCov_100k$TEcov),colorcut = seq(0,1,length.out=25), colramp = function(n) cols(25))
plot(hexbin(Ctorrendii_mergedCov_100k$genecov,Ctorrendii_mergedCov_100k$TEcov),colorcut = seq(0,1,length.out=25), colramp = function(n) cols(25))
dev.off()

write.table(Acepistipes_mergedCov_50k, file="Acepistipes_mergedCov_50k.tab",quote=F,sep="\t",row.names = F)
write.table(Agallica_mergedCov_50k, file="Agallica_mergedCov_50k.tab",quote=F,sep="\t", row.names = F)
write.table(Aostoyae_mergedCov_50k, file="Aostoyae_mergedCov_50k.tab",quote=F,sep="\t", row.names = F)
write.table(Asolidipes_mergedCov_50k, file="Asolidipes_mergedCov_50k.tab",quote=F,sep="\t", row.names = F)
write.table(Gnecrorhiza_mergedCov_50k, file="Gnecrorhiza_mergedCov_50k.tab",quote=F,sep="\t", row.names = F)
write.table(Ctorrendii_mergedCov_50k, file="Ctorrendii_mergedCov_50k.tab",quote=F,sep="\t", row.names = F)