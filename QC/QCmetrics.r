
#####################
## Get Coverage plots  
# options(width=40)
# 
# library(CoverageView)

# setwd("~/Documents/Master_thesis/downstream_analysis")
# bigWigFile<-system.file("./bw", "IN_MN_2i_ALL.bw", package="CoverageView") 
# bamFile<-system.file("./bam", "IN_MN_2i_ALL.mm9.bam", package="CoverageView")
# 
# G4_bw<-CoverageBigWigFile("./bw/IN_MN_2i_ALL.bw", reads_mapped=54283037)
# G4_bam<-CoverageBamFile("./bam/IN_MN_2i_ALL.mm9.bam", reads_mapped=54283037)
# 
# G4_cov<-cov.matrix(G4_bw, coordfile="./MACS2/MN_2i_ALL_summits.bed", extend=1000, num_cores=4, bin_width=10)
# 
# draw.profile(G4_cov, ylab="avg coverage", outfile="./figures/IN_MN_2i_ALL.png", main="IN_MN_2i_ALL")
# draw.heatmap(G4_cov, outfile="./figures/IN_MN_2i_ALL.png")

######################
## Get strand cross-correlation plots 
# library(spp)
# library(csaw)
# library(imbforge)
# library(encodeChIPqc)
# library(stringr)
# 
# setwd("~/Documents/Master_thesis/downstream_analysis")
# # Specify path to treatment and input file"
# chipFile<-"./bam/G4_MN_2i_ALL.mm9.bam"
# inputFile<-"./bam/IN_MN_2i_ALL.mm9.bam"
# 
# chip<-read.bam.tags(chipFile)
# 
# input<-read.bam.tags(inputFile)
# 
# binding.characteristics<-get.binding.characteristics(chip,srange=c(50,500),bin=5,accept.all.tags=T)
# bc <- calculateBindingCharacteristics(chip, input, read.len=45)
# 
# print(paste("binding peak separation distance =",bc$peak$x))
# pp_B<-phantomPeak(binding.characteristics$chip, binding.characteristics$input,binding.characteristics$binding.characteristics)
# pp<-phantomPeak(bc$chip,bc$input,bc$bc)
# n<-500
# par(mfrow=c(2,2))
# corr <- correlateReads(chip, param=readParam(pe="both"), max.dist = n)
# 
# dev.off()
# plot(0:n, corr,  type='l', xlab="strand shift (bp)",ylab="cross-correlation")
# abline(v=bc$peak$x,lty=2,col=2)
# 
# 
############

## Get list of bamfiles in directory 
# dataDir<-"./bam"
# bamFiles<-dir(file.path(getwd(), dataDir), pattern="*.bam$", full.name=T)
# names(bamFiles)<-gsub(".bam","",basename(bamFiles))
# bamList<-BamFileList(bamFiles)
# 
# path(bamList)

## Check QC for a single bam file
# expFile<-"./bam"
# bamFile<-file.path(getwd(), expFile)
# exampleExp=ChIPQCsample(bamFile, peaks=mypeaks, blacklist=blacklist.mm9, annotation="mm9", chromosome="chr1")
# QCmetrics(exampleExp)
# RelativeCrossCoverage(exampleExp)
# plotCC(exampleExp)
# plotFrip(exampleExp)
# frip(exampleExp)
# regi(exampleExp)
# plotRegi(exampleExp)
# coveragehistogram(exampleExp)[1:10]
# plotCoverageHist(exampleExp)
# plotSSD(exampleExp)

################### 

## Get QCmetrics 
library(GenomicAlignments)
library(ChIPQC)
library(GenomicRanges)
library(rtracklayer)
library("TxDb.Mmusculus.UCSC.mm9.knownGene")

dir<-"~/Documents/Master_thesis/downstream_analysis"
setwd(dir)
blacklist.mm9<-("./bed_files/Blacklist.mm9.bed")

# Sample sheet
sampleSheet<-"sample_sheet.csv"
samples=read.csv(sampleSheet)
QCsamples = ChIPQC(samples, annotation="mm9")

save(QCsamples, file="ChIPQC_G4.RData")
RData<-file.path(getwd(), "ChIPQC_G4.RData")
load(RData)
# Show QCmetrics report 
QCtable<-QCmetrics(QCsamples)
write.table(QCtable, file="QCmetrics.xls")

FragmentLengthCrossCoverage(QCsamples)
RelativeCrossCoverage(QCsamples)
# Coverage histogram 
pdf("./figures/CoverageHistogram.pdf")
plotCoverageHist(QCsamples)
dev.off()

# Cross-coverage plot
pdf("./figures/CrossCoverage.pdf")
ccplot<-plotCC(QCsamples, colourBy = "Tissue", facetBy = "Factor")
ccplot+facet_wrap("Factor", scales="free_y")
dev.off()

# Relative enrichment of reads in annotated genomic intervals   
pdf("./figures/EnrichmentGenomicIntervals.pdf")
plotRegi(QCsamples)
dev.off()

# Plot peak profiles
pdf("./figures/PeakProfile.pdf")
plotPeakProfile(QCsamples)
dev.off() 

# Plot relative number of reads that overlap peaks vs. background reads
pdf("./figures/FractionReadsinPeaks.pdf")
plotRap(QCsamples)
dev.off()

# Plot fraction of reads in blacklisted regions
pdf("./figures/FractionReadsinBlacklistedRegions.pdf")
plotFribl(QCsamples)
dev.off()

# Sample clustering
pdf("./figures/CorrelationHeatmap.pdf")
plotCorHeatmap(QCsamples)
dev.off() 

# PCA
pdf("./figures/PCA.pdf")
plotPrincomp(QCsamples)
dev.off() 

