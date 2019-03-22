################### 

library(Rmpi)
library(snow)
cl<-makeCluster(8, type="MPI") # Do not use

###################

library(GenomicAlignments)
library(ChIPQC)
library(GenomicRanges)
library("TxDb.Mmusculus.UCSC.mm9.knownGene")
library(BiocParallel)
register(SerialParam())


dir<-"/proj/uppstore2018208/KRISTINA/PROJECT_FILES"
setwd(dir)
blacklist.mm9<-read.table("./bed_files/Blacklist.mm9.bed")

# Sample sheet
sampleSheet<-"./sample_sheet_G4.csv" # Input sample sheet
samples=read.csv(file=sampleSheet, sep=",", na.strings="NA")

# Run ChIPQC function
QCsamples = clusterCall(cl, ChIPQC(samples, annotation="mm9", blacklist=blacklist.mm9, chromosomes=paste0("chr", c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,"X","Y"))))

# Show QCmetrics  
QCmetrics(QCsamples)
# Save QCmetrics to xls 
QCtable<-QCmetrics(QCsamples)
write.table(QCtable, file="QCmetrics.xls", sep="\t", col.names=NA)

FragmentLengthCrossCoverage(QCsamples)
RelativeCrossCoverage(QCsamples)

# Make QC report 
ChIPQCreport(QCsamples)

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

