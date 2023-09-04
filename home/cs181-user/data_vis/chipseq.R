library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)  # Genome annotation package for D. melanogaster

# Load the ChIP-seq peaks from a BED file
peaks <- readPeakFile('chip_peaks.bed')

# Create a GRanges object from the peaks
gr_peaks <- GRanges(seqnames = peaks$seqnames, 
                    ranges = IRanges(start = peaks$start, end = peaks$end))

# Create a genome-wide profile plot
plotProfile(gr_peaks, annotations = TxDb.Dmelanogaster.UCSC.dm6.ensGene,
            ylab = "Peak Density", main = "Genome-wide ChIP-seq Profile")

# Load the necessary libraries
library(ChIPpeakAnno)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)  # Genome annotation package for D. melanogaster

# Load the ChIP-seq peaks from a BED file
peaks <- readPeakFile('chip_peaks.bed')

# Annotate ChIP-seq peaks with nearby genes
annotated_peaks <- annotatePeakInBatch(peaks, AnnotationData = TxDb.Dmelanogaster.UCSC.dm6.ensGene)

# Plot gene annotation summary
plotAnnoBar(annotated_peaks)

# Plot distribution of distances to TSS
plotDistToTSS(annotated_peaks)