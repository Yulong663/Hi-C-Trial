BiocManager::install("GenomicRanges")
BiocManager::install("IRanges")
BiocManager::install("InteractionSet")
BiocManager::install("diffHic")
BiocManager::install("rlang")
BiocManager::install("SummarizedExperiment")
BiocManager::install("statmod")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")

library("GenomicRanges")
library("IRanges")
library("InteractionSet")
library("diffHic")
library("edgeR")
library("SummarizedExperiment")
library("statmod")
library("BSgenome.Hsapiens.UCSC.hg19")

setwd("D:/R studio/diffHic")
# below are test regions

#example of IRanges object
IRanges(0:9*10+1, 1:10*10)

#example of GRanges object
all.regions <- GRanges("chrA", IRanges(0:9*10+1, 1:10*10))
all.regions
index.1 <- c(1,5,10)
index.2 <- c(3,2,6)
region.1 <- all.regions[index.1]
region.2 <- all.regions[index.2]
region.1
region.2

## more complicate example
gr <- GRanges(
  seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
  ranges = IRanges(101:110, end = 111:120, names = head(letters, 10)),
  strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
  score = 1:10,
  GC = seq(1, 0, length=10))

gr<-c(gr[1:6],gr)      # combine GRanges Object
gr<-c(GRanges(),gr)


#example of GInteractions object
gi.1 <- GInteractions(region.2, region.2)
gi.1
gi.2 <- GInteractions(all.regions, all.regions)
gi.2
gi<-c(gi.1,gi.2)
gi
regions(gi)

#example of ContactMatrix object
row.indices <- 1:10
col.indices <- 9:10
row.regions <- all.regions[row.indices]
col.regions <- all.regions[col.indices]
Nr <- length(row.indices)
Nc <- length(col.indices)
counts <- matrix(rpois(Nr*Nc, lambda=10), Nr, Nc)
cm <- ContactMatrix(counts, row.regions, col.regions)
cm
cm[1:5,]
as.matrix(rbind(cm,cm[1:5,]))

#example run of diffHic 
mono_chr1<-read.table("mono_chr1.txt",header = FALSE,sep="\t")
mono1<-mono_chr1[,4:6235]
data1<-as.matrix(mono1)


mac_chr1<-read.table("mac_chr1.txt",header=FALSE,sep="\t")
mac1<-mac_chr1[,4:6235]
data2<-as.matrix(mac1)

print(dim(data1))
print(dim(data2))
data2[1,1:3]
data1[1,1]
dim(data1)

chr1_regions <- GRanges("chr1", IRanges(0:6231*40000, 1:6232*40000-1))
chr1_regions
cm1<-ContactMatrix(data1,chr1_regions,chr1_regions)
cm2<-ContactMatrix(data2,chr1_regions,chr1_regions)
regions(cm1)
mergeCMs(cm,cm,filter=5)

# example of RangedSummarizedExperiment object

## construction of RangedSummarizedExperiment object
nrows <- 200
ncols <- 6
count <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
rowranges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                     IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                     strand=sample(c("+", "-"), 200, TRUE),
                     feature_id=sprintf("ID%03d", 1:200))
colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                     row.names=LETTERS[1:6])

sn<-SummarizedExperiment(assays=list(count=counts),rowRanges=rowRanges, colData=colData)


sn$counts
LETTERS[1:6]

# construction of Design Matrix
Condition <- factor(c("flox", "flox", "ko", "ko"))
Sample <- factor(seq_along(input))
Condition <- rep(Condition, 2)
Sample <- rep(Sample, 2)
Direction <- rep(c("Up", "Down"), each=length(input))
design <- model.matrix(~0 + Sample + Direction:Condition)



