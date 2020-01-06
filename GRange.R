########################################################
#Construction of RangedSummarizedExperiment object
########################################################

#construction of whole genome GRanges Object
hj1<-1:22
hj2<-LETTERS[24:25]
as.character(hj1)
hj<-paste("chr",c(hj1,hj2),sep="")
hj
bin_number<-read.table("chr_bin.txt",header = FALSE,sep = "\t") # the first and last 10 bins of each chromosome are not included
print(bin_number)
bin_size=bin_number$V2

#gr=GRanges(seqnames = rep("chr1",bin_size[1]-20),
            #ranges = IRanges(10:(bin_size[1]-10)*100000,11:(bin_size[1]-9)*100000-1)
           #)

print(bin_size)
gr=GRanges()
for (i in 1:24){
  size=bin_size[i]-20
  print(size)
  print(i)
  rowRanges<-GRanges(seqnames = Rle(rep(hj[i],size)),
                     ranges = IRanges(10:(size+9)*100000,11:(size+10)*100000-1),  #be careful of the number in IRanges 
                     strand = Rle(rep("*",size)))
  print(rowRanges)          
  gr<-c(gr,rowRanges)
}
gr

######
#construction of whole genome column information
######

coldata <- DataFrame(Treatment=c("mono3", "mono6","mac3","mac4"))    
coldata

######
#construction of whole genome RangedSummarizedExperiment Object
######

up_count=read.table(file ="raw_up.txt",header = FALSE,sep = "\t")
down_count=read.table(file ="raw_down.txt",header = FALSE,sep = "\t")
up_matrix<-as.matrix(up_count)
down_matrix<-as.matrix(down_count)


finder<-SummarizedExperiment(assays=list(up=up_matrix,down=down_matrix),rowRanges=gr, colData=coldata)
colnames(up_matrix)<-NULL
colnames(down_matrix)<-NULL        #Remember when import matrix,set colnames to NULL before construction
type(assays(finder)$up)
is.matrix(assays(finder)$up)



########################################################
#Matrix conbination
########################################################

all_counts<-cbind(assay(finder,"up"),assay(finder,"down"))  #equal to assays(finder)$up 
dim(all_counts)
total1<-c(49682034,26941899,24072245,43353691)
total2<-c(49182858,26679140,23824160,42985576)     # total counts that don't include the first and last 1MB regions

colnames(all_counts)<-LETTERS[1:8] 
ydom<-DGEList(all_counts,lib.size=rep(total1, 2))
ydom





########################################################
#Constructing the design matrix 
########################################################

input<-c("mono", "mono", "mac", "mac")              # The order c("mono3", "mono6", "mac3", "mac4")
Condition <- factor(c("mono", "mono", "mac", "mac"))
Sample <- factor(seq_along(input))
Condition <- rep(Condition, 2)
Sample <- rep(Sample, 2)
Direction <- rep(c("Up", "Down"), each=length(input))
Direction
design <- model.matrix(~0 + Sample + Direction:Condition)
design
design <- design[,!grepl("DirectionDown", colnames(design))]
design
colnames(design) <- sub("DirectionUp:", "", colnames(design))
design

########################################################
#Pre-processing of the count data
########################################################

ab <- aveLogCPM(ydom)
keep <- ab > 0
ydom <- ydom[keep,]
summary(keep)
cur.regions <- rowRanges(finder)[keep,]
summary(cur.regions)
cur.regions

########################################################
#Testing for significant differences with edgeR
########################################################
ydom <- estimateDisp(ydom, design)
plotBCV(ydom)
ydom

fitdom <- glmQLFit(ydom, design, robust=TRUE)
plotQLDisp(fitdom)

con <- makeContrasts(Conditionmono - Conditionmac, levels=design)
resdom <- glmQLFTest(fitdom, contrast=con)
topTags(resdom)

output <- data.frame(as.data.frame(cur.regions)[,1:3],
  Mono=fitdom$coefficients[,"Conditionmono"]/log(2),
  Mac=fitdom$coefficients[,"Conditionmac"]/log(2),
  resdom$table)
output$FDR <- p.adjust(resdom$table$PValue, method="BH")
head(output)

o <- order(output$PValue)
o
output <- output[o,]
head(output)                            #sort the output table with the order of PValue

head(output[,1:5], 10)
sig <- output$FDR <= 0.05
summary(sig)

change.type <- abs(output$Mac) > abs(output$Mono)                       # To show that interaction strength of Mono bigger than Mac
summary(change.type[sig])

tophit <- cur.regions[o[1]]                                             # Most dynamic region between Mono anc Mac
tophit
expanded <- resize(tophit, fix="center", width=width(tophit)*10)        # 1Mb region aroun Top1 dynamic region
expanded


hg.frag<-cutGenome(BSgenome.Hsapiens.UCSC.hg19,"GATC",4)                # Produce the fragment of MboI
hg.param<-pairParam(hg.frag)                                                














        