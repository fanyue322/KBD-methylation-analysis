#####violin plot for age correlation##############
library(ggplot2)
bmp("figure.tif",width=600,height=600)
ggplot(dat1, aes(x=age, y=methylation,fill=trait)) + 
  geom_violin() +
  geom_boxplot(width=0.05, fill="white")+
  labs(x="", y = "Correlations with age")+
   scale_color_manual(values=c(myCol[1],myCol[2]))+
   scale_fill_manual(name = "",values=c(myCol[1],myCol[2]))+
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed")+
  theme_classic(base_size = 15)+
  theme( axis.text = element_text(size = rel(2),color = "black"),
         axis.title = element_text(size = rel(2)),
         legend.text=element_text(size = rel(2)))
dev.off()

library(ChAMP)
beta=assays(rgSet)[["Beta"]]
RSobject <- RatioSet(beta, annotation = c(array = "IlluminaHumanMethylation450k", 
                                          annotation = "ilmn12.hg19"))
RSanno <- getAnnotation(RSobject)[, c("chr", "pos", "Name", 
                                      "UCSC_RefGene_Name","Enhancer")]
data(illumina450Gr)

############extract CpG sites from DMR#############
cpg.idx <- unique(unlist(apply(myDMR[[1]], 1, function(x) rownames(RSanno)[which(RSanno$chr == 
                                                                                         x[1] & RSanno$pos >= as.numeric(x[2]) & RSanno$pos <= 
                                                                                         as.numeric(x[3]))])))
############CpG sites distribution plot################
myResultsGr <- illumina450Gr[match(cpg.idx_KBD, names(illumina450Gr))]

idx=which(myResultsGr$feature=='5\'UTR')
idx=union(idx,which(myResultsGr$feature=='TSS200'))
idx=union(idx,which(myResultsGr$feature=='TSS1500'))
idx=union(idx,which(myResultsGr$feature=='1stExon'))
levels(myResultsGr$feature)=c(levels(illumina450Gr$feature),'Promoter')
myResultsGr[idx]$feature='Promoter'
cgiNames <- levels(illumina450Gr$cgi)
featureNames <- levels(illumina450Gr$feature)
featureNames=featureNames[c(2,4,5)]
featureNames=c('Promoter',featureNames)
sfo <- c(1,3,2,4)
scgio <- c(1, 4, 3, 2)
sfcgio <- rep((sfo - 1) * 4, each = 4) + rep(scgio,  4)
cout=1
for(i in 1:length(featureNames))
{
  for(j in 1:length(cgiNames))
  {
    lassoSizes[cout]=length(intersect(which(myResultsGr$cgi==cgiNames[j]),which(myResultsGr$feature==featureNames[i])))
    cout=cout+1
  }
}

bmp("figure.tif",width=900,height=600)
par(mar = c(7, 4, 4, 3) + 0.5)
plot(c(1, 16), y = c(range(0.3 * sqrt(lassoSizes))[1] * 
                       0.8, range(0.3 * sqrt(lassoSizes))[2] * 1.2), 
     type = "n", xaxt = "n", xlab = "", yaxt = "n", 
     ylab = "Number of CpG sites", main = "Kashin-Beck Disease",cex.main=rel(2),cex.lab=rel(2),font.main=1)
segments(1:16, rep(0, 16), 1:16, 0.3 * sqrt(lassoSizes[sfcgio]), 
         lty = 3, col = "grey")
points(1:16, 0.3 * sqrt(lassoSizes[sfcgio]), pch = 16, 
       cex = 0.3 * sqrt(lassoSizes[sfcgio]), col = rep(rainbow(4, 
                                                               alpha = 0.5)[sfo], each = 4))
text(1:16, 0.3 * sqrt(lassoSizes[sfcgio]), lassoSizes[sfcgio], 
     pos = 3, cex = 1.5)  
axis(1, at = 1:16, labels = rep(cgiNames[scgio], 
                                4), las = 2,cex.axis=1.5)
par(xpd = T)
# segments(seq(1, 28, 4), rep(-2.5, 7), seq(4, 28, 
#                                            4), rep(-2.5, 7))
mtext(text = featureNames[sfo], side = 1, at = seq(2.5, 
                                                   16, 4), line = 5.5, las = 1, cex=1.5)
axis(2, at = c(0, max(0.3 * sqrt(lassoSizes))), labels = F)
dev.off()
############veen plot###############
library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")
myCol=myCol[1:2]
venn.diagram(
  x = list(cpg.idx_KBD,cpg.idx_OA),
  category.names = c("KBD" , "OA"),
  filename = 'figure.png',
  #circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  #setnames
  main='Differential Methylated sites',
  main.cex=.6,
  main.fontface=1,
  cex = .6,
  fontface = 1,
  fontfamily = "sans",
  cat.cex = 0.6,
  cat.fontface = 1,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans",
  #    rotation = 1,
  output=TRUE,
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  resolution = 300,
  compression = "lzw"
)

##########DMR plot###############
library(DMRcate)
tmp=GRanges(seqnames='chr6',ranges=IRanges(start=168925150,end=169002604),strand='+')
lev2=c(2,4,6,8,11)###Control index
gcase='KBD'
caseind = which(grouplev %in% gcase)
lev1=caseind
CpGs1=rowSums(assays(rgSet)[["Beta"]][,lev1])-rowSums(assays(rgSet)[["Beta"]][,lev2])

gcase='OA'
caseind = which(grouplev %in% gcase)
lev1=caseind
CpGs2=rowSums(assays(rgSet)[["Beta"]][,lev1])-rowSums(assays(rgSet)[["Beta"]][,lev2])

phenotype=c('KBD','OA')
CpGs=cbind(CpGs1,CpGs2)
colnames(CpGs)=c('KBD','OA')
CpGs=abs(CpGs)

library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")
myCol=myCol[1:2]
cols=myCol

pdf("rplot.pdf",width=6,height=10)
DMR.plot(ranges=tmp,dmr=1,CpGs=CpGs,what="Beta",arraytype="450K",
         phen.col=cols,genome="hg19")
dev.off()


########GO plot##############
library(DOSE)
library(dplyr)
library(reshape2)
library(ggplot2)
KBD=data.frame(Concept.ID=go_KBD$ID,Concept.name=go_KBD$Description,ontology=go_KBD$ONTOLOGY,
               n.genes=go_KBD$Count,odds.ratio=go_KBD$GeneRatio,p.value=go_KBD$p.adjust,FDR=go_KBD$qvalue)
OA=data.frame(Concept.ID=go_OA$ID,Concept.name=go_OA$Description,ontology=go_OA$ONTOLOGY,
               n.genes=go_OA$Count,odds.ratio=go_OA$GeneRatio,p.value=go_OA$p.adjust,FDR=go_OA$qvalue)
dots <- function(x){
  if(x <= 0.05 & x > 0.0001){"0.0001 - 0.05"}
  else if(x <= 0.0001 & x > 1e-9){"1e-9 - 0.0001"}
  else if(x <= 1e-9){"< 1e-9"}
  else{"> 0.05"}
}

KBD$dots <- sapply(KBD$FDR, dots)
OA$dots <- sapply(OA$FDR, dots)
KBD$dots <- factor(KBD$dots, levels=c("> 0.05", "0.0001 - 0.05", "1e-9 - 0.0001",  "< 1e-9"))
OA$dots <- factor(OA$dots, levels=c("> 0.05", "0.0001 - 0.05", "1e-9 - 0.0001",  "< 1e-9"))
KBD$trait='KBD'
OA$trait='OA'

final=rbind(KBD,OA)
final$odds.ratio=as.numeric(as.character(final$odds.ratio))
names(final)[5]='gene.ratio'
names(final)[7]='qvalue'
final$trait=factor(final$trait)
bmp('test.tif',width=1200,height=600)
ggplot(final,aes(y=Concept.name, x=trait))+
  geom_tile(aes(fill=gene.ratio), colour='white')+
  scale_fill_gradient2(midpoint=0.125, low="blue", high="red")+
  theme(axis.text.x = element_text(size=rel(2),angle = -40, hjust=0, vjust=1),axis.title.x=element_blank(), 
        axis.text.y = element_text(size=rel(2)), axis.title.y=element_blank(),
        legend.title=element_text(size = rel(2), face = "bold", hjust = 0))+
  geom_point(data=final[final$dots != "> 0.05",], aes(size=dots), shape=21)+
  scale_size_discrete(range=c(0,3), name="qvalue")+
  coord_flip()+
scale_x_discrete(limits=rev(levels(final$trait)))
dev.off()
