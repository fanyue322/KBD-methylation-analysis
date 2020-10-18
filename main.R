library(ChAMP)
library('data.table')
source('./champcode.R')
####data processing##################
rgSet <- read.metharray.exp("./idatfile")

detP <- detectionP(rgSet)
remove <- apply(detP, 1, function (x) any(x > 0.01))
rgSet <- rgSet[!remove,]
rgSet <- preprocessFunnorm(rgSet)
rgSetms <- getM(rgSet)
rgSetms.noSNPs <- rmSNPandCH(rgSetms, dist=2, mafcut=0.05)
pheno=phenotype ### assigning the phenotype
pheno=as.factor(pheno)
rgSet$replicate=pheno
rgSet<-rgSet[rownames(rgSetms.noSNPs),]

age<-age ### assigning age
sex<-age ### assigning sex
age=scale(age)
rgSet$age=age
rgSet$sex=sex

adjPVal=0.05
adjust.method='BH'
arraytype='450K'
pheno=rgSet$replicate
compare.group=c("KBD","N")  ##according to your phenotypes, KBD/OA/Control
Compare <- list(x1 = sort(compare.group))
idx=which(pheno %in% Compare[[1]])
age<-rgSet$age[which(pheno %in% Compare[[1]])]
sex<-rgSet$sex[which(pheno %in% Compare[[1]])]
age=scale(age)
age=as.numeric(age)
myNorm=getM(rgSet)  # Mvalue
myNorm2=assays(rgSet)[["Beta"]] # Beta value
myNorm=myNorm[,idx]
myNorm2=myNorm2[,idx]
pheno=pheno[which(pheno %in% Compare[[1]])]
pheno=as.character(pheno)
pheno=as.factor(pheno)
covariates=data.frame(sex=sex)
myDMP_KBD<-champ.DMP2(beta=myNorm,pheno=pheno,covariates=covariates)       ###DMP analysis
myDMR_KBD<-champ.DMR2(beta=myNorm,pheno=pheno,arraytype="450K",method="ProbeLasso",
                      minProbes=3,adjPvalDmr=0.05,cores=3,meanLassoRadius=375,minDmrSep=1000,
                      minDmrSize=50,adjPvalProbe=1,Rplot=T,PDFplot=T,resultsDir="./CHAMP_ProbeLasso",covariates=covariates,beta2=myNorm2)  ##DMR analysis
######enrichment analysis#################
myDMR_KBD=myDMR_KBD[[1]]
myDMR_KBD$diff=abs(tmp1$betaAv_KBD-tmp1$betaAv_N)
myDMR_KBD_PMD=tmp3[which(myDMR_KBD$diff>0.2),]  

sig=data.frame()
for(i in 1:nrow(myDMR_KBD_PMD))
{
  tmp=myDMR_KBD_PMD[i,]
  gene=tmp$ensemblID
  gene=strsplit(gene,';')
  res=data.frame(id=gene,Pvalue=tmp$dmrP)
  names(res)=c('id','Pvalue')
  sig=rbind(sig,res)
}         #########Mapping the DMR to the corresponding genes
########GO enrichment analysis############
id <- bitr(sig$ENSEMBL,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
idx=match(id$ENSEMBL,sig$id)
names(sig)=c('ENSEMBL','Pvalue')
data <- merge(id,sig,by="ENSEMBL")
data <- data[order(data$Pvalue,decreasing=T),]
data <- subset(data,select=c("ENTREZID","Pvalue"))
geneList = data[,1]
names(geneList) = as.character(data[,1])
geneList = sort(geneList, decreasing = TRUE)
go_padj <- 0.05
go_result <- enrichGO(OrgDb="org.Hs.eg.db", 
                      gene = names(geneList), #a vector of entrez gene id.
                      pvalueCutoff = go_padj,
                      keyType = "ENTREZID",
                      pAdjustMethod = "BH", #one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none
                      qvalueCutoff = 0.2,## qvalueCutoff
                      ont = "ALL",
                      minGSSize = 10, #minimal size of genes annotated by Ontology term for testing
                      maxGSSize = 500, #maximal size of genes annotated for testing
                      pool=T,
                      readable=TRUE)#whether mapping gene ID to gene Name
########KEGG enrichment analysis############  
kegg_padj <- 0.05
kegg_result <- enrichKEGG(gene = names(geneList), organism ="hsa",
                          pvalueCutoff = kegg_padj,
                          pAdjustMethod = "BH",
                          qvalueCutoff =0.2,
                          minGSSize = 3,
                          use_internal_data=F)
pv.out <- pathview(gene.data = names(geneList), pathway.id =kegg_result$ID[1], species ="hsa",Map.null=TRUE)

          
          
