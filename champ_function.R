champ.DMP2<-function (beta = myNorm, pheno = myLoad$pd$Sample_Group, compare.group = NULL, 
                      adjPVal = 0.05, adjust.method = "BH", arraytype = "450K",covariates=NULL) 
{
  message("[===========================]")
  message("[<<<<< ChAMP.DMP START >>>>>]")
  message("-----------------------------")
  CalculateDMP <- function(beta, pheno, tmp_compare, adjPVal = adjPVal, 
                           adjust.method = adjust.method) {
    message("  -----------------------------")
    message("  Start to Compare : ", tmp_compare[1], ", ", 
            tmp_compare[2])
    p <- pheno[which(pheno %in% tmp_compare)]
    tmpbeta <- beta[, which(pheno %in% tmp_compare)]
    design <- model.matrix(~0 + p)
    design <- cbind(design,covariates)
    contrast.matrix <- makeContrasts(contrasts = paste(colnames(design)[2:1], 
                                                       collapse = "-"), levels = colnames(design))
    message("  Contrast Matrix")
    print(contrast.matrix)
    fit <- lmFit(tmpbeta, design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    tryCatch(fit3 <- eBayes(fit2), warning = function(w) {
      stop("limma failed, No sample variance.\n")
    })
    DMP <- topTable(fit3, coef = 1, number = nrow(tmpbeta), 
                    adjust.method = adjust.method, p.value = adjPVal)
    message("  You have found ", sum(DMP$adj.P.Val <= adjPVal), 
            " significant MVPs with a ", adjust.method, " adjusted P-value below ", 
            adjPVal, ".")
    message("  Calculate DMP for ", tmp_compare[1], " and ", 
            tmp_compare[2], " done.")
    return(DMP)
  }
  message("!!! Important !!! New Modification has been made on champ.DMP(): \n")
  message("    (1): In this version champ.DMP() if your pheno parameter contains more than two groups of phenotypes, champ.DMP() would do pairewise differential methylated analysis between each pair of them. But you can also specify compare.group to only do comparasion between any two of them.\n")
  message("    (2): champ.DMP() now support numeric as pheno, and will do linear regression on them. So covariates like age could be inputted in this function. You need to make sure your inputted \"pheno\" parameter is \"numeric\" type.\n")
  Compare <- NULL
  message("--------------------------------")
  if (is.null(pheno) | length(unique(pheno)) <= 1) {
    stop("pheno parameter is invalid. Please check the input, pheno MUST contain at least two phenotypes.")
  }
  else {
    if (length(pheno) != ncol(beta)) 
      stop("Your Pheno's length is not in accord with your beta value's ncol.")
    message("\n[ Section 1:  Check Input Pheno Start ]\n")
    if (class(pheno) == "numeric") {
      message("  You pheno is numeric type.")
      message("    pheno parameter contains :", length(pheno), 
              " values.")
      message("    pheno parameter ranges from ", min(pheno), 
              " to ", max(pheno))
    }
    else {
      message("  You pheno is ", class(pheno), " type.")
      message("    Your pheno information contains following groups. >>")
      sapply(unique(pheno), function(x) message("    <", 
                                                x, ">:", sum(pheno == x), " samples."))
      message("    [The power of statistics analysis on groups contain very few samples may not strong.]")
      if (length(unique(pheno)) == 2) {
        message("    pheno contains only 2 phenotypes")
        if (is.null(compare.group)) {
          message("    compare.group parameter is NULL, two pheno types will be added into Compare List.")
          Compare <- list(x1 = unique(pheno))
        }
        else if (sum(compare.group %in% unique(pheno)) == 
                 2) {
          message("    Your compare.group parameter is in accord with your pheno. Two pheno types has been added into Compare List.")
          Compare <- list(x1 = unique(pheno))
        }
        else {
          stop(" You have specified compare.group, but it's not in accord with your pheno parameter. Please recheck your compare.group or pheno.")
        }
      }
      else if (length(unique(pheno)) > 2) {
        message("    pheno contains ", length(unique(pheno)), 
                " phenotypes")
        if (is.null(compare.group)) {
          message("    compare.group parameter is NULL, EACH PAIR of phenotypes will be added into Compare List.")
          Compare <- as.list(data.frame(combn(unique(pheno), 
                                              2)))
        }
        else if (sum(compare.group %in% unique(pheno)) == 
                 2) {
          message("    Your compare.group parameter is in accord with your pheno. Two pheno types has been added into Compare List.")
          Compare <- list(x1 = sort(compare.group))
        }
        else {
          stop("    Your pheno parameter contains multiple phenotypes, but values in your compare.group parameter are not all find in them. Please recheck your compare.group or pheno.")
        }
      }
      else {
        stop("    !!! Something wrong with your pheno. Please check your input.")
      }
      tmpnamelist <- vector()
      for (i in 1:length(Compare)) {
        tmpname <- paste(Compare[[i]][1], Compare[[i]][2], 
                         sep = "_to_")
        message("    ", tmpname, " compare group : ", 
                Compare[[i]][1], ", ", Compare[[i]][2])
        tmpnamelist <- c(tmpnamelist, tmpname)
      }
      names(Compare) <- tmpnamelist
    }
    message("\n[ Section 1:  Check Input Pheno Done ]\n")
  }
  DMPs <- list()
  if (is.null(Compare)) {
    message("\n[ Section 2:  Find Numeric Covariates Linear Regression CpGs Start ]\n")
    df <- data.frame(pheno = pheno)
    model.matrix <- model.matrix(~pheno, data = df)
    model.matrix <- cbind(model.matrix,covariates)
    fit1 <- lmFit(beta, model.matrix)
    fit2 <- eBayes(fit1)
    DMP <- topTable(fit2, coef = 2, number = nrow(beta), 
                    adjust.method = adjust.method, p.value = adjPVal)
    message("  You have found ", sum(DMP$adj.P.Val <= adjPVal), 
            " significant MVPs with a ", adjust.method, " adjusted P-value below ", 
            adjPVal, ".")
    if (sum(DMP$adj.P.Val <= adjPVal) != 0) 
      DMPs[["NumericVariable"]] <- DMP
  }
  else {
    message("\n[ Section 2:  Find Differential Methylated CpGs Start ]\n")
    for (i in names(Compare)) {
      DMP <- CalculateDMP(beta, pheno, Compare[[i]], adjPVal, 
                          adjust.method)
      if (sum(DMP$adj.P.Val <= adjPVal) != 0) 
        DMPs[[i]] <- DMP
    }
  }
  message("\n[ Section 2:  Find Numeric Vector Related CpGs Done ]\n")
  if (length(DMPs) == 0) 
    stop("ChAMP.DMP Have not detected even one significant CpGs. You may try other threshold.")
  message("\n[ Section 3:  Match Annotation Start ]\n")
  if (arraytype == "EPIC") 
    data(probe.features.epic)
  else data(probe.features)
  for (i in names(DMPs)) {
    com.idx <- intersect(rownames(DMPs[[i]]), rownames(probe.features))
    if (!is.null(Compare)) {
      avg <- cbind(rowMeans(beta[com.idx, which(pheno == 
                                                  Compare[[i]][1])]), rowMeans(beta[com.idx, which(pheno == 
                                                                                                     Compare[[i]][2])]))
      avg <- cbind(avg, avg[, 2] - avg[, 1])
      colnames(avg) <- c(paste(Compare[[i]], "AVG", sep = "_"), 
                         "deltaBeta")
      DMPs[[i]] <- data.frame(DMPs[[i]][com.idx, ], avg, 
                              probe.features[com.idx, ])
    }
    else {
      DMPs[[i]] <- data.frame(DMPs[[i]][com.idx, ], probe.features[com.idx, 
                                                                   ])
    }
  }
  message("\n[ Section 3:  Match Annotation Done ]\n")
  message("[<<<<<< ChAMP.DMP END >>>>>>]")
  message("[===========================]")
  message("[You may want to process DMP.GUI() or champ.GSEA() next.]\n")
  return(DMPs)
}

champ.DMR2<-function (beta = myNorm, pheno = myLoad$pd$Sample_Group, compare.group = NULL, 
                      arraytype = "450K", method = "Bumphunter", minProbes = 7, 
                      adjPvalDmr = 0.05, cores = 3, maxGap = 300, cutoff = NULL, 
                      pickCutoff = TRUE, smooth = TRUE, smoothFunction = loessByCluster, 
                      useWeights = FALSE, permutations = NULL, B = 250, nullMethod = "bootstrap", 
                      meanLassoRadius = 375, minDmrSep = 1000, minDmrSize = 50, 
                      adjPvalProbe = 0.05, Rplot = T, PDFplot = T, resultsDir = "./CHAMP_ProbeLasso/", 
                      rmSNPCH = T, fdr = 0.05, dist = 2, mafcut = 0.05, lambda = 1000, 
                      C = 2,covariates=NULL,beta2=myNorm2) 
{
  message("[===========================]")
  message("[<<<<< ChAMP.DMR START >>>>>]")
  message("-----------------------------")
  message("!!! important !!! We just upgrate champ.DMR() function, since now champ.DMP() could works on multiple phenotypes, but ProbeLasso can only works on one DMP result, so if your pheno parameter contains more than 2 phenotypes, and you want to use ProbeLasso function, you MUST specify compare.group=c(\"A\",\"B\"). Bumphunter and DMRcate should not be influenced.")
  message("\n[ Section 1:  Check Input Pheno Start ]\n")
  if (length(which(is.na(beta))) > 0) 
    message(length(which(is.na(beta))), " NA are detected in your beta Data Set, which may cause fail or uncorrect of SVD analysis. You may want to impute NA with champ.impute() function first.")
  if (!class(pheno) %in% c("character", "factor", "numeric")) 
    stop("pheno parameter must be a category vector, could be character, factor or numeric (numeric is not recommended).")
  message("  You pheno is ", class(pheno), " type.")
  message("    Your pheno information contains following groups. >>")
  sapply(unique(pheno), function(x) message("    <", x, ">:", 
                                            sum(pheno == x), " samples."))
  if (method == "ProbeLasso") {
    message("  ProbeLasso Method can only be done between two phenotypes. So we need to do more check here...")
    if (length(unique(pheno)) > 2) {
      message("    Your pheno contains more than two phenotypes.")
      message("    You may specify compare.group to do comparision between certain two phenotypes")
      if (is.null(compare.group)) {
        stop("    You did not specifically compare.group parameter to specify which two phenotypes you want to analysis.")
      }
      else if (sum(compare.group %in% unique(pheno)) == 
               2) {
        message("    Your compare.group is in accord with your pheno parameter, which is good.")
        message("    Now champ.DMR() would extract values for only these two phenotypes to analysis.")
        beta <- beta[, which(pheno %in% compare.group)]
        beta2<-beta2[, which(pheno %in% compare.group)]
        pheno <- pheno[which(pheno %in% compare.group)]
      }
      else {
        stop("    Seems you specified compare.group, but elements in your compare.group are not all found in your pheno parameter. Please recheck your pheno or compare.group.")
      }
    }
    else if (length(unique(pheno)) == 2) {
      message("    Your pheno parameter contains extactly two phenotypes, which is good and compare.group is not needed, champ.DMR() would proceed with your whole data set.")
    }
    else {
      stop("    Seems something wrong with your pheno data. champ.DMR() can not proceed. Please recheck your pheno information.")
    }
  }
  message("\n[ Section 1:  Check Input Pheno Done ]\n")
  message("\n[ Section 2:  Run DMR Algorithm Start ]\n")
  if (arraytype == "EPIC") {
    RSobject <- RatioSet(beta, annotation = c(array = "IlluminaHumanMethylationEPIC", 
                                              annotation = "ilm10b4.hg19"))
  }
  else {
    RSobject <- RatioSet(beta, annotation = c(array = "IlluminaHumanMethylation450k", 
                                              annotation = "ilmn12.hg19"))
  }
  probe.features <- getAnnotation(RSobject)
  if (cores > detectCores()) 
    cores <- detectCores()
  if (method == "Bumphunter") {
    message("<< Find DMR with Bumphunter Method >>")
    message(cores, " cores will be used to do parallel Bumphunter computing.")
    registerDoParallel(cores = cores)
    cpg.idx <- intersect(rownames(beta), rownames(probe.features))
    Anno <- probe.features[cpg.idx, ]
    Anno <- Anno[order(Anno$chr, Anno$pos), ]
    cpg.idx <- rownames(Anno)
    cl <- clusterMaker(Anno$chr, Anno$pos, maxGap = maxGap)
    names(cl) <- cpg.idx
    bumphunter.idx <- cpg.idx[which(cl %in% names(which(table(cl) > 
                                                          minProbes)))]
    message("According to your data set, champ.DMR() detected ", 
            sum(table(cl) > minProbes), " clusters contains MORE THAN ", 
            minProbes, " probes within", maxGap, " maxGap. These clusters will be used to find DMR.\n")
    X <- cbind(1, (as.numeric(as.factor(pheno)) - 1))
    Beta <- beta[bumphunter.idx, ]
    Beta <- replace(Beta, which(Beta <= 0.001), 0.001)
    Beta <- replace(Beta, which(Beta >= 0.999), 0.999)
    Y <- log((Beta/(1 - Beta)), 2)
    Bumps <- bumphunter(Y, design = X, chr = Anno[bumphunter.idx, 
                                                  ]$chr, pos = Anno[bumphunter.idx, ]$pos, cluster = cl[bumphunter.idx], 
                        cutoff = cutoff, pickCutoff = pickCutoff, smooth = smooth, 
                        smoothFunction = smoothFunction, useWeights = useWeights, 
                        permutations = permutations, verbose = TRUE, B = B, 
                        nullMethod = nullMethod)
    message("<< Calculate DMR success. >>")
    DMR <- Bumps$table[which(Bumps$table$p.valueArea <= adjPvalDmr), 
                       ]
    message("Bumphunter detected ", nrow(DMR), " DMRs with P value <= ", 
            adjPvalDmr, ".")
    if (nrow(DMR) == 0) 
      stop("No DMR detected.")
    rownames(DMR) <- paste("DMR", 1:nrow(DMR), sep = "_")
    DMR <- data.frame(DMR[, 1:3], width = DMR[, 3] - DMR[, 
                                                         2], strand = "*", DMR[, 4:14])
    colnames(DMR)[1:3] <- c("seqnames", "start", "end")
    OutputDMR <- list(BumphunterDMR = DMR)
  }
  else if (method == "ProbeLasso") {
    if (!file.exists(resultsDir)) 
      dir.create(resultsDir)
    message("champ.DMR Results will be saved in ", resultsDir)
    message("<< Find DMR with ProbeLasso Method >>")
    gc()
    DMP <- champ.DMP2(beta = beta, pheno = pheno, adjPVal = 1, covariates=covariates)
    if (length(DMP) > 1) 
      stop("Your pheno parameter seems contains more than 2 phenotypes. champ.DMR() only take covariates with only 2 phenotypes. Please manually extract your sample and covariates, then retry champ.DMR()")
    DMP <- DMP[[1]]
    if (arraytype == "EPIC") 
      data(illuminaEPICGr)
    else data(illumina450Gr)
    if (length(which(DMP$adj.P.Val < adjPvalProbe)) == 0) 
      stop("There is no probe show significant difference from champ.DMP() function.")
    myResultsGr <- illumina450Gr[match(rownames(DMP), names(illumina450Gr))]
    myResultsGr$P.Value <- DMP$P.Value[match(names(myResultsGr), 
                                             rownames(DMP))]
    myResultsGr$adj.P.Val <- DMP$adj.P.Val[match(names(myResultsGr), 
                                                 rownames(DMP))]
    seqlevels(myResultsGr) <- sort(seqlevels(myResultsGr))
    myResultsGr <- sort(myResultsGr, ignore.strand = T)
    myResultsGr$adj.P.Val <- p.adjust(mcols(myResultsGr)$P.Value, 
                                      method = "BH")
    message("<< Get closestProbe for each Probe >>")
    closestProbe <- as.data.frame(distanceToNearest(myResultsGr, 
                                                    ignore.strand = T))$distance
    closestProbeSp <- split(closestProbe, mcols(myResultsGr)$featureCgi)
    rm(closestProbe)
    message("<< Get lassoQuantileThreshold for each featureCgi >>")
    lassoQuantileDeviation <- abs(meanLassoRadius - rowMeans(as.data.frame(lapply(closestProbeSp, 
                                                                                  function(x) quantile(x, (1:1000)/1000)))))
    lassoQuantileThreshold <- which.min(lassoQuantileDeviation)/1000
    lassoSizes <- lapply(closestProbeSp, function(x) quantile(x, 
                                                              lassoQuantileThreshold, na.rm = T))
    message("<< Get expend ranges for each probe >>")
    myResultsGrSp <- split(myResultsGr, myResultsGr$featureCgi)
    lassoGr <- mapply(function(x, y) promoters(x, upstream = y, 
                                               downstream = y), x = myResultsGrSp, y = lassoSizes)
    lassoGr <- unlist(GRangesList(lassoGr))
    rm(myResultsGrSp)
    myResultsSigGr <- myResultsGr[which(mcols(myResultsGr)$adj.P.Val < 
                                          adjPvalProbe)]
    lassoProbeCountOverlap <- countOverlaps(lassoGr, myResultsSigGr, 
                                            ignore.strand = T)
    rm(myResultsSigGr)
    message("<< Get DMR from overlapped probes >>")
    dmrGr <- reduce(lassoGr[which(lassoProbeCountOverlap >= 
                                    minProbes)], min.gapwidth = minDmrSep, ignore.strand = TRUE)
    rm(lassoProbeCountOverlap, lassoGr)
    strand(dmrGr) <- "*"
    dmrGr <- dmrGr[which(width(dmrGr) > minDmrSize)]
    probeIndex <- as.data.frame(findOverlaps(dmrGr, myResultsGr))
    pValuesGr <- myResultsGr[probeIndex$subjectHits, "P.Value"]
    myBetas <- beta[match(names(pValuesGr), rownames(beta)), ]
    myBetas2<-beta2[match(names(pValuesGr), rownames(beta)), ]
    myBetas <- split(as.data.frame(myBetas), probeIndex$queryHits)
    myBetas2 <- split(as.data.frame(myBetas2), probeIndex$queryHits)
    message("<< Get adjusted P value for DMR >>")
    correl <- lapply(myBetas, function(x) cor(t(x)))
    weights <- lapply(correl, function(x) 1/apply(x^2, 1, 
                                                  sum))
    rm(correl)
    dmrQP <- qnorm(mcols(pValuesGr)$P.Value)
    dmrQP <- split(dmrQP, probeIndex$queryHits)
    dmrQPW <- mapply("*", dmrQP, weights)
    rm(dmrQP)
    if (class(dmrQPW) == "matrix") 
      dmrStat <- sum(dmrQPW)
    else dmrStat <- lapply(dmrQPW, sum)
    rm(dmrQPW)
    dmrSd <- lapply(weights, function(x) sqrt(sum(x^2)))
    rm(weights)
    dmrP <- mapply(function(x, y) pnorm(x, 0, sd = y), dmrStat, 
                   dmrSd)
    rm(dmrStat, dmrSd)
    dmrP <- p.adjust(dmrP, method = "BH")
    goodDmr <- which(dmrP < adjPvalDmr)
    dmrGr <- dmrGr[goodDmr]
    dmrP <- dmrP[goodDmr]
    dmrpRank <- rank(dmrP, ties.method = "min")
    rm(goodDmr)
    message("<< Get Start-End Ranges for each DMR >>")
    probeIndex <- as.data.frame(findOverlaps(dmrGr, myResultsGr))
    dmrProbesGr <- myResultsGr[probeIndex$subjectHits]
    myBetas <- beta[match(names(dmrProbesGr), rownames(beta)), ]
    myBetas2 <- beta2[match(names(dmrProbesGr), rownames(beta)), ]
    myBetas <- as.data.frame(myBetas)
    myBetas2 <- as.data.frame(myBetas2)
    dmrCoreStart <- start(dmrProbesGr)
    dmrCoreEnd <- end(dmrProbesGr)
    myBetas <- split(myBetas, probeIndex$queryHits)
    myBetas2 <- split(myBetas2, probeIndex$queryHits)
    dmrCoreStart <- split(dmrCoreStart, probeIndex$queryHits)
    dmrCoreStart <- sapply(dmrCoreStart, min)
    dmrCoreEnd <- split(dmrCoreEnd, probeIndex$queryHits)
    dmrCoreEnd <- sapply(dmrCoreEnd, max)
    message("<< Calculate Methylation Scores for each DMR >>")
    groupIndex <- pheno
    dmrGroupMeans <- do.call(rbind, lapply(myBetas, function(x) sapply(split(t(x), 
                                                                             groupIndex), mean)))
    dmrGroupMeans2 <- do.call(rbind, lapply(myBetas2, function(x) sapply(split(t(x), 
                                                                               groupIndex), mean)))
    colnames(dmrGroupMeans) <- paste("betaAv", colnames(dmrGroupMeans), 
                                     sep = "_")
    colnames(dmrGroupMeans2) <- paste("betaAv", colnames(dmrGroupMeans2), 
                                      sep = "_")
    probeGroupMeans <- lapply(myBetas, function(x) split(as.data.frame(t(x)), 
                                                         groupIndex))
    probeGroupMeans2<- lapply(myBetas2, function(x) split(as.data.frame(t(x)), 
                                                          groupIndex))
    rm(groupIndex, myBetas,myBetas2)
    probeGroupMeans <- lapply(probeGroupMeans, function(x) lapply(x, 
                                                                  colMeans))
    probeGroupMeans <- do.call(rbind, lapply(probeGroupMeans, 
                                             function(x) t(do.call(rbind, x))))
    probeGroupMeans2 <- lapply(probeGroupMeans2, function(x) lapply(x, 
                                                                    colMeans))
    probeGroupMeans2 <- do.call(rbind, lapply(probeGroupMeans2, 
                                              function(x) t(do.call(rbind, x))))
    colnames(probeGroupMeans) <- paste("betaAv", colnames(probeGroupMeans), 
                                       sep = "_")
    colnames(probeGroupMeans2) <- paste("betaAv", colnames(probeGroupMeans2), 
                                        sep = "_")
    message("<< Generate Probe-level Data >>")
    myDmrProbesGr <- myResultsGr[probeIndex$subjectHits]
    myDmrProbesGr <- as(cbind(as.data.frame(myDmrProbesGr), 
                              probeGroupMeans2), "GRanges")
    rm(probeGroupMeans)
    myDmrProbesGr$dmrNo <- probeIndex$queryHits
    myDmrProbesGr$dmrP <- dmrP[probeIndex$queryHits]
    myDmrProbesGr$dmrpRank <- dmrpRank[probeIndex$queryHits]
    myDmrProbesGr$dmrChrom <- seqnames(dmrGr[probeIndex$queryHits])
    myDmrProbesGr$dmrStart <- start(dmrGr[probeIndex$queryHits])
    myDmrProbesGr$dmrEnd <- end(dmrGr[probeIndex$queryHits])
    myDmrProbesGr$dmrSize <- width(dmrGr[probeIndex$queryHits])
    myDmrProbesGr$dmrCoreStart <- dmrCoreStart[probeIndex$queryHits]
    myDmrProbesGr$dmrCoreEnd <- dmrCoreEnd[probeIndex$queryHits]
    myDmrProbesGr$dmrCoreSize <- myDmrProbesGr$dmrCoreEnd - 
      myDmrProbesGr$dmrCoreStart + 1
    message("<< Generate DMR metadata >>")
    myDmrGr <- dmrGr
    myDmrGr$dmrNo <- unique(probeIndex$queryHits)
    myDmrGr$dmrP <- dmrP
    rm(dmrP)
    myDmrGr$dmrpRank <- dmrpRank
    rm(dmrpRank)
    myDmrGr$dmrChrom <- seqnames(dmrGr)
    myDmrGr$dmrStart <- start(dmrGr)
    myDmrGr$dmrEnd <- end(dmrGr)
    myDmrGr$dmrSize <- width(dmrGr)
    rm(dmrGr)
    myDmrGr$dmrCoreStart <- dmrCoreStart
    myDmrGr$dmrCoreEnd <- dmrCoreEnd
    myDmrGr$dmrCoreSize <- myDmrGr$dmrCoreEnd - myDmrGr$dmrCoreStart + 
      1
    genes <- split(as.data.frame(myResultsGr)[probeIndex$subjectHits, 
                                              c("ensemblID", "geneSymbol")], probeIndex$queryHits)
    rm(probeIndex)
    myDmrGr$ensemblID <- sapply(genes, function(x) paste(unique(unlist(strsplit(x$ensemblID, 
                                                                                ";"))), collapse = ";"))
    myDmrGr$geneSymbol <- sapply(genes, function(x) paste(unique(unlist(strsplit(x$geneSymbol, 
                                                                                 ";"))), collapse = ";"))
    rm(genes)
    myDmrGr <- as(cbind(as.data.frame(myDmrGr), dmrGroupMeans2), 
                  "GRanges")
    rm(dmrGroupMeans)
    interplot <- function(illumina450Gr, lassoSizes) {
      cgiNames <- levels(illumina450Gr$cgi)
      featureNames <- levels(illumina450Gr$feature)
      sfo <- c(6, 7, 3, 1, 4, 2, 5)
      scgio <- c(1, 4, 3, 2)
      sfcgio <- rep((sfo - 1) * 4, each = 4) + rep(scgio, 
                                                   7)
      lassoSizes <- round(unlist(lassoSizes))
      par(mar = c(7, 4, 4, 3) + 0.5)
      plot(c(1, 28), y = c(range(0.3 * sqrt(lassoSizes))[1] * 
                             0.8, range(0.3 * sqrt(lassoSizes))[2] * 1.2), 
           type = "n", xaxt = "n", xlab = "", yaxt = "n", 
           ylab = "lasso radius [bp]", main = paste("lasso quantile = ", 
                                                    round(lassoQuantileThreshold, 2), "\nmean lasso radius = ", 
                                                    meanLassoRadius, "bp", sep = ""), bty = "n")
      segments(1:28, rep(0, 28), 1:28, 0.3 * sqrt(lassoSizes[sfcgio]), 
               lty = 3, col = "grey")
      points(1:28, 0.3 * sqrt(lassoSizes[sfcgio]), pch = 16, 
             cex = 0.3 * sqrt(lassoSizes[sfcgio]), col = rep(rainbow(7, 
                                                                     alpha = 0.5)[sfo], each = 4))
      text(1:28, 0.3 * sqrt(lassoSizes[sfcgio]), lassoSizes[sfcgio], 
           pos = 3, cex = 0.8)
      axis(1, at = 1:28, labels = rep(cgiNames[scgio], 
                                      7), las = 2)
      par(xpd = T)
      segments(seq(1, 28, 4), rep(-2.5, 7), seq(4, 28, 
                                                4), rep(-2.5, 7))
      mtext(text = featureNames[sfo], side = 1, at = seq(2.5, 
                                                         28, 4), line = 5.5, las = 1, cex.axis = 1)
      axis(2, at = c(0, max(0.3 * sqrt(lassoSizes))), labels = F)
    }
    if (Rplot) 
      interplot(illumina450Gr, lassoSizes)
    if (PDFplot) {
      pdf(paste(resultsDir, "myLassos.pdf", sep = "/"), 
          width = 9, height = 9)
      interplot(illumina450Gr, lassoSizes)
      dev.off()
    }
    DMRProbes <- as.data.frame(myDmrProbesGr)
    DMRProbes <- data.frame(probe.features[rownames(DMRProbes), 
                                           ], DMRProbes[, which(colnames(DMRProbes) == "P.Value"):which(colnames(DMRProbes) == 
                                                                                                          "dmrNo")])
    DMRProbes <- split(DMRProbes, DMRProbes$dmrNo)
    DMR <- as.data.frame(myDmrGr)
    message("ProbeLasso detected ", nrow(DMR), " DMRs with P value <= ", 
            adjPvalDmr, ".")
    if (nrow(DMR) == 0) 
      stop("No DMR detected.")
    rownames(DMR) <- paste("DMR", DMR$dmrNo, sep = "_")
    names(DMRProbes) <- rownames(DMR)
    if (arraytype == "EPIC") 
      DMR[, 1] <- paste("chr", DMR[, 1], sep = "")
    OutputDMR <- list(ProbeLassoDMR = DMR)
  }
  else if (method == "DMRcate") {
    message(cores, " cores will be used to do parallel DMRcate computing.")
    message("<< Find DMR with DMRcate Method >>")
    myMs <- logit2(beta)
    if (rmSNPCH) 
      myMs <- rmSNPandCH(myMs, dist = dist, mafcut = mafcut)
    design <- model.matrix(~pheno)
    if (arraytype == "450K") {
      myannotation <- cpg.annotate(datatype = "array", 
                                   fdr = fdr, myMs, design = design, coef = ncol(design), 
                                   analysis.type = "differential", annotation = c(array = "IlluminaHumanMethylation450k", 
                                                                                  annotation = "ilmn12.hg19"), what = "M")
    }
    else {
      myannotation <- cpg.annotate(datatype = "array", 
                                   fdr = fdr, myMs, design = design, coef = ncol(design), 
                                   analysis.type = "differential", annotation = c(array = "IlluminaHumanMethylationEPIC", 
                                                                                  annotation = "ilm10b4.hg19"), what = "M")
    }
    M <- do.call("cbind", lapply(myannotation, as.data.frame))
    colnames(M) <- names(myannotation)
    dmrcoutput <- dmrcate(myannotation, min.cpgs = minProbes, 
                          lambda = lambda, C = C, mc.cores = cores)
    data(dmrcatedata)
    DMR <- as.data.frame(extractRanges(dmrcoutput, genome = "hg19"))
    message("DMRcate detected ", nrow(DMR), " DMRs with mafcut as= ", 
            adjPvalDmr, ".")
    if (nrow(DMR) == 0) 
      stop("No DMR detected.")
    if (nrow(DMR) != 0) {
      rownames(DMR) <- paste("DMR", 1:nrow(DMR), sep = "_")
      OutputDMR <- list(DMRcateDMR = DMR)
    }
    else {
      OutputDMR <- NULL
    }
  }
  else {
    stop("Please assign correct DMR method: 'Bumphunter' or 'ProbeLasso'")
  }
  message("\n[ Section 2:  Run DMR Algorithm Done ]\n")
  message("[<<<<<< ChAMP.DMR END >>>>>>]")
  message("[===========================]")
  message("[You may want to process DMR.GUI() or champ.GSEA() next.]\n")
  return(OutputDMR)
}
