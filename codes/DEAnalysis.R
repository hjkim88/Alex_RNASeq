###
#   File name : DEAnalysis.R
#   Author    : Hyunjin Kim
#   Date      : Feb 7, 2020
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Perform DE analysis on the raw counts
#
#   * The comparisons would be:
#   1. Stenosis vs Sclerosis
#
#   Instruction
#               1. Source("DEAnalysis.R")
#               2. Run the function "dea_hans" - specify raw count path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_DEAnalysis.R/DEAnalysis.R")
#               > dea_hans(rCntPath="./data/raw_counts.rda",
#                          outputDir="./results/differential_expression/")
###

dea_hans <- function(rCntPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Ferrari/RNASeq_Analysis/data/raw_counts.rda",
                     outputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Ferrari/RNASeq_Analysis/results/differential_expression/") {
  
  ### load library
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  
  ### load raw count data
  load(rCntPath)
  
  
  ####################################################
  ### A function to perform DE analysis with limma ###
  ####################################################
  #' @title limmaWithComparisons
  #' @param rawData A numeric matrix where rows are genes and columns are samples and where
  #'		normCnt[i, j] is an expression (raw count or log-scale normalized value) of the
  #'		i-th gene on the j-th sample. colnames(rawData) provides the sample names and
  #'		rownames(rawData) contains the gene ids (either entrez ids or gene symbols).
  #' @param grp A character vector of length ncol(rawData) where names(grp) = colnames(rawData). 
  #' 		It assigns each sample in rawData to a class. Or a list of two character vectors
  #' 		representing two different independent classes.
  #' @param exp_class A chracter string, must be a member of the string vector unique(grp).
  #'		Specifies which annotation class to use as the "case" in the comparison.
  #'		Or a vector with two character strings that are from the grp[[1]] for a comparison.
  #' @param ctrl_class A chracter string, must be a member of the string vector unique(grp).
  #'		Specifies which annotation class to use as the "reference" in the comparison.
  #'		Or a vector with two character strings that are from the grp[[2]] for a comparison.
  #'		
  #'		* If exp_class = "A" and ctrl_class = "B", the direction of
  #'      differential expression is A - B. If exp_class = c("A", "B") and
  #'      and ctrl_class = c("C", "D"), then the direction of
  #'      differential expression is (AC - AD) - (BC - BD).
  #'      
  #' @param analysis_type A type of differential expression analysis.
  #'    Should be either "Normalized" or "RawCounts".
  #' @param bat_eff If NULL, no batch effect is used in the analysis. Otherwise, it is a 
  #'		a character vector of length ncol(rawData) specifying batch membership for each 
  #'		sample. Specifically:
  #'		- names(bat_eff) = colnames(rawData)
  #'		- The values bat_eff[i] represent batch ids 
  #'		- If bat_eff[i] == bat_eff[j], this means that the samples names(bat_eff)[i]
  #'			and names(bat_eff)[j] come from the same batch.
  #' @param thresh Numeric. Filters out from the results genes with adjusted
  #' 		p-value larger than this value
  #' @adj_method: A method for adjusting the p-values to account for multiple testing.
  #'   	Must be either "BH" (Benjamini-Hochberg) or "BY" (Bonferroni) or "holm" (HOLM).
  #' @return data.frame Data frame comprising the results of the differential gene 
  #' 		expression analysis. Each row represents a gene and contains comparison 
  #' 		statistics.
  #' @export
  #' @author Hyunjin Kim
  ####################################################
  limmaWithComparisons <- function(rawData, grp, exp_class, ctrl_class,
                                   analysis_type = c("RawCounts", "Normalized"), bat_eff = NULL, 
                                   thresh = 1, adj_method = c("BH", "BY", "holm")) {
    
    ### load library
    if(!require(edgeR)) {
      source("https://bioconductor.org/biocLite.R")
      biocLite("edgeR")
      library(edgeR)
    }
    if(!require("checkmate", quietly = TRUE)) {
      install.packages("checkmate")
      require("checkmate", quietly = TRUE)
    }
    
    ### argument checking
    assertDataFrame(rawData)
    assert(checkVector(grp), checkList(grp))
    if(testList(grp)) {
      assert(length(grp) == 2)
      assert(length(exp_class) == 2)
      assert(length(ctrl_class) == 2)
      assert(length(which(grp[[1]] == exp_class[1])) > 0)
      assert(length(which(grp[[1]] == exp_class[2])) > 0)
      assert(length(which(grp[[2]] == ctrl_class[1])) > 0)
      assert(length(which(grp[[2]] == ctrl_class[2])) > 0)
    } else {
      assert(length(exp_class) == 1)
      assert(length(ctrl_class) == 1)
      assert(length(which(grp == exp_class)) > 0, length(which(grp == ctrl_class)) > 0)
    }
    assertCharacter(exp_class)
    assertCharacter(ctrl_class)
    assertChoice(analysis_type, c("RawCounts", "Normalized"))
    assert(checkNull(bat_eff), checkCharacter(bat_eff))
    assertNumeric(thresh)
    assertChoice(adj_method, c("BH", "BY", "holm"))
    
    ### sometimes, there are some variables which can not be transformed into R variable names
    ### so, just to be safe, change all the variables to R-usuable ones
    exp_class <- make.names(exp_class)
    ctrl_class <- make.names(ctrl_class)
    
    ### there are two different grp options
    if(testList(grp)) {
      sampleType1 <- make.names(as.character(grp[[1]]))
      sampleType2 <- make.names(as.character(grp[[2]]))
      
      ### there are two options: considering batch effect or not
      if(is.null(bat_eff)) {
        ### make a data frame for design matrix
        Coldata <- data.frame(sampleType1, sampleType2)
        
        ### set names and the refererence group
        rownames(Coldata) <- colnames(rawData)
        Coldata$sampleType1 <- relevel(Coldata$sampleType1, ref = exp_class[2])
        Coldata$sampleType2 <- relevel(Coldata$sampleType2, ref = ctrl_class[2])
        
        ### prepare design matrix
        design <- model.matrix(~0+sampleType1*sampleType2, data = Coldata)
      } else {
        ### make a data frame for design matrix
        batch_eff <- as.character(bat_eff)
        Coldata <- data.frame(sampleType1, sampleType2, batch_eff)
        
        ### set names and the refererence group
        rownames(Coldata) <- colnames(rawData)
        Coldata$sampleType1 <- relevel(Coldata$sampleType1, ref = exp_class[2])
        Coldata$sampleType2 <- relevel(Coldata$sampleType2, ref = ctrl_class[2])
        
        ### prepare design matrix
        design <- model.matrix(~0+sampleType1*sampleType2+batch_eff, data = Coldata)
      }
      
      ### if it is a RNA-Seq data, we need to normalize it 
      if(analysis_type == "RawCounts") {
        ### raw counts to DGEList object
        d <- DGEList(rawData)
        
        ### remove 0 or low counts
        keep <- filterByExpr(d, design)
        d <- d[keep,,keep.lib.sizes=FALSE]
        
        ### calculate normalization factors
        d <- calcNormFactors(d)
        
        ### voom
        d <- voomWithQualityWeights(d, design)
      } else {
        d <- rawData
      }
      
      ### fit the linear model
      fit <- lmFit(d, design)
      fit2 <- eBayes(fit)
      
      ### get the differentially expressed genes
      result <- topTable(fit2, coef = paste0("sampleType1", exp_class[1], ":sampleType2", ctrl_class[1]),
                         adjust.method=adj_method[1], number=Inf)
    } else {
      sampleType <- relevel(as.factor(make.names(as.character(grp))), ref = ctrl_class)
      
      ### there are two options: considering batch effect or not
      if(is.null(bat_eff)) {
        ### make a data frame for design matrix
        Coldata <- data.frame(sampleType)
        
        ### set names and the refererence group
        rownames(Coldata) <- colnames(rawData)
        Coldata$sampleType <- relevel(Coldata$sampleType, ref = ctrl_class)
        
        ### prepare design matrix
        design <- model.matrix(~0+sampleType, data = Coldata)
        colnames(design) <- levels(sampleType)
      } else {
        ### make a data frame for design matrix
        batch_eff <- as.character(bat_eff)
        Coldata <- data.frame(sampleType, batch_eff)
        
        ### set names and the refererence group
        rownames(Coldata) <- colnames(rawData)
        Coldata$sampleType <- relevel(Coldata$sampleType, ref = ctrl_class)
        
        ### prepare design matrix
        design <- model.matrix(~0+sampleType+batch_eff, data = Coldata)
        colnames(design) <- c(levels(sampleType), levels(bat_eff)[-1])
      }
      
      ### if it is a RNA-Seq data, we need to normalize it 
      if(analysis_type == "RawCounts") {
        ### raw counts to DGEList object
        d <- DGEList(rawData)
        
        ### remove 0 or low counts
        keep <- filterByExpr(d, design)
        d <- d[keep,,keep.lib.sizes=FALSE]
        
        ### calculate normalization factors
        d <- calcNormFactors(d)
        
        ### voom
        d <- voomWithQualityWeights(d, design)
      } else {
        d <- rawData
      }
      
      ### fit the linear model
      fit <- lmFit(d, design)
      
      ### extract specific comparison of interest
      contrastMat <- makeContrasts(contrasts = paste(exp_class, ctrl_class, sep = "-"), levels = design)
      
      ### fit the contrasts
      fit2 <- contrasts.fit(fit, contrastMat)
      fit2 <- eBayes(fit2)
      
      ### get the differentially expressed genes
      result <- topTable(fit2, adjust.method=adj_method[1], number=Inf)
    }
    
    ### order based on adj.p.val and filter out low significance genes
    result <- result[order(result$adj.P.Val),]
    result <- result[result$adj.P.Val <= thresh, ,drop = FALSE]
    
    ### get normalized gene expression data
    if(analysis_type == "RawCounts") {
      ge <- d$E
    } else {
      ge <- d
    }
    
    ### there are two different grp options
    if(testList(grp)) {
      ### add baseMean for each group
      exp1_ctrl1_rowMeans <- apply(ge[,intersect(which(Coldata$sampleType1 == exp_class[1]),
                                                 which(Coldata$sampleType2 == ctrl_class[1])),
                                      drop=FALSE], 1, mean)
      exp1_ctrl2_rowMeans <- apply(ge[,intersect(which(Coldata$sampleType1 == exp_class[1]),
                                                 which(Coldata$sampleType2 == ctrl_class[2])),
                                      drop=FALSE], 1, mean)
      exp2_ctrl1_rowMeans <- apply(ge[,intersect(which(Coldata$sampleType1 == exp_class[2]),
                                                 which(Coldata$sampleType2 == ctrl_class[1])),
                                      drop=FALSE], 1, mean)
      exp2_ctrl2_rowMeans <- apply(ge[,intersect(which(Coldata$sampleType1 == exp_class[2]),
                                                 which(Coldata$sampleType2 == ctrl_class[2])),
                                      drop=FALSE], 1, mean)
      result <- data.frame(V1=exp1_ctrl1_rowMeans[rownames(result)],
                           V2=exp1_ctrl2_rowMeans[rownames(result)],
                           V3=exp2_ctrl1_rowMeans[rownames(result)],
                           V4=exp2_ctrl2_rowMeans[rownames(result)],
                           result[,1:6],
                           stringsAsFactors = FALSE, check.names = FALSE)
      colnames(result)[1:4] <- c(paste0("normMean_", exp_class[1], "_", ctrl_class[1]),
                                 paste0("normMean_", exp_class[1], "_", ctrl_class[2]),
                                 paste0("normMean_", exp_class[2], "_", ctrl_class[1]),
                                 paste0("normMean_", exp_class[2], "_", ctrl_class[2]))
    } else {
      ### add baseMean for each group
      exp_rowMeans <- apply(ge[,which(Coldata$sampleType == exp_class), drop=FALSE], 1, mean)
      ctrl_rowMeans <- apply(ge[,which(Coldata$sampleType == ctrl_class), drop=FALSE], 1, mean)
      result <- data.frame(V1=exp_rowMeans[rownames(result)],
                           V2=ctrl_rowMeans[rownames(result)],
                           result[,1:6],
                           stringsAsFactors = FALSE, check.names = FALSE)
      colnames(result)[1:2] <- c(paste0("normMean_", exp_class), paste0("normMean_", ctrl_class))  
    }
    
    return(result)
    
  }
  
  ### A function to print volcano plot of DE analysis with DE result
  volPlotWithDeRes <- function(deresult, outputFilePath, src=c("DESeq", "limma"), fdr=0.05, lfc=0.6) {
    
    ### load library
    if(!require(ggplot2)) {
      install.packages("ggplot2")
      library(ggplot2)
    }
    
    if(src == "DESeq") {
      deresult$padj[which(is.na(deresult$padj))] <- 1
      volcanoData <- as.data.frame(cbind(deresult$log2FoldChange, -log10(deresult$padj), as.character(deresult$padj < fdr & abs(deresult$log2FoldChange) > lfc)))
    } else if(src == "limma") {
      deresult$adj.P.Val[which(is.na(deresult$adj.P.Val))] <- 1
      volcanoData <- as.data.frame(cbind(deresult$logFC, -log10(deresult$adj.P.Val), as.character(deresult$adj.P.Val < fdr & abs(deresult$logFC) > lfc)))
    } else {
      stop("ERROR: the src parameter should be either \"DESeq\" or \"limma\".")
    }
    
    colnames(volcanoData) <- c("logFC", "logFDR", "Significance")
    volcanoData$logFC <- as.numeric(as.character(volcanoData$logFC))
    volcanoData$logFDR <- as.numeric(as.character(volcanoData$logFDR))
    
    ggplot(data=volcanoData, aes(x=logFC, y=logFDR, colour=Significance)) +
      xlab("log2_Fold_Change") +
      ylab("-log10(FDR)") +
      ggtitle(paste("Significant (adj.p < ", fdr, "& |logFC| > ", lfc, ") DE genes -", sum(volcanoData$Significance == "TRUE"))) +
      theme_classic(base_size = 16) +
      geom_point(alpha=0.4)
    ggsave(filename = outputFilePath, width = 12, height = 10)
  }
  
  
  ### 1. Stenosis vs Sclerosis
  ### limma-voom
  deresult <- limmaWithComparisons(rawData = rawCnt[,-c(1,2)], grp = sample_info,
                                   exp_class = "Stenosis",
                                   ctrl_class = "Sclerosis",
                                   analysis_type = "RawCounts", adj_method = "BH",
                                   bat_eff = NULL, thresh = 1)
  ### write out the DE result tables
  write.xlsx2(data.frame(Entrez_ID=rownames(deresult), Gene_Symbol=rawCnt[rownames(deresult),"Gene_Symbol"], deresult,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "DE_result_limma_voom_Stenosis_vs_Sclerosis.xlsx"),
              sheetName = "DE_Result", row.names = FALSE)
  ### draw a volcano plot
  volPlotWithDeRes(deresult, outputFilePath = paste0(outputDir, "Volcano_plot_Stenosis_vs_Sclerosis.png"),
                   src = "limma", fdr = 0.2, lfc = 0.6)
  
}
