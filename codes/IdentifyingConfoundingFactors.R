###
#   File name : IdentifyingConfoundingFactors.R
#   Author    : Hyunjin Kim
#   Date      : Feb 6, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Run PCA & TSNE to find possible confounding factors
#
#   Instruction
#               1. Source("IdentifyingConfoundingFactors.R")
#               2. Run the function "identifyCF" - specify an input file path (raw counts) and output directory
#               3. The PCA plots will be generated under the output directory
#
#   Example
#               > source("The_directory_of_IdentifyingConfoundingFactors.R/IdentifyingConfoundingFactors.R")
#               > identifyCF(rawCntPath="./data/raw_counts.rda",
#                            outputDir="./results/pca/")
###

identifyCF <- function(rawCntPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Ferrari/RNASeq_Analysis/data/raw_counts.rda",
                       outputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Ferrari/RNASeq_Analysis/results/pca") {
  
  ### load dataset
  load(rawCntPath)
  
  ### remove duplicated gene symbols
  rawCnt <- rawCnt[!duplicated(rawCnt$Gene_Symbol),]
  rownames(rawCnt) <- rawCnt$Gene_Symbol
  
  ### A function to transform RNA-Seq data with VST in DESeq2 package
  normalizeRNASEQwithVST <- function(readCount, filtering=TRUE) {
    
    ### load library
    if(!require(DESeq2)) {
      source("https://bioconductor.org/biocLite.R")
      biocLite("DESeq2")
      library(DESeq2)
    }
    
    ### make a design matrix for DESeq2 data
    condition <- data.frame(factor(rep("OneClass", ncol(readCount))))
    
    ### Data preparation for DESeq2 format
    deSeqData <- DESeqDataSetFromMatrix(countData=readCount, colData=condition, design= ~0)
    
    if(filtering == TRUE) {
      ### Remove rubbish rows - this will decrease the number of rows
      deSeqData <- deSeqData[rowSums(counts(deSeqData))>1,]
    }
    
    ### VST
    vsd <- vst(deSeqData)
    transCnt <- data.frame(assay(vsd), check.names = FALSE)
    
    return (transCnt)
    
  }
  
  ### A function to perform 2D PCA and save a plot
  ### normalizedMat: rows are genes and columns are samples
  ### grp: group information of the samples
  ### num: the number of top genes to be used based on variance (-1 [default]: use all the genes)
  ### component: to draw a plot with PC1 & PC2 or PC2 & PC3
  ### title: title of the plot
  ### outDir: output directory for the plot
  pca_plot <- function(normalizedMat, grp,
                       num = -1, component=c("PC1&PC2", "PC2&PC3"),
                       title="PCA_Plot", outDir="./") {
    
    ### load library
    if(!require(ggfortify, quietly = TRUE)) {
      install.packages("ggfortify")
      library(ggfortify, quietly = TRUE)
    }
    if(!require(FactoMineR, quietly = TRUE)) {
      install.packages("FactoMineR")
      library(FactoMineR, quietly = TRUE)
    }
    if(!require(factoextra, quietly = TRUE)) {
      install.packages("factoextra")
      library(factoextra, quietly = TRUE)
    }
    if(!require(xlsx, quietly = TRUE)) {
      install.packages("xlsx")
      require(xlsx, quietly = TRUE)
    }
    
    ### select the top genes based on variance
    if(num >= 0 && num <= nrow(normalizedMat)) {
      v <- apply(normalizedMat, 1, var)
      v <- v[order(-v)]
      top_genes <- names(v)[1:num]
    } else {
      top_genes <- rownames(normalizedMat)
    }
    
    ### PCA
    pca_result <- PCA(t(normalizedMat[top_genes,]), graph = FALSE)
    colnames(pca_result$ind$coord) <- paste0("PC", 1:ncol(pca_result$ind$coord))
    colnames(pca_result$var$contrib) <- paste0("PC", 1:ncol(pca_result$var$contrib))
    pca_group <- data.frame(pca_result$ind$coord, group=grp)
    
    colors = topo.colors(length(unique(grp)))
    names(colors) = unique(grp)
    
    ### save as png
    if(component[1] == "PC1&PC2") {
      ggplot(pca_group,aes(x=PC1,y=PC2,col=group)) +
        labs(title=paste0(title, "_PC1-2")) +
        geom_text(aes(label=colnames(normalizedMat)),hjust="inward", vjust="inward") +
        scale_color_manual(values = colors) +
        theme_classic(base_size = 16)
      ggsave(filename = paste0(outDir, title, "_PC1-2", ".png"), width = 10, height = 8)
      
      fviz_contrib(pca_result, choice = "var", axes = 1, top = 30)
      ggsave(filename = paste0(outDir, title, "_PC1_contribution.png"), width = 12, height = 8)
      
      fviz_contrib(pca_result, choice = "var", axes = 2, top = 30)
      ggsave(filename = paste0(outDir, title, "_PC2_contribution.png"), width = 12, height = 8)
    } else if(component[1] == "PC2&PC3") {
      ggplot(pca_group,aes(x=PC2,y=PC3,col=group)) +
        labs(title=paste0(title, "_PC2-3")) +
        geom_text(aes(label=colnames(normalizedMat)),hjust="inward", vjust="inward") +
        scale_color_manual(values = colors) +
        theme_classic(base_size = 16)
      ggsave(filename = paste0(outDir, title, "_PC2-3", ".png"), width = 10, height = 8)
      
      fviz_contrib(pca_result, choice = "var", axes = 2, top = 30)
      ggsave(filename = paste0(outDir, title, "_PC2_contribution.png"), width = 12, height = 8)
      
      fviz_contrib(pca_result, choice = "var", axes = 3, top = 30)
      ggsave(filename = paste0(outDir, title, "_PC3_contribution.png"), width = 12, height = 8)
    } else {
      stop("\"component\" parameter should be \"PC1&PC2\" or \"PC2&PC3\"")
    }
    
    write.xlsx2(data.frame(Gene_Symbol=rownames(pca_result$var$contrib), pca_result$var$contrib,
                           stringsAsFactors = FALSE, check.names = FALSE),
                file = paste0(outDir, title, "_PC_contribution.xlsx"),
                sheetName = "PCA_contribution", row.names = FALSE)
    
  }
  
  ### A function to perform t-SNE and save a plot
  ### normalizedMat: rows are genes and columns are samples
  ### grp: group information of the samples
  ### showNames: if TRUE, sample names will be presented in the plot
  ###            if FALSE, only dots will be presented (default: FALSE)
  ### num: the number of top genes to be used based on variance (-1 [default]: use all the genes)
  ### title: title of the plot
  ### outDir: output directory for the plot
  tsne_plot <- function(normalizedMat, grp, showNames = FALSE, num = -1,
                        title="TSNE_Plot", outDir="./") {
    ### load library
    if(!require(Rtsne)) {
      install.packages("Rtsne")
      library(Rtsne)
    }
    
    ### select the top genes based on variance
    if(num >= 0 && num <= nrow(normalizedMat)) {
      v <- apply(normalizedMat, 1, var)
      v <- v[order(-v)]
      top_genes <- names(v)[1:num]
    } else {
      top_genes <- rownames(normalizedMat)
    }
    
    ### TSNE
    set.seed(1234)
    tryCatch({
      writeLines("Perplexity = 30")
      t <- Rtsne(t(normalizedMat[top_genes,]), perplexity = 30)
    }, error = function(err) {
      tryCatch({
        writeLines("Perplexity = 10")
        t <- Rtsne(t(normalizedMat[top_genes,]), perplexity = 10)
      }, error = function(err) {
        tryCatch({
          writeLines("Perplexity = 5")
          t <- Rtsne(t(normalizedMat[top_genes,]), perplexity = 5)
        }, error = function(err) {
          tryCatch({
            writeLines("Perplexity = 3")
            t <- Rtsne(t(normalizedMat[top_genes,]), perplexity = 3)
          }, error = function(err) {
            tryCatch({
              writeLines("Perplexity = 2")
              t <- Rtsne(t(normalizedMat[top_genes,]), perplexity = 2)
            }, error = function(err) {
              writeLines("Perplexity = 1")
              t <- Rtsne(t(normalizedMat[top_genes,]), perplexity = 1)
            })
          })
        })
      })
    })
    
    colors = rainbow(length(unique(grp)))
    names(colors) = unique(grp)
    
    ### save a plot
    png(filename=paste0(outDir, title, ".png"), width = 1000, height = 800)
    if(showNames) {
      plot(t$Y, col="white", xlab="tsne1", ylab="tsne2", main = title)
      text(t$Y, labels=colnames(normalizedMat), col=colors[grp])
    } else {
      plot(t$Y, col=colors[grp], pch=19, xlab="tsne1", ylab="tsne2", main = title)
    }
    legend("topright", legend = unique(grp), col = colors[unique(grp)], pch = 15)
    dev.off()
  }
  
  ### make a pca plot
  pca_plot(normalizeRNASEQwithVST(rawCnt[,-c(1,2)]), grp = sample_info, num = 1000,
           title = "PCA_plot_AV", outDir = outputDir)
  
  ### make a tsne plot
  tsne_plot(normalizeRNASEQwithVST(rawCnt[,-c(1,2)]), grp = sample_info, num = 1000,
            showNames = TRUE, title = "TSNE_plot_AV", outDir = outputDir)
  
}
