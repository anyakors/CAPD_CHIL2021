#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(c("rhdf5", "dendextend", "sva", "edgeR", "preprocessCore"))

#install.packages(c("tools", "corrplot", "colorspace", "data.table"))

library("tools")#
library("dendextend")#
library("preprocessCore")#
library("sva")
library("edgeR")
library("corrplot")
library("colorspace")
library("data.table")

main <- function() {

  args <- commandArgs(trailingOnly = TRUE)
  dir.create(file.path(args[2]), showWarnings = FALSE)

  filename <- args[1]
  #res <- read.csv(file = filename)
  #r <- res[,1]#
  #res <- res[1:nrow(res),2:ncol(res)]
  #res <- t(res)
  #colnames(res) <- r

  res <- fread(filename)
  r <- res[,1]
  r <- t(r)
  res <- res[1:nrow(res),2:ncol(res)]
  res <- t(res)
  colnames(res) <- r
  cat("Imported counts")
  cat("\n")

  #jpeg(file.path(args[2], "plots", substr(r[1], 1, 3), "boxplot_before_norm.jpg"))
  #boxplot(log2(1+res), main="read count distribution by sample")
  #dev.off()

  # NORMALIZE

  res <- cpm(res)
  #cat("Summary per sample:", summary(res), sep = "\n")

  f <- calcNormFactors(res, method="TMM")
  res <- sweep(res, MARGIN=2, f, `*`)

  # CALCULATE CORR, CLUSTER AND PLOT DENDROGRAM

  cc = cor(res)

  dir.create(file.path(args[2], "plots"), showWarnings = FALSE)
  dir.create(file.path(args[2], "plots", substr(r[1], 1, 3)), showWarnings = FALSE)
  jpeg(file.path(args[2], "plots", substr(r[1], 1, 3), "corrplot.jpg"))
  corrplot(cc, tl.col = "black", order = "hclust", hclust.method = "average", addrect = 2, tl.cex = 0.5)
  dev.off()

  dend <- as.dendrogram(hclust(as.dist(1-cc)))
  useries = unique(r)
  series_match = useries[match(r, useries)]
  colos <- colorspace::rainbow_hcl(length(useries), c=160, l=50)
  names(colos) = useries
  series_color <- colos[series_match]

  clu = cutree(dend, h=0.25)
  labels_colors(dend) <- series_color[order.dendrogram(dend)]
  dend <- color_branches(dend, h=0.25)

  jpeg(file.path(args[2], "plots", substr(r[1], 1, 3), "dendrogram.jpg"))
  par(mar=c(4,1,1,12), cex=0.6)
  plot(dend, horiz=TRUE)
  colored_bars(cbind(clu, series_color), dend, rowLabels = c("Cluster", "Sample"), horiz=TRUE, y_scale=0.05, text_shift=0.1)
  legend("topleft", legend=useries, fill=colos, bg="white", cex=0.6)
  dev.off()

  meancc = apply(cc,2,mean)
  sdcc = sd(meancc)
  numbersd = (meancc-mean(meancc))/sdcc

  #sdout = -2
  sdout = as.numeric(args[3])
  cat(paste("Outlier threshold:", sdout, sep=" "))
  cat("\n")

  jpeg(file.path(args[2], "plots", substr(r[1], 1, 3), "outliers.jpg"))
  plot(numbersd)
  text(1:ncol(cc), numbersd, labels=colnames(cc), cex=0.9, pos=3)
  abline(h=sdout) 
  dev.off()

  outliers = dimnames(res)[[2]][numbersd<sdout] 
  cat(paste("Dropped", length(outliers), "outliers:", outliers, sep=" "))
  cat("\n")
  res <- res[ , -which(colnames(res) %in% outliers)]
  res <- t(res)
  write.csv(res, file.path(args[2], paste(substr(r[1], 1, 3), ".csv", sep="")))

  cat("Saved normalized .csv count table in the out dir")
  cat("\n")
  
}

main()