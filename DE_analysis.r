library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)

# load the counts
x_AML <- read.delim("/Volumes/Backup Plus/archs4_retrieval/counts_AML.tsv", header=TRUE, sep="\t")
x_RAN <- read.delim("/Volumes/Backup Plus/archs4_retrieval/counts_RAN.tsv", header=TRUE, sep="\t")
x_RAN <- x_RAN[,2:1001]

# concat the samples together
x_all <- cbind(x_AML,x_RAN)

# create an edgeR DGElist object
d0 <- DGEList(counts= x_all[,2:1679], genes= x_all[,1])

# calculate normalising factors for each sample: 
# in some samples certain genes will be oversampled
# so the other genes will appear downregulated;
# but they're just downsamples -> taking this into acc
d0 <- calcNormFactors(d0)

# not taking into account lowly expressed transcripts
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]

# how many transcripts left
dim(d)

# extracting sample names
snames <- colnames(x_all[,2:1679])
group <- substr(snames, 1, 3)

# creating a 'model' matrix indicating the experiment design
# in this case we only have 2 groups of samples: AML and random
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)

# limma's linear fit
fit <- lmFit(y, mm)

# Bayes smoothening
fit_B <- eBayes(fit)
res <- decideTests(fit_B)

head(coef(fit))

# find the most contrasted transcripts using Bayes model again
contr <- makeContrasts(groupRAN - groupAML, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)

# top DE transcripts, sorted by p-value
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 50)

glXYPlot(x=tmp$coefficients, y= tmp$lods, xlab="logFC", ylab="logodds")

#qvals <- p.adjust(fit$p.value[,2], method = 'BH')