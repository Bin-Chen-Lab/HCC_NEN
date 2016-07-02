#visualize disease signatures using volcano plot

#install.packages("devtools")
#devtools::install_github("stephenturner/Tmisc")

# Load the Tmisc library
library(Tmisc)
library(calibrate)
library(DESeq)


cancer = "LIHC"
load(paste( cancer, "/", cancer, '_sig_deseq_final.RData', sep=''))

# Read the data directly from the Github Gist 
# https://gist.github.com/stephenturner/806e31fce55a8b7175af
#res <- read.gist("806e31fce55a8b7175af", header=TRUE)
#head(res)

# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pval), pch=20, main="Volcano plot", xlim=c(-8,15)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(res, padj<.0001 ), points(log2FoldChange, -log10(pval), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>2), points(log2FoldChange, -log10(pval), pch=20, col="orange"))
with(subset(res, padj<.0001 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pval), pch=20, col="green"))

# Label points with the textxy function from the calibrate plot
with(subset(res, padj<.005 & abs(log2FoldChange)>4), textxy(log2FoldChange, -log10(pval), labs=Gene, cex=.8))

res = subset(res, log2FoldChange != -Inf)
plotMA(res, col = ifelse( res$padj < 1E-1 & abs(res$log2FoldChange) > 2, "red", "gray"))
