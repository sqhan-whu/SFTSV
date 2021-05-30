library('edgeR')
library('limma')
library('locfit')
library('splines')
library('KernSmooth')
library('statmod')
library('MASS')
library('dplyr')
library('data.table')
library('reshape2')
library('DESeq2')
#Set working directory
mainDir = '' #CHANGE this to reflect the main directory you've selected
setwd(mainDir)

#Get raw data
rawdata = read.table("gene_expression_4_vs7.txt",sep='\t',header=TRUE)
group = factor(c(rep(1,4), rep(2,7)))  #1 = negative, 2 = positive
y = DGEList(counts=rawdata[,2:ncol(rawdata)], genes=rawdata[,1],group=group)
keep = rowSums(cpm(y)>1) >= 2
y = y[keep,]
y$samples$lib.size = colSums(y$counts)

delivery = group
data.frame(Sample=colnames(y),delivery)

design = model.matrix(~0+delivery)
rownames(design) = colnames(y)

y = estimateGLMCommonDisp(y, design, verbose=TRUE)
y = estimateGLMTrendedDisp(y,design)
y = estimateGLMTagwiseDisp(y,design)

#Exact Test
et = exactTest(y)
et_tags = topTags(et, n=Inf)$table %>% filter(PValue < 0.001 & abs(logFC) >1)

#GLM
fit = glmFit(y, design)
contrast.design= c(1,-1)
lrt = glmLRT(fit, contrast=contrast.design)
glm_tags = topTags(lrt, n = Inf)$table %>% filter(PValue < 0.001 & abs(logFC) >1)

#QL F-test
fit = glmQLFit(y, design)
contrast.design=c(1,-1)
qlf = glmQLFTest(fit, contrast=contrast.design)
qlf_tags = topTags(qlf, n=Inf)$table  %>% filter(PValue < 0.001 & abs(logFC) >1)

#Merge tests
top_genes = Reduce(intersect, list(et_tags$genes, glm_tags$genes, qlf_tags$genes))
length(top_genes)




#Plot hierachical clustering of samples 

colData <- data.frame(condition = group,row.names = colnames(rawdata)[2:ncol(rawdata)])
dds <- DESeqDataSetFromMatrix(countData = rawdata[,2:ncol(rawdata)], 
                              colData = colData,
                              design = ~ condition )
dds <- dds[ rowSums(counts(dds)) > 1, ]

rld <- rlog(dds, blind = FALSE)
#rld <- vst(dds, blind = FALSE)
sampleDists <- dist(t(assay(rld)))
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf('fig_3a.pdf')
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()



#Normalize data
rawdata_norm = as.data.frame(t(t(rawdata[,-1])/colSums(rawdata[,-1])*1000000))
data= cbind(rawdata$Geneid, rawdata_norm)
#Filter data and save only genes of interest
idx = match(top_genes, rawdata$Geneid)
data_filtered= cbind(rawdata$Geneid[idx], rawdata_norm[idx,])
names(data_filtered)[names(data_filtered)=='id'] = 'genes'
write.csv(data_filtered, 'gene_expression_4_vs7.diff955.csv', row.names = FALSE)


#Plot heatmap of difference expression genes

library(pheatmap)
library(RColorBrewer)

pdf("fig_3b.pdf",width=8, height=10)
data=read.csv("gene_expression_4_vs7.diff955.csv",head=T,row.names=1)
data=log(data+1,2)
tmp=data.frame(group=c(rep('negative',4),rep('positive',7))) 
rownames(tmp)=colnames(data) 

pheatmap(data,scale = "row",cluster_cols = TRUE,
          cluster_rows = TRUE,
                   color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(seq(-2, 2, by = .1))),
                   breaks = seq(-2, 2, by = .1),annotation_col =tmp ,
show_rownames = F,
                   )#gaps_col = c(4,9))
dev.off()



#Plot Difference gene-expression in SFTSV virus (positive) and SFTSV virus (negative) with P < 0.001,exact test, likelihood ratio test, and quasi-likehood F test. 

total_table <- merge(merge(topTags(et, n=Inf)$table,topTags(lrt, n = Inf)$table,by='genes'),topTags(qlf, n = Inf)$table,by='genes')
detable=total_table
detable<- data.table(detable)
pal<- colorRampPalette("transparent", space = "Lab") # Do not colour NS genes
pdf('fig_3c.pdf', w= 12/2.54, h= 12/2.54, pointsize= 10)
par(las= 1, mgp= c(1.75, 0.5, 0), bty= 'l', mar= c(3, 3, 3, 0.5))
smoothScatter(x= detable$logCPM, y= detable$logFC, xlab= 'logCPM', ylab= 'logFC',
    main= "Differential gene expression [SFTSV - Healthy]", colramp= pal, col= 'blue', nrpoints= 0.1)
#lines(loess.smooth(x= detable$logCPM, y= detable$logFC, span= 0.1), lwd= 2, col= 'grey60')
abline(h= c( 0), col= 'grey30')
points(x= detable$logCPM, y= detable$logFC, col= ifelse(detable$PValue.x < 0.001 & detable$PValue.y < 0.001 & detable$PValue < 0.001, '#FF000080', 'grey75'), cex= 2, pch= '.') # Mark DE genes
mtext(side= 3, line= -1.2, text= sprintf('PValue < 0.001: %s', nrow(detable[PValue.x < 0.001 & PValue.y < 0.001 & PValue < 0.001 & logFC > 0])), adj= 1)
mtext(side= 1, line= -1.2, text= sprintf('PValue < 0.001: %s', nrow(detable[PValue.x < 0.001 & PValue.y < 0.001 & PValue < 0.001 & logFC < 0])), adj= 1)
grid(col= 'grey50')
dev.off()

write.csv(total_table, 'final.total.fold.csv', row.names = FALSE)

