library(tidyverse)
library(pheatmap)

p <- read.table('phenotype.txt',header=T)
g <- read.table('./fl_genotype.txt',header=T)
sample <- p$IID
g <- g[,c('ID',sample)]
rownames(g) <- g$ID
g_matrix <- as.matrix(g[,2:377])

# mean(p$UHML)
# p1 <- p %>% arrange(UHML)
# mean(head(p1,50)$UHML)

p1 <- pheatmap(g_matrix,filename='376_genotype_heatmap.pdf',width=80,height=10,show_rownames = F)
sample_order <- colnames(g_matrix)[p1$tree_col[["order"]]]
fl <- p[match(sample_order,p$IID),]$UHML
fl <- as.matrix(fl)
rownames(fl) <- sample_order

# test 2
long_sample <- c('S302','S279','S003')
median_sample <- c('S379','S348','S397')
short_sample <- c('S293','S130','S360')

l <- fl[long_sample,]
m <- fl[median_sample,]
s <- fl[short_sample,]

# heatmap
options(repr.plot.width=3.5, repr.plot.height=3)
g_select <- g_matrix[,c(long_sample,median_sample,short_sample)]
apply(g_select,2,function(x)sum(x,na.rm = T))

color <- c('#4dbbd5','#4ab598','#fec551')

# 调整分支顺序
# https://blog.csdn.net/woodcorpse/article/details/109733806
exprTable_t <- as.data.frame(t(g_select))
col_dist = dist(exprTable_t)
hclust_1 <- hclust(col_dist)
manual_order = c("S360", "S130", "S293",'S397','S348','S379', "S279", "S302",  "S003")
dend = reorder(as.dendrogram(hclust_1), wts=order(match(manual_order, rownames(exprTable_t))))
col_cluster <- as.hclust(dend) 

annotation_col = data.frame(
                    SampleType = rep(c("short", "median","long"), each=3)
                )
rownames(annotation_col) = manual_order

pheatmap(g_select, cluster_cols = col_cluster, show_rownames = F,legend = F,
        color = color,annotation_col = annotation_col,cutree_cols = 3)


pheatmap(g_select, cluster_cols = col_cluster, show_rownames = F,legend = F,color = color,annotation_col = annotation_col, ,cutree_cols = 3,filename = 'sample_breeding_1.pdf',width=3.5,height=3)

mean(s)
mean(m)
mean(l)
