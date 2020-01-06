# R code adapted from cummeRbund vignette source
library(cummeRbund)
cuff <- readCufflinks() # in a CUFFDIFF directory, not cufflinks!
cuff # to display it

conditions.vec <- unique(cummeRbund::samples(cuff)[,2])
runDetails <- runInfo(cuff)
infoTable   <- replicates(cuff)

# ==============================


pdf("fig-heatmap-combined.pdf")
h<-csHeatmap(myGenes,cluster='both', replicates=FALSE, labRow=TRUE)
h
dev.off()

pdf("fig-heatmap-replicates.pdf")
h.rep<-csHeatmap(myGenes,cluster='both',replicates=TRUE, labRow=TRUE)
h.rep
dev.off()

# ==============================
pdf("fig-expression-barplot.pdf")
b<-expressionBarplot(myGenes)
b
print(b)
dev.off()

# ==============================



###################################################
### code chunk number 52: gene_level_1
###################################################
myGeneId <- "ENSMUSG00000039376"

myGene<-getGene(cuff,myGeneId)
myGene
head(fpkm(myGene))
head(fpkm(isoforms(myGene)))


###################################################
### code chunk number 53: gene_plots_1
###################################################
pdf("fig-single-gene-plot-expression.pdf")
gl<-expressionPlot(myGene)
gl
dev.off()

pdf("fig-single-gene-plot-expression-reps.pdf")
gl.rep<-expressionPlot(myGene,replicates=TRUE)
gl.rep
dev.off()

pdf("fig-single-gene-plot-isoforms.pdf")
gl.iso.rep<-expressionPlot(isoforms(myGene),replicates=T)
gl.iso.rep
dev.off()

pdf("fig-single-gene-plot-isoforms-reps.pdf")
gl.cds.rep<-expressionPlot(CDS(myGene),replicates=T)
gl.cds.rep
dev.off()


###################################################
### code chunk number 58: gene_plots_2
###################################################
pdf("fig-single-gene-plot-bar.pdf")
gb<-expressionBarplot(myGene)
gb
dev.off()

pdf("fig-single-gene-plot-bar-rep.pdf")
gb.rep<-expressionBarplot(myGene,replicates=T)
gb.rep
dev.off()


###################################################
### code chunk number 62: gene_plots_bar_isoforms
###################################################
pdf("fig-single-gene-plot-bar-isoform.pdf")
igb<-expressionBarplot(isoforms(myGene),replicates=T)
igb
print(igb)
dev.off()

###################################################
### code chunk number 63: features_1
###################################################
head(features(myGene))

###################################################
### code chunk number 64: features_2
###################################################
#pdf("fig-single-gene-track.pdf")
#genetrack<-makeGeneRegionTrack(myGene)
#plotTracks(genetrack)
#dev.off()

pdf("fig-dist-heat.pdf");
myDistHeat<-csDistHeat(genes(cuff))
print(myDistHeat)
dev.off()

pdf("fig-dist-heat-reps.pdf");
myDistHeat<-csDistHeat(genes(cuff), replicates=TRUE)
print(myDistHeat)
dev.off()

###################################################
### code chunk number 76: dim_reduction_1
###################################################

pdf("fig-pca-1.pdf")
genes.PCA<-PCAplot(genes(cuff),"PC1","PC2")
genes.PCA
dev.off()

pdf("fig-mds-1.pdf")
genes.MDS<-MDSplot(genes(cuff)) # this seems to often fail with a "cmdscale" error
genes.MDS
dev.off()

pdf("fig-pca-rep.pdf")
genes.PCA.rep<-PCAplot(genes(cuff),"PC1","PC2",replicates=T)
genes.PCA.rep
dev.off()

pdf("fig-mds-rep.pdf")
genes.MDS.rep<-MDSplot(genes(cuff),replicates=T)
genes.MDS.rep
dev.off()


###################################################
ic<-csCluster(myGenes,k=4)
head(ic$cluster)

pdf("fig-icp.pdf")
icp<-csClusterPlot(ic)
icp
dev.off()

###################################################
### code chunk number 83: specificity_1
###################################################
myGenes.spec<-csSpecificity(myGenes)
head(myGenes.spec)

print("Note: similarity search is VERY SLOW, but is super cool")
mySimilar<-findSimilar(cuff, myGeneId,n=20)
pdf("fig-similar.pdf")
mySimilar.expression<-expressionPlot(mySimilar,logMode=T,showErrorbars=F)
mySimilar.expression
dev.off()

###################################################
### code chunk number 86: similar_2
###################################################
myProfile<-c(500,0,400)
mySimilar2<-findSimilar(cuff,myProfile,n=10)
mySimilar2.expression<-expressionPlot(mySimilar2,logMode=T,showErrorbars=F)


###################################################
### code chunk number 87: similar_plots_2
###################################################
myProfile<-c(500,0,400)
mySimilar2<-findSimilar(cuff,myProfile,n=10)
mySimilar2.expression<-expressionPlot(mySimilar2,logMode=T,showErrorbars=F)
print(mySimilar2.expression)

# ==============================

# "SCV visualization"

pdf("fig-genes-scv.pdf");
genes.scv<-fpkmSCVPlot(genes(cuff))
genes.scv
dev.off()

pdf("fig-isoforms-scv.pdf");
isoforms.scv<-fpkmSCVPlot(isoforms(cuff))  # useless???
isoforms.scv
dev.off()

