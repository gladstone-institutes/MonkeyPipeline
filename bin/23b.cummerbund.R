# You need cummerbund to run this stuff
# source("http://bioconductor.org/biocLite.R"); biocLite("cummeRbund")

# You also need Cuffdiff input files (the whole directory!)
#myGenes<-getGenes(cuff,myGeneIds)
#hMap <- csHeatmap(myGenes,cluster='both')
# See docs at: http://www.bioconductor.org/packages/release/bioc/vignettes/cummeRbund/inst/doc/cummeRbund-manual.pdf
# And a "walkthrough" at: http://compbio.mit.edu/cummeRbund/manual_2_0.html

# R code adapted from cummeRbund vignette source

require("cummeRbund")

CUFFDIFF_DIR <- "./"

PNG_DEFAULT_SIZE  <- 3840
PNG_SIZE_PER_ITEM <- 500 # in pixels
PNGMARGIN         <- 400 # default margin for each figure
PNGRES            <- 300 # 300 dpi is a good PNG resolution

setwd(CUFFDIFF_DIR)

if (!exists("cuff")) { cuff <- readCufflinks() }

conditions.vec  <- unique(cummeRbund::samples(cuff)[,2])
runDetails      <- runInfo(cuff)
infoTable       <- replicates(cuff)
totalNumGroups  <- ncol(genes(cuff))      # Number of distinct experimental groups
totalNumSamples <- nrow(replicates(cuff)) # Total number of samples IN ALL, counting all replicates separately

gene.features<-annotation(genes(cuff))
head(gene.features)

gene.fpkm<-fpkm(genes(cuff))
head(gene.fpkm)

gene.repFpkm<-repFpkm(genes(cuff))
head(gene.repFpkm)

gene.counts<-count(genes(cuff))
head(gene.counts)

isoform.fpkm<-fpkm(isoforms(cuff))
head(isoform.fpkm)

gene.diff<-diffData(genes(cuff))
head(gene.diff)

sample.names<-samples(genes(cuff))
head(sample.names)
gene.featurenames<-featureNames(genes(cuff))
head(gene.featurenames)

gene.matrix<-fpkmMatrix(genes(cuff))
head(gene.matrix)

###################################################
### code chunk number 34: data_access_4
###################################################

print("Making a replicate FPKM matrix...")
gene.rep.matrix<-repFpkmMatrix(genes(cuff))
head(gene.rep.matrix)

###################################################
### code chunk number 35: data_access_5
###################################################
print("Making a COUNT matrix...")
gene.count.matrix<-countMatrix(genes(cuff))
head(gene.count.matrix)

# =================================================

sigs.ids <- cummeRbund::getSig(cuff,alpha=0.05,level='genes') # Significant GENE IDs at cutoff 0.001
# "all genes that reject the null hypothesis in any condition tested"
# mySigTable<-getSigTable(cuff,alpha=0.01,level='genes')
print(paste("Num significant genes at this cutoff is ", length(sigs.ids), sep=''))

sigs <- cummeRbund::getGenes(cuff, head(sigs.ids, n=5))
warning("Only picking the top 5 genes for testing here")

#Isoform-level FPKMs for gene set
head(fpkm(isoforms(sigs)))

#Replicate FPKMs for TSS groups within gene set
head(repFpkm(TSS(sigs)))

print("Replicates are:")
print(replicates(cuff))

print("Option: can add annotation at this step, but we don't right now.")
#annot<-read.table("gene_annotation.tab",sep="\t",header=T,na.string="-")
#addFeatures(cuff,annot,level="genes")

head(fpkm(sigs))

# Show the ISOFORMS of these genes
pdf("test-a.pdf"); igb<-expressionBarplot(isoforms(sigs),replicates=T); print(igb); dev.off()

# Show the genes
pdf("test-ba.pdf"); igb<-expressionBarplot(sigs,replicates=T); print(igb); dev.off()

# get annotation with:
annotation(sigs)

# getGene or getGenes gets a gene object from the ID
pdf("test-expr.pdf")
gl.rep<-expressionPlot(sigs,replicates=TRUE)
print(gl.rep)
dev.off()

pdf("test-dend.pdf")
den<-csDendro(sigs)
print(den)
dev.off()

head(fpkm(sigs))
head(fpkm(isoforms(sigs)))

# heatmap

pdf("test-c.pdf");
ih<-csHeatmap(isoforms(sigs),cluster='both',labRow=F)
print(ih)
dev.off()

pdf("test-d.pdf");
ih<-csHeatmap(sigs,cluster='both',labRow=F, replicates=TRUE)
print(ih)
dev.off()

pdf("test-cx.pdf");
th<-csHeatmap(TSS(sigs),cluster='both',labRow=F)
print(th)
dev.off()

pdf("test-cx3.pdf");
ch<-csHeatmap(CDS(sigs),cluster='both',labRow=F)
print(ch);
dev.off()

pdf("test-cxdist3.pdf");
myRepDistHeat <- csDistHeat(sigs,replicates=T)
print(myRepDistHeat)
dev.off()

############## DONE WITH SIGNIFICANT-ONLY GENES ###############
# ==============================================================================================
print("---------------------------------------------------------------------------------------")
print("Now plotting all the 'general' plots for ALL genes/isoforms, not just significant ones.")
print("---------------------------------------------------------------------------------------")

# =================================================

disp<-dispersionPlot(genes(cuff))
print("Generating the dispersion scatterplots...")
pdf("Cuffdiff_Fig_01_Dispersion_Scatterplots.pdf", width=20, height=20)
print(disp)
dev.off()
png("Cuffdiff_Fig_01_Dispersion_Scatterplots.png", width=PNG_DEFAULT_SIZE, height=PNG_DEFAULT_SIZE, res=PNGRES)
print(disp) # PNG version of the (large) PDF above
dev.off()

# =================================================

FIG_2_INCH_SIZE <- 20
print("2a: Gene COV vs FPKM...")
pdf("Cuffdiff_Fig_02a_Gene_COV_vs_FPKM.pdf", width=FIG_2_INCH_SIZE, height=FIG_2_INCH_SIZE)
genes.scv    <- fpkmSCVPlot(genes(cuff))
print(genes.scv)
dev.off()

print("2b: Isoform COV vs FPKM...")
pdf("Cuffdiff_Fig_02b_Isoform_COV_vs_FPKM.pdf", width=FIG_2_INCH_SIZE, height=FIG_2_INCH_SIZE)
isoforms.scv <- fpkmSCVPlot(isoforms(cuff))
print(isoforms.scv)
dev.off()

# =================================================

print("3a: FPKM Density...")
pdf("Cuffdiff_Fig_03a_Gene_FPKM_Density_per_group.pdf", width=18, height=18)
dens <- cummeRbund::csDensity(genes(cuff))
print(dens)
dev.off()

densRep <- cummeRbund::csDensity(genes(cuff),replicates=T)
print("3b: FPKM Density, showing each individual sample...")
pdf("Cuffdiff_Fig_03b_Gene_FPKM_Density_per_sample.pdf", width=18, height=18)
print(densRep)
dev.off()

# =================================================

print("4a")
pdf("Cuffdiff_Fig_04a_Gene_FPKM_Boxplot_per_group.pdf", width=18, height=18)
box <- cummeRbund::csBoxplot(genes(cuff))
print(box)
dev.off()

print("4b")
pdf("Cuffdiff_Fig_04b_Gene_FPKM_Boxplot_per_sample.pdf", width=18, height=18)
boxRep <- cummeRbund::csBoxplot(genes(cuff), replicates=T)
print(boxRep)
dev.off()

# =================================================

print("SKIPPING the scatterplot calculations.")
if (FALSE) {
	print("5 - scatterplot per-group - this is slow... (5+ minutes per group)")
	png("Cuffdiff_Fig_05a_Gene_FPKM_Scatterplot_per_group.png", width=(PNGMARGIN + totalNumGroups*PNG_SIZE_PER_ITEM), height=(PNGMARGIN + totalNumGroups*PNG_SIZE_PER_ITEM), res=PNGRES)
	sMat <-csScatterMatrix(genes(cuff), replicates=FALSE, hexbin=FALSE)
	print(sMat)
	dev.off()
	print("5 - scatterplot per-sample - this is slow... (5+ minutes per replicate)")
	png("Cuffdiff_Fig_05b_Gene_FPKM_Scatterplot_per_sample.png", width=(PNGMARGIN + totalNumSamples*PNG_SIZE_PER_ITEM), height=(PNGMARGIN + totalNumSamples*PNG_SIZE_PER_ITEM), res=PNGRES)
	sMatRep <-csScatterMatrix(genes(cuff), replicates=TRUE, hexbin=FALSE)
	print(sMatRep)
	dev.off()

        pdf("fig-scatter-first-two-groups-only.pdf");
        s2<-csScatter(genes(cuff), conditions.vec[1], conditions.vec[2],smooth=T)
        print(s2)
        dev.off()
}
# =================================================

print("6 Volcano plot (slow-ish, ~3 minutes)")
vol <- cummeRbund::csVolcanoMatrix(genes(cuff), showSignificant=TRUE, alpha=0.01)
png("Cuffdiff_Fig_06_Volcano_plot_Matrix_red_pval_cutoff_0_01.png", width=(PNGMARGIN + totalNumGroups*PNG_SIZE_PER_ITEM), height=(PNGMARGIN + totalNumGroups*PNG_SIZE_PER_ITEM), res=PNGRES)
print(vol)
dev.off()

if (FALSE) {
     print("Volcano plot between two SPECIFIC conditions!")
     pdf("fig-volcano-conditions.pdf");
     vcon<-csVolcanoMatrix(genes(cuff), conditions[1], conditions[2])
     print(vcon)
     dev.off()
}



# =================================================

PLOT_7_INCH_SIZE <- 20

print("7a Dendrogram per-group")
pdf("Cuffdiff_Fig_07a_Dendrogram_by_group.pdf", width=PLOT_7_INCH_SIZE, height=PLOT_7_INCH_SIZE)
dendro <- cummeRbund::csDendro(genes(cuff), replicates=FALSE)
print(dendro)
dev.off()
print("7b Dendrogram per-sample")
pdf("Cuffdiff_Fig_07b_Dendrogram_by_sample.pdf", width=PLOT_7_INCH_SIZE, height=PLOT_7_INCH_SIZE)
dendroRep <- cummeRbund::csDendro(genes(cuff), replicates=TRUE)
print(dendroRep)
dev.off()

# =================================================

PLOT_8_INCH_SIZE <- 20

print("8a MDS (multi-dimensional scaling) plot (similar to PCA)")
pdf("Cuffdiff_Fig_08a_MDS_Multi_Dimensional_Scaling_by_group.pdf", width=PLOT_8_INCH_SIZE, height=PLOT_8_INCH_SIZE)
genes.MDS <- cummeRbund::MDSplot(genes(cuff))
print(genes.MDS)
dev.off()

print("8b MDS (multi-dimensional scaling) plot (similar to PCA) by replicate")
pdf("Cuffdiff_Fig_08b_MDS_Multi_Dimensional_Scaling_by_sample.pdf", width=PLOT_8_INCH_SIZE, height=PLOT_8_INCH_SIZE)
genes.MDS.rep <- cummeRbund::MDSplot(genes(cuff), replicates=TRUE)
print(genes.MDS.rep)
dev.off()

print("8c PCA plot")
pdf("Cuffdiff_Fig_08c_PCA_by_group.pdf", width=PLOT_8_INCH_SIZE, height=PLOT_8_INCH_SIZE)
genes.PCA <- cummeRbund::PCAplot(genes(cuff), "PC1", "PC2")
print(genes.PCA)
dev.off()

print("8d PCA plot by replicate")
pdf("Cuffdiff_Fig_08d_PCA_by_sample.pdf", width=PLOT_8_INCH_SIZE, height=PLOT_8_INCH_SIZE)
genes.PCA <- cummeRbund::PCAplot(genes(cuff), "PC1", "PC2", replicates=TRUE)
print(genes.PCA)
dev.off()

# =================================================

DIST_MAT_INCH_SIZE <- 20

print("9a Distance Heatmap (by group)")
pdf("Cuffdiff_Fig_09a_Distance_Heatmap_by_group.pdf", width=DIST_MAT_INCH_SIZE, height=DIST_MAT_INCH_SIZE)
distHeat <- cummeRbund::csDistHeat(genes(cuff))
print(distHeat)
dev.off()

print("9b Distance Heatmap (by sample)")
pdf("Cuffdiff_Fig_09b_Distance_Heatmap_by_sample.pdf", width=DIST_MAT_INCH_SIZE, height=DIST_MAT_INCH_SIZE)
distHeatRep <- cummeRbund::csDistHeat(genes(cuff), replicates=TRUE)
print(distHeatRep)
dev.off()

# =================================================

SIG_MAT_INCH_SIZE <- 20
print("10a Significant Genes Matrix")
pdf("Cuffdiff_Fig_10a_Significant_Gene_Matrix_0.0500.pdf", width=SIG_MAT_INCH_SIZE, height=SIG_MAT_INCH_SIZE)
sig <- cummeRbund::sigMatrix(cuff, level='genes', alpha=0.05)
print(sig)
dev.off()

print("10b Significant Genes Matrix")
pdf("Cuffdiff_Fig_10b_Significant_Gene_Matrix_0.0100.pdf", width=SIG_MAT_INCH_SIZE, height=SIG_MAT_INCH_SIZE)
sig <- cummeRbund::sigMatrix(cuff, level='genes', alpha=0.01)
print(sig)
dev.off()

print("10c Significant Genes Matrix")
pdf("Cuffdiff_Fig_10c_Significant_Gene_Matrix_0.0050.pdf", width=SIG_MAT_INCH_SIZE, height=SIG_MAT_INCH_SIZE)
sig <- cummeRbund::sigMatrix(cuff, level='genes', alpha=0.005)
print(sig)
dev.off()

print("10d Significant Genes Matrix")
pdf("Cuffdiff_Fig_10d_Significant_Gene_Matrix_0.0010.pdf", width=SIG_MAT_INCH_SIZE, height=SIG_MAT_INCH_SIZE)
sig <- cummeRbund::sigMatrix(cuff, level='genes', alpha=0.001)
print(sig)
dev.off()

print("10e Significant Genes Matrix")
pdf("Cuffdiff_Fig_10e_Significant_Gene_Matrix_0.0005.pdf", width=SIG_MAT_INCH_SIZE, height=SIG_MAT_INCH_SIZE)
sig <- cummeRbund::sigMatrix(cuff, level='genes', alpha=0.0005)
print(sig)
dev.off()

print("10f Significant Genes Matrix")
pdf("Cuffdiff_Fig_10f_Significant_Gene_Matrix_0.0001.pdf", width=SIG_MAT_INCH_SIZE, height=SIG_MAT_INCH_SIZE)
sig <- cummeRbund::sigMatrix(cuff, level='genes', alpha=0.0001)
print(sig)
dev.off()

# =================================================
if (FALSE) {
     FIG_11_INCH_SIZE <- 20
     print("Making an MAPlot between two conditions. These are two ARBITRARY conditions!")
     pdf("fig-maplot.pdf")
     m<-MAplot(genes(cuff), conditions.vec[1], conditions.vec[2])
     print(m)
     dev.off()
     pdf("fig-maplot-count.pdf")
     mCount<-MAplot(genes(cuff), conditions.vec[1], conditions.vec[2], useCount=T)
     print(mCount)
     dev.off()
}

# =================================================

print("Generating the success file 'z.success.cummerbund'")
file.create("z.success.cummerbund", showWarnings=TRUE) # Generate a file to show that this ran successfully. DO NOT CHANGE THIS NAME without also changing the matching name in "23a.cummerbund.pl"

if (file.exists("cuffData.db")) {
     print("Now removing the (huge) 'cuffData.db' file, since cummeRbund is finished with it.")
     file.remove("cuffData.db")
} else {
     print("Weirdly, the 'cuffData.db' file did not appear to exist...")
}
