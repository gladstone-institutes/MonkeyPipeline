#!/wynton/group/gladstone/third_party/Rscript  --no-save  --no-restore  --no-site-file  --no-init-file

require(methods) # <-- just in case we use RScript to run this

if (interactive()) { options(error=recover) } else { options(error=traceback) } # useful even non-interactively
options(stringsAsFactors=FALSE)

# Warning: single quotes are NOT ALLOWED on any of the PBS lines above! They trigger a "qsub: unmatched" error and exit code 1.

# This file gets, as input, a *NAIVE* FPKM matrix file, (not one with the only-expressed-genes values) from edgeR. Then it makes a bunch of figures using that information.
# That is why it is called "edgeRNaiveFpkmFile"

if (FALSE) {
     # Here are the commands to install everything!
     install.packages(c('data.table','gapmap','gplots','pheatmap','ggplot2'));
     source('http://bioconductor.org/biocLite.R'); biocLite(c('DESeq2','genefilter'));
}
#require(reshape) # reshape for 'melt'
if (!require(data.table))  { print("You need to: install.packages('data.table'); "); stopifnot(FALSE == "data.table was not installed"); }
if (!require(gapmap))  { print("You need to: install.packages('gapmap'); "); stopifnot(FALSE == "gapmap was not installed"); }
if (!require(gplots))  { print("You need to: install.packages('gplots'); "); stopifnot(FALSE == "gplots was not installed"); }
if (!require(pheatmap))   { print("You need to install 'pheatmap': install.packages('pheatmap'); "); stopifnot(FALSE == "pheatmap was not installed") }
if (!require(ggplot2))    { print("You need to install 'ggplot2': install.packages('ggplot2'); "); stopifnot(FALSE == "ggplot2 was not installed") }
if (!require(DESeq2))     { print("You need to install 'DESeq2': source('http://bioconductor.org/biocLite.R'); biocLite('DESeq2'); "); stopifnot(FALSE == "DESeq2 was not installed") }
# =================
rowVars.copied_from_genefilter <- function (x, ...) {
     n         <- rowSums(!is.na(x)) # Copied from genefilter to reduce the dependencies by 1: source('http://bioconductor.org/biocLite.R'); biocLite('genefilter')
     n[n <= 1] <- NA
     return(rowSums((x - rowMeans(x, ...))**2, ...)/(n-1))
}

saved.par.default <- par(no.readonly=TRUE)

looks_true                <- function(string) { return(!is.null(string) && !(toupper(string) %in% c("", "F","FALSE","UNDEF","0"))) }
print_stdout_and_stderr   <- function(...) { write(paste(..., sep=''), stderr()); write(paste(..., sep=''), stdout()); }
VERBOSE                   <- looks_true(Sys.getenv("verbose"))
FORCE                     <- looks_true(Sys.getenv("force"))
DEBUG                     <- looks_true(Sys.getenv("debug"))
FPKM_MATRIX_FILE          <- Sys.getenv("edgeRNaiveFpkmFile")  # must be output from edgeR! Usually named "fpkm_naive.txt" in the "50.edgeR_diff_expr" directory
SUBREAD_RAW_COUNTS_FILE   <- Sys.getenv("subreadRawCountsFile") #  # must be output from subread! Usually named "subread.featureCounts.counts.txt" in the "20.subread_counts" directory
OUTDIR                    <- Sys.getenv("summaryFigureDir") # Should be the full path

DENDROGRAM_PREFIX        <- "Fig_D1_dendrogram"
DENDROGRAM_FLAT_PREFIX   <- "Fig_D2_dendrogram_with_flat_edge"
GAPMAP_GLOBAL_FIG_PREFIX <- "Fig_G1"
PCA_PREFIX               <- "Fig_PCA"

if (interactive() && Sys.getenv("DISPLAY") != "" && Sys.getenv("RSTUDIO_USER_IDENTITY") == "") {
     # It's interactive, AND the display isn't set, AND ALSO we are not in RSTUDIO!
     print("[ERROR]: This stuff totally fails to run, and pops up an X11 window annoyingly, if 'export DISPLAY='' ' has not been run beforehand to prevent X11 and Cairo from running. Note that you CANNOT use 'Sys.setenv('DISPLAY','')' in the code--- by then, it's too late to properly affect Cairo support, and it still tries to pop up an X11 window.")
     print("[ERROR]: Note also that the behaviors are DIFFERENT when using 'Rscript <script name>' and sourcing the script from within R. Solution: run export DISPLAY='' before you try this script.")
     stopifnot(Sys.getenv("DISPLAY") == "")
}

AGWDEBUG             <- Sys.getenv("AGWDEBUG") # To use the testing code, do the following in the terminal:  export AGWDEBUG=1 ; export DISPLAY='' ; R, then source("/data/work/Code/alexgw/monkey_agw/bin/51b.summary.R", echo=FALSE)
#AGWDEBUG = TRUE; warning("DANGER! Manual setting of this value.")
#FPKM_MATRIX_FILE        = "fpkm.txt"; warning("DANGER! Manual setting of this value.")
#SUBREAD_RAW_COUNTS_FILE = "subread_counts.txt"; warning("DANGER! Manual setting of this value.")

if (!dir.exists(OUTDIR)) { dir.create(OUTDIR, recursive=TRUE) }
print_stdout_and_stderr("Changing directory to the (possibly newly-created) output directory ", OUTDIR, "...")
setwd(OUTDIR) # Let's go to the output directory! Wherever it is.
print("Generating summary figures now...")
everything           <- read.table(FPKM_MATRIX_FILE, sep="\t", header=T, row.names=1, check.names=FALSE)
fpkm_matrix          <- data.matrix(everything)
groupNamesPerSample  <- gsub("[.].*", "", colnames(fpkm_matrix), perl=T) # Remove everything after the first "."
#          Example of what groupNamesPerSample looks like:    c("GeneDel", "GeneDel", "WT", "WT", "DRUGX, "DRUGX")

resetPar <- function(prevPar) {
     #saved.par.default <- par(no.readonly=TRUE)
     par(prevPar)
}

ourPdf <- function(fileprefix, width=12, height=12) {
     fileprefix = sub("[.]pdf$", "", fileprefix, perl=T, ignore.case=T)
     pdf(file=paste0(fileprefix,".pdf"), width=width, height=height, onefile=FALSE) # 'onefile=FALSE' prevents the blank page in PDFs
     par(pty='s')
}
ourPng <- function(fileprefix, width=2400, height=2400, res=144) {
     fileprefix = sub("[.]png$", "", fileprefix, perl=T, ignore.case=T)
     png(file=paste0(fileprefix,".png"), width=width, height=height, res=res)
     par(pty='s')
}

# Just for fixing a problem for samples with super long names in GapMap
agw_single_linebreak <- function(names.vec, len=32) {
     enough_chars <- paste0("^(.{", len, ",}[-._ ])")
     # Note that we add two spaces at the beginning of the new (wrapped) line
     return( sub(enough_chars, "\\1\n  ", names.vec, perl=T) ) # just add a newline and some spaces to indicate line wrapping after this many characters AND some kind of space/hyphen/etc. But only break it ONCE!
}

# Takes the species, quality score, and ".bam" off the sample name, so that it fits on the graph/plot.
agw_clip_bam_paths <- function(names.vec) {
     names.vec <- basename(names.vec)
     names.vec <- sub("[._][sb]am", "", names.vec, perl=T, ignore.case=T) # remove the SAM / BAM suffix
     names.vec <- sub("[-._](galGal[0-9]+|mm[0-9]+|hg[0-9]+)[-._]q[0-9]+", "", names.vec, perl=T, ignore.case=T) # remove the species-and-then-quality suffix (e.g., _mm10_q30, if any). Only do this if it's a recognized species.
     return(names.vec)
}
fpkm_log_matrix <- log2(1+fpkm_matrix)

n_samples <- ncol(fpkm_log_matrix)
# To make a gapped cluster heatmap, you need to pass a matrix object for heatmap, and dendrogram class objects for drawing dendrograms and ordering.
corr_color_scale = gplots::colorpanel(64,"yellow","red","black")  #corr_color_scale = gray.colors(32, start=0.0, end=1.0, gamma=1.0)


print("Clustering for generating dendrograms...")
hc2 <- hclust(stats::dist(t(fpkm_log_matrix), method="minkowski"), "ward.D2") ## transpose it so we cluster ARRAYS and not genes

dendro_height = 16 # inches
dendro_width  = 12 + n_samples/3.0 # give the samples enough horizontal "breathing room" when there are a ton of them

ourPdf(DENDROGRAM_PREFIX, height=dendro_height, width=dendro_width)
par(mar=c(10, 4, 4, 10), pty='s')
plot(hc2, lwd=4, lty=1, col="black", col.lab="red", xlab="Samples", main="Samples clustered by log2(FPKM+1) using Ward's clustering & Euclidean distance")
resetPar(prevPar=saved.par.default)
dev.off()  ##

ourPdf(DENDROGRAM_FLAT_PREFIX, height=dendro_height, width=dendro_width)
par(mar=c(10, 4, 4, 10), pty='s')
plot(hc2, lwd=4, lty=1, hang=-1, xlab="Samples", main="Samples clustered by expression value: log2(FPKM+1)")
resetPar(prevPar=saved.par.default)
dev.off()  ##

agw_do_gapmap <- function(mmm, uniquePrefix="") {
     # ============ GAPMAP of sample similarity ==============
     print("Clustering for a GapMap...")
     colnames(mmm) <- paste0("  ", agw_single_linebreak(colnames(mmm), len=36), "  ") # gapmap is really dumb and doesn't put enough spacing between the labels and the matrix
     colnames(mmm) <- agw_clip_bam_paths(colnames(mmm))

     gapmat <- t(mmm)
     if (!is.null(rownames(gapmat))) { rownames(gapmat) <- paste(" ", rownames(gapmat), " ", sep='') }
     if (!is.null(colnames(gapmat))) { colnames(gapmat) <- paste(" ", colnames(gapmat), " ", sep='') }
     distxy_gapmap <- stats::dist(gapmat, method="minkowski")      #calculate distance matrix. default is Euclidean distance
     hc_gapmap     <- stats::hclust(distxy_gapmap, "ward.D2") #perform hierarchical clustering. default is complete linkage.
     dend_gapmap   <- as.dendrogram(hc_gapmap)
     print("Now actually generating a GapMap...")
     HUGE_PDF_DIM = 30  # has to be HUGE width and height to keep the text from going off the edge
     
     with_gap_filename <- paste0(GAPMAP_GLOBAL_FIG_PREFIX, "_a", uniquePrefix, "_GapMap_Minkowski_Dist")
     no_gap_filename   <- paste0(GAPMAP_GLOBAL_FIG_PREFIX, "_b", uniquePrefix, "_GapMap_Minkowski_Dist_No_Gaps")
     
     ourPdf(with_gap_filename, width=HUGE_PDF_DIM, height=HUGE_PDF_DIM)
     gapmap::gapmap(m=as.matrix(distxy_gapmap), mode="quantitative", d_row=rev(dend_gapmap), d_col=dend_gapmap, col=corr_color_scale, show_legend=TRUE, main="Gapmap of sample similarity", label_size=2.8)
     dev.off()
     ourPdf(no_gap_filename, width=HUGE_PDF_DIM, height=HUGE_PDF_DIM)  # in gapmap below: mode=threshold with no actual threshold means "no gaps"  # has to be HUGE width and height to keep the text from going off the edge
     gapmap::gapmap(m=as.matrix(distxy_gapmap), mode="threshold", d_row=rev(dend_gapmap), d_col=dend_gapmap, col=corr_color_scale, show_legend=TRUE, main="Gapmap of sample similarity", label_size=2.8) # label_size is hard to use properly
     dev.off()
     # ============ GAPMAP of sample similarity ==============
}



agw_do_pca <- function(mmm, groups, uniquePrefix="", min_reads_per_row=25, ntop_most_interesting_genes=1500, additional_main="") {
     # min_reads_per_row:             At least THIS many counts is required in a row for it to be examined at all. Otherwise there is essentially no data
     # ntop_most_interesting_genes:   DESEQ2 by default only uses the top **500** most expressed  genes!
     # additional_main:   Additional 'main' text
     
     # Also makes a couple of heatmaps
     # ============ Make a PCA plot ==============
     PDFWIDTH  <- 14
     PDFHEIGHT <- 14
     
     # Heatmap uses a different set from the PCA plot... but why
     HEATMAP_PREFIX <- "High_variance_genes_row_means_subtracted"
     HEATMAP_HEIGHT <- 45 # inches
     HEATMAP_WIDTH  <- 15
     MAX_NUM_GENES_FOR_TOP_VARIANCE_HEATMAP <- ntop_most_interesting_genes
     
     gapless_heatmap_filename    <- paste0(HEATMAP_PREFIX, "_a", uniquePrefix, "_no_gaps")
     visual_gap_heatmap_filename <- paste0(HEATMAP_PREFIX, "_b", uniquePrefix, "_visual_gaps")
     heatmap_matrix_filename     <- paste0(HEATMAP_PREFIX, "_c", uniquePrefix, "_raw_data_matrix.txt.tab")
     
     barplot_filename                 <- paste0(PCA_PREFIX, "_a", uniquePrefix)
     pca_prefix_filename              <- paste0(PCA_PREFIX, "_b", uniquePrefix)
     pca_gene_subset_matrix_filename  <- paste0(PCA_PREFIX, "_c", uniquePrefix, "_genes_used_in_PCA_plot--matrix.txt")  # only genes that made our cutoff and were used in the PCA plot
     
     conditionFactor  <- factor(groupNamesPerSample)
     cold             <- data.frame("condition"=conditionFactor, "sample_name"=colnames(mmm))
     ddsMat           <- DESeq2::DESeqDataSetFromMatrix(countData=mmm, colData=cold, design=~condition);
     colnames(ddsMat) <- colnames(mmm)
     
     # Note: there is probably some weird issue here still (try monkey_agw --test1 to replicate it)
     # You should see this error message, generated by one of the calls to DESeq:
     #         Error in lfproc(x, y, weights = weights, cens = cens, base = base, geth = geth,  :
     #         newsplit: out of vertex space
     #         Calls: <Anonymous> ... estimateDispersionsFit -> localDispersionFit -> locfit -> lfproc -> .C
     #         In addition: There were 17 warnings (use warnings() to see them)
     # Cause is currently unknown to me. It doesn't seem like a fatal error, though... (?)
     
     # "blind" - Ignore the sample labels and compute a gene's empirical dispersion value as if all samples were replicates of a single condition. This can be done even if there are no biological replicates. This method can lead to loss of power; see the vignette for details. The single estimated dispersion condition is called "blind" and used for all samples.
     # Blind = true will work even when there are NO REPLICATES for ANY EXPERIMENTAL GROUP.
     # Blind = false will work if at least ***ONE*** experimental group has replicates, however.
     # "If many genes have large differences in counts due to the experimental design, it is important to set blind=FALSE for downstream analysis."
     
     any_condition_has_replicates <- (length(levels(ddsMat@colData$condition)) < length(ddsMat@colData$condition)) # At least ONE condition has replicates (i.e. there are MORE rows than unique experimental conditions)
     if (any_condition_has_replicates) { blind_status <- FALSE } else { blind_status <- TRUE }
     rld    <- DESeq2::rlog(ddsMat, blind=blind_status) # LOG TRANSFORMED.
     # to get the values out of 'rld', use the 'assay' function: assay(rld)
     
     print("Handling barplot...")
     # ========= BARPLOT OF PC1 through PC16 ==============
     rv               <- rowVars.copied_from_genefilter(assay(rld)) # <-- This code (through the pca line) is verbatim from the function "DESeq2::plotPCA"
     select           <- order(rv, decreasing = TRUE)[seq_len(min(ntop_most_interesting_genes, length(rv)))]   # <-- select the top genes
     our.subset       <- assay(rld)[select, ]  # jsut a subset...
     write.table(our.subset, file=pca_gene_subset_matrix_filename, sep="\t", row.names=F, col.names=T, quote=F) # Save the list of genes we used for PCA plotting (ones that had enough expression/variance/etc...)
     pca              <- prcomp(t(our.subset))
     allPer.vec       <- 100 * summary(pca)$importance[2,]
     allPer.text      <- paste(format(allPer.vec,digits=0, nsmall=1, scientific=FALSE), "%", sep='')
     allPerLabels.vec <- paste(names(allPer.vec), "\n", allPer.text, sep='')
     
     ourPdf(barplot_filename, width=PDFWIDTH, height=PDFHEIGHT)
     barplot(allPer.vec, main=paste0("PCA: % variance explained by each component\nCalculated with the ", length(select), " highest-variance genes", additional_main), ylim=c(0,100), names.arg=allPerLabels.vec)
     dev.off()
     # ========= BARPLOT OF PC1 through PC16 ==============
     print("Barplot done...")
     
     # ========= PCA PLOT ===================
     pca.modified_version_of_plotPCA_from_DESeq2 = function(object, intgroup="condition", topn=500, returnData=TRUE) { # always return data
          # This version returns ALL components, not just 1 and 2. It's a modified version of the DESeq2 code from: # https://github.com/Bioconductor-mirror/DESeq2/blob/master/R/plots.R
          rv         <- rowVars.copied_from_genefilter(assay(object))     # use highest VARIANCE genes.
          select     <- order(rv, decreasing=TRUE)[seq_len(min(topn, length(rv)))] # select the top N genes by the ordering in 'rv'
          pca        <- prcomp(t(assay(object)[select,])) # perform a PCA on the data in assay(x) for the selected genes
          if (!all(intgroup %in% names(colData(object)))) { stop("Failure in trying to calculate PCA---the argument 'intgroup' should specify columns of colData(ddsMat)") }
          intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
          # add the intgroup factors together to create a new grouping factor
          group <- if (length(intgroup) > 1) { factor(apply( intgroup.df, 1, paste, collapse=" : ")) }
                   else { colData(object)[[intgroup]] }
          d <- data.frame(pca$x, group=group, intgroup.df, name=colnames(object)) # returns PC1, PC2, ... PC8... etc
          percentVar <- pca$sdev^2 / sum( pca$sdev^2 ) # the contribution to the total variance for each component
          names(percentVar) <- colnames(pca$x)
          attr(d, "percentVar")    <- percentVar          # <-- the percent variance for PC1... PCwhatever
          attr(d, "numComponents") <- length(percentVar)
          return(d)
     }
     
     #pcaDat         <- DESeq2::plotPCA(rld, intgroup = c( "condition", "sample_name"), returnData=TRUE, ntop=ntop) # <-- returns data, does NOT plot it
     pcaDat          <- pca.modified_version_of_plotPCA_from_DESeq2(rld, intgroup=c("condition","sample_name"), topn=ntop_most_interesting_genes, returnData=TRUE) # <-- returns data, does NOT plot it.
     print(pcaDat) # Print a summary of the data...
     
     agw.plot_pca_as_pdf <- function(thePcaDat, filename_prefix, pcx_string, pcy_string) {
          lets <- c(LETTERS, letters, paste(LETTERS, sort(rep(c(2:9), times=length(LETTERS))), sep=''))[1:ncol(rld)] # lets = letters. Picks one letter per sample, A...Z, then a..z, then A2..A2... etc
          shps <- unlist(lapply(lets, utf8ToInt)) # shps = shapes for the PCA plot
          ppp <- ggplot(thePcaDat, aes_string(x=pcx_string, y=pcy_string, shape="sample_name", color="condition"))
          ppp <- ppp + geom_point(size=6)
          ppp <- ppp + scale_shape_manual(values=shps)
          percentVar.text.vec <- round(100 * attr(thePcaDat, "percentVar"),1) # one decimal place
          ppp <- ppp + xlab(paste0(pcx_string," (explains ",percentVar.text.vec[pcx_string],"% variance)")) + ylab(paste0(pcy_string," (explains ",percentVar.text.vec[pcy_string],"% variance)"))
          ppp <- ppp + theme_bw()
          ppp <- ppp + coord_fixed() # preserve aspect ratio
          ppp <- ppp + theme(aspect.ratio=1)
          ppp <- ppp + ggtitle(paste0("PCA: % variance explained by each component\nCalculated with the ", length(select), " highest-variance genes", additional_main))
          pdf(paste0(filename_prefix,"---",pcy_string,"_vs_",pcx_string,".pdf"), width=PDFWIDTH, height=PDFHEIGHT)
          print(ppp) # <-- plot it
          dev.off()
     }
     # See the documentation here: http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf
     MAX_PCA_TO_PLOT <- 5 # do not go beyong PCA5. Minimum possible value for this is 2
     pcount <- 1
     for (i in seq(from=1, to=min(MAX_PCA_TO_PLOT, attr(pcaDat,"numComponents")))) {
          for (j in seq(from=1, to=min(MAX_PCA_TO_PLOT, attr(pcaDat,"numComponents")))) {
               if (j <= i) { next; }
               agw.plot_pca_as_pdf(thePcaDat=pcaDat, filename_prefix=paste0(pca_prefix_filename,"_",pcount), pcx_string=paste0("PC",i), pcy_string=paste0("PC",j)) # Also makes a pdf!
               pcount <- pcount+1
          }
     }
     
     
     
     
     
     
     
     # ==================== Below: This part is UNRELEATED to the PCA part. It's only here because we've parsed the relevant data objects already. ===============
     
     # ======== HEATMAP OF TOP GENES by VARIANCE ==========
     submat_for_clustering    <- assay(rld)
     row_sum_above_zero.b.vec <- rowSums(submat_for_clustering)        > 0
     row_has_variance.b.vec   <- apply(submat_for_clustering, 1, var) != 0 # variance for each row
     submat_for_clustering    <- submat_for_clustering[row_sum_above_zero.b.vec & row_has_variance.b.vec, , drop=F]

     topVarGenes <- head(order(rowVars.copied_from_genefilter(submat_for_clustering), decreasing=TRUE), MAX_NUM_GENES_FOR_TOP_VARIANCE_HEATMAP)
     mat         <- submat_for_clustering[topVarGenes, ]
     
     annot.df  <- as.data.frame(colData(rld)$condition)
     rownames(annot.df) <- colnames(mat)
     colnames(annot.df) <- 1:ncol(annot.df)
     
     print(paste0("Writing the heatmap matrix data to <", heatmap_matrix_filename, ">, in case we need it later"))
     write.table(mat, file=heatmap_matrix_filename, quote=F, sep="\t", row.names=TRUE, col.names=NA)
     
     print("Note: if you happen to get the message: 'gpar' element 'fill' must not be length 0   <-- that means the rownames in annotation_col don't match up with the colnames in 'mat'")
     ourPdf(gapless_heatmap_filename, width=HEATMAP_WIDTH, height=HEATMAP_HEIGHT)
     pheatmap::pheatmap(mat, width=HEATMAP_WIDTH, height=HEATMAP_HEIGHT
                      , annotation_col=annot.df
                      , scale="row"  # actually scale each row
                      , clustering_distance_rows="correlation" # or: euclidean
                      , clustering_distance_cols="correlation"
                      , clustering_method="complete"
                      , cluster_rows=TRUE, cluster_cols=TRUE #cutree_rows=3, cutree_cols=3
                      , show_colnames=T, show_rownames=T  #, kmeans_k=6 # aggregates rows into a total of k clusters
                      , display_numbers=FALSE
                      , border_color="#00000000" # note: the last two digits here are the alpha channel/ transparency
                      , main=paste("Top ", nrow(mat), " highest-variance genes, row-scaled. Clust dist is pearson corr", sep='')
                      , fontsize_row=3.5)
     dev.off()
     
     # second heatmap, same as above but with visually apparent GAPS
     ourPdf(visual_gap_heatmap_filename, width=HEATMAP_WIDTH, height=HEATMAP_HEIGHT)
     n_row_clusters <- ceiling(log2(ncol(mat)))+1 # Note: this is based on the number of SAMPLES, not actually the number of genes/rows. This is somewhat arbitrary.
     n_col_clusters <- ceiling(log2(ncol(mat)))
     pheatmap::pheatmap(mat, width=HEATMAP_WIDTH, height=HEATMAP_HEIGHT
                      , annotation_col=annot.df
                      , scale="row"  # actually scale each row
                      , clustering_distance_rows="correlation" # or: euclidean
                      , clustering_distance_cols="correlation"
                      , clustering_method="complete"
                      , cluster_rows=TRUE, cluster_cols=TRUE, cutree_rows=n_row_clusters, cutree_cols=n_col_clusters # by default, divide it up by log2(# items)
                      , show_colnames=T, show_rownames=T
                      , border_color="#00000000" # note: the last two digits here are the alpha channel/ transparency
                      , main=paste("Top ", nrow(mat), " highest-variance genes, row-scaled. Clust dist is pearson corr\nNumber of row clusters specified: ", n_row_clusters, sep='')
                      , fontsize_row=3.5)
     dev.off()
}

agw_do_gapmap(fpkm_matrix    , uniquePrefix="_linear_FPKM");
agw_do_gapmap(fpkm_log_matrix, uniquePrefix="_log2_FPKM");

# ======================= READ A SUBREAD FEATURECOUNTS 'COUNTS' MATRIX =================
# Read the SUBREAD RAW COUNTS file into the matrix "all.counts.mat"
all.ddd  <- read.table(gzfile(SUBREAD_RAW_COUNTS_FILE), sep="\t", header=T, check.names=F, stringsAsFactors=FALSE, fill=FALSE, row.names=1, quote='', comment.char='#') # <-- comment.char is important: the first line of a subread file is a '#' comment (can also use nskip = 1)
colnames(all.ddd) <- agw_clip_bam_paths(colnames(all.ddd))  #all.dt  <- data.table::data.table(ids=rownames(all.ddd), original.data=all.ddd) # all.dt includes the ANNOTATION ROWS
firstDataColIdx   <- (1+which(toupper(colnames(all.ddd))=="LENGTH")); print("Assuming that the last column of annotation is named 'Length' (case-insensitive)...")
all.counts.mat    <- data.matrix(all.ddd[, firstDataColIdx:ncol(all.ddd), drop=F]) # all.counts.mat = only DATA
stopifnot(all(!is.na(all.counts.mat))) # There should not be any NA values in the RAW COUNTS matrix!
# ======================= DONE READING A SUBREAD FEATURECOUNTS 'COUNTS' MATRIX =================


MIN_READS_FOR_EXPR = max(10, 2*ncol(all.counts.mat)) # require at least 10 reads total, OR an average of 2 reads per sample (whichever is greater)
expressed.genes.mat  <- all.counts.mat[ base::rowSums(all.counts.mat) >= MIN_READS_FOR_EXPR, ]      # Subset of only "expressed" genes

agw_do_pca(mmm=expressed.genes.mat, uniquePrefix=paste0("--A-top--500--genes_with_count_above_", MIN_READS_FOR_EXPR), ntop_most_interesting_genes=500 , additional_main=paste0(" with at least ", MIN_READS_FOR_EXPR, " reads."))
agw_do_pca(mmm=expressed.genes.mat, uniquePrefix=paste0("--B-top-1000--genes_with_count_above_", MIN_READS_FOR_EXPR), ntop_most_interesting_genes=1000, additional_main=paste0(" with at least ", MIN_READS_FOR_EXPR, " reads."))
agw_do_pca(mmm=expressed.genes.mat, uniquePrefix=paste0("--C-top-1500--genes_with_count_above_", MIN_READS_FOR_EXPR), ntop_most_interesting_genes=1500, additional_main=paste0(" with at least ", MIN_READS_FOR_EXPR, " reads."))
agw_do_pca(mmm=expressed.genes.mat, uniquePrefix=paste0("--D-all_genes_with_count_above_", MIN_READS_FOR_EXPR), ntop_most_interesting_genes=999999, additional_main=paste0(" with at least ", MIN_READS_FOR_EXPR, " reads."))
######### Do not actually use this!  agw_do_pca(mmm=all.counts.mat     , uniquePrefix=paste0("--top-1500--all-genes")                                 , ntop_most_interesting_genes=1500, additional_main=paste0(" (with no minimum read cutoff)."))

if (FALSE) {
     # Showing the individual genes that contributed to each PCA.
     eee <- rbind("SNAKES"=seq(from=0, to=15000, by=1000)
                 ,"3CAKES"=seq(from=0, to=15000, by=1000)
                 ,"4CAKES"=seq(from=0, to=15000, by=1000)
                 ,"5CAKES"=seq(from=0, to=15000, by=1000)
                 ,"6CAKES"=seq(from=0, to=15000, by=1000)
                 ,"7CAKES"=seq(from=0, to=15000, by=1000)
                 ,"3xCAKES"=seq(from=0, to=45000, by=3000)
                 ,"2xCAKES"=seq(from=0, to=30000, by=2000)
                 ,"NEGATIVE"=seq(from=15000, to=0, by=-1000)
                 ,"NEGATIVE3"=seq(from=0, to=-45000, by=-3000)
                , head(all.counts.mat, n=80))
     
     pca.object <- prcomp(eee, center=TRUE, scale=TRUE, retx=TRUE)  # PCA with centering and scaling
     pca.object$rotation  # The loadings are here
     
     pca.object$scores
     z <- stats::predict(pca.object) # <------- this is showing the individual genes that contributed to each PCA. 
     
     pheatmap(log2(1+eee[ names(head(sort(abs(z[,1]), decreasing=TRUE),n=30)), ]))
     #prcomp(USArrests, scale = TRUE)
     #z <- prcomp(~ Murder + Assault + Rape, data = USArrests, scale = TRUE, retx=T)
     #plot(prcomp(USArrests))
     #summary(prcomp(USArrests, scale = TRUE))
     #biplot(prcomp(USArrests, scale = TRUE))
}
