#!/wynton/group/gladstone/third_party/Rscript  --no-save  --no-restore  --no-site-file  --no-init-file

# Warning: single quotes are NOT ALLOWED on any of the PBS lines above! They trigger a "qsub: unmatched" error and exit code 1.

print0   <- function(...) { print(paste0(...)) } # Print + paste0
print_stdout_and_stderr <- function(...) { write(paste(..., sep=''), stderr()); write(paste(..., sep=''), stdout()); }

options(error=traceback) # useful even non-interactively
#options(error=recover)
options(stringsAsFactors=FALSE)

if (!require(gplots))    { print("You need to install 'gplots'--install.packages('gplots')"); stopifnot(FALSE == "gplots was not installed") }
print("Beware of QUOTES if you are loading excel data! Then you need quote='\"' ")

# See below for filenames

print("Note: this script is unusual, because you can run it both via monkey OR ALSO as a command-line program!")
print("How to run this program on the commmand line (way #1):   Rscript /work/monkey_development/bin/51c.cluster.R YOUR_FPKM_FILE.fpkm.txt")
print("How to run this program via monkey (way #2):             Set a bunch of environment variables (mostly edgeRNaiveFpkmFile and summaryFigureDir) and just call it directly with no arguments.")
print("The output will be named based on your input file... maybe.")

looks_true <- function(string) { return(!is.null(string) && (length(string) >= 1) && !(toupper(string) %in% c("", "F","FALSE","UNDEF","0"))) }
VERBOSE               <- looks_true(Sys.getenv("verbose"))
FORCE                 <- looks_true(Sys.getenv("force"))
FPKM_MATRIX_FILE      <- Sys.getenv("edgeRNaiveFpkmFile")
OUTDIR                <- Sys.getenv("summaryFigureDir")
N_GROUPS              <- Sys.getenv("nGroups")
##DEBUG                 <- looks_true(Sys.getenv("debug"))

if (nchar(FPKM_MATRIX_FILE) > 0 || nchar(OUTDIR) > 0) {
     IS_CMD_LINE_INVOCATION = FALSE # We are getting it from environment variables instead
     print("Since FPKM_MATRIX_FILE and/or OUTDIR were defined via environment variables, it looks like we are running this script with the inputs as environment variables, instead of as passed-in command line arguments.")
     stopifnot(nchar(FPKM_MATRIX_FILE) > 0 && nchar(OUTDIR) > 0)
     if (!dir.exists(OUTDIR)) { dir.create(OUTDIR, recursive=TRUE) }
     print_stdout_and_stderr("Changing directory to the (possibly newly-created) output directory ", OUTDIR, "...")
     setwd(OUTDIR) # Let's go to the output directory! Wherver it is.
} else {
     IS_CMD_LINE_INVOCATION = TRUE # Are we getting the filename from the command line ("Rscript script.R FILENAME") or environemnt variables? ('export FPKM_MATRIX_FILE')
     scriptExecuter <- commandArgs()[1] # First item
     if (scriptExecuter == "RStudio" || looks_true(Sys.getenv("AGWDEBUG"))) {
          # To use the testing code, do the following in the terminal:  export AGWDEBUG=1 ; R, then source("/data/work/Code/alexgw/monkey_agw/bin/51b.summary.R", echo=FALSE)
          # Or just run the code from RStudio!
          warning("Note: we are assuming that since you are running this from RStudio, you want INTERACTIVE DEBUGGING.")
          AGW_DEBUG <- TRUE
     } else {
          AGW_DEBUG <- FALSE 
     }
     if (AGW_DEBUG) {
          FPKM_MATRIX_FILE = "/wynton/group/gladstone/biocore/MonkeyPipeline/test_suite/test_fpkm_naive.txt" } # just a test file from another project
     else {
          inputArgs = commandArgs(trailingOnly=TRUE); stopifnot(length(inputArgs) > 0 )
          FPKM_MATRIX_FILE <- inputArgs[1]
     }
}

stopifnot(file.exists(FPKM_MATRIX_FILE))
basefile = basename(FPKM_MATRIX_FILE)
print0("We got the following input FPKM file: <", FPKM_MATRIX_FILE, ">")

# =========================
# FUNCTIONS FOR COMPUTING CV on a vector (cv) or 2D-matrix/data frame (cv.2d)
cv    <- function(x) ( sd(x)/mean(x) )             # calculate coefficient of variation
cv.2d <- function(x) ( apply(x,1,sd)/rowMeans(x) ) # calc coeff of variation on a matrix or data frame
# =========================

# Function for re-ordering a matrix or data frame in the same order that it clusters
dfInClusterOrder <- function(x, clusterObj, clusterAssignments=NULL, includeClusterAssignmentColumn=FALSE) {
  stopifnot(class(clusterObj) == "hclust") # Expects an hclust object
  stopifnot(is.null(clusterAssignments) || (is.integer(clusterAssignments) && length(clusterAssignments) == nrow(x)))
  stopifnot(is.matrix(x) || is.data.frame(x)) # requires a 2d data frame or matrix currently
  # Re-order a matrix or data frame in the SAME ROW ORDER AS IN THE input DENDROGRAM!!!
  # Useful for writing data to an output file in the same order that it appears in a heatmap!
  row_idx = rev(stats::order.dendrogram(as.dendrogram(clusterObj))) # Re-order in dendrogram order
  col_idx = seq(from=1, to=ncol(x)) # Do NOT reorder the columns!
  if (includeClusterAssignmentColumn) {
    return(data.frame("CLUSTER_ASSIGNMENT"=clusterAssignments[ rownames(x)[row_idx] ],   x[row_idx, col_idx])) # Always returns a data frame if cluster assignments are included
  } else {
    return(x[row_idx, col_idx]) # Could be a matrix OR data frame, depending on the input
  }
}


print("Assuming that the input file is: FPKM values, with one ID column on the left and one column header.")
# File should look something like the lines below:
#ids	            ctrl.1 exp.1	other.1
#ENSG00000223972	974.4	 849.4	454.0
#ENSG00000227232	7.6  	 7.9    4.0
#ENSG00000243485	2.7  	 2.8  	1.9

x.all.df <- read.table(FPKM_MATRIX_FILE, sep="\t", header=T, check.names=F, stringsAsFactors=F, row.names=NULL, quote='', comment.char='')
rownames(x.all.df) <- make.names(x.all.df[,1, drop=T], unique=TRUE)
x.mat <- data.matrix(x.all.df[, 2:ncol(x.all.df)])

# THRESHOLDS FOR REMOVING ROWS
MIN_ALLOWED_AVG <- 1.0 # Filter out anything with very low expression
MAX_ALLOWED_AVG <- 5000 # Filter out a few mitochondrial genes

too_hi <- rowMeans(x.mat) > MAX_ALLOWED_AVG
too_lo  <- rowMeans(x.mat) < MIN_ALLOWED_AVG

print0(sum(too_hi), " genes/features were TOO HIGH (average > ", MAX_ALLOWED_AVG, ")")
print0(sum(too_lo), " genes/features were TOO LOW  (average < ", MIN_ALLOWED_AVG, ")")

print("Getting the values from x.mat that are in the 'gold' range and not too high or too low.")
x.gold <- x.mat[ !too_hi & !too_lo, , drop=F]

stopifnot(nrow(x.gold) > 2)


print("Calculating coefficient of variation...")
cv.vec <- cv.2d(x.gold)

NUM_INTERESTING_TO_PLOT <- 1000

ok.vec <- is.finite(cv.vec)  &  (order(cv.vec, na.last=T, decreasing=FALSE) <= NUM_INTERESTING_TO_PLOT)
stopifnot(sum(ok.vec) <= NUM_INTERESTING_TO_PLOT) # Make sure we winnowed it down to the top N

x.everything.log <- log2(1.0 + x.mat)

x.top <- x.gold[ok.vec, , drop=F]
x.top.log <- log2( 1.0 + x.top ) # Note that we add a PSEUDOCOUNT of 1.0 to the FPKM here!

# The final value to actually plot
xp <- x.top.log

#ac  <- gplots::colorpanel(26, low="#00FFFF","#000000","#FFFF00") # cyan-black-yellow
ac  <- gplots::colorpanel(26, low="#0099FF","#000000","#FF3300") # blue-black-red
brr <- seq(from=0.0, to=4.0, length.out=(length(ac)+1)) # note: log2 scale!
    
dist.pear  <- function(zzz) as.dist(1-cor(t(zzz)))
hclust.ave <- function(zzz) hclust(zzz, method="average")
#hclustfun=function(x){hclust(x,method='ward.D2')},
#distfun=function(x){dist(x,method='minkowski')}

DIST_METHOD <- "minkowski"
AGG_METHOD  <- "ward.D2"

MAX_ALLOWED_CLUSTER_GROUPS <- 20

colClustEverything <- hclust(dist(t(x.everything.log), DIST_METHOD), AGG_METHOD) ## transpose to cluster SAMPLES, not genes
colClust           <- hclust(dist(t(xp)              , DIST_METHOD), AGG_METHOD) ## transpose to cluster SAMPLES, not genes
rowClust           <- hclust(dist(  xp               , DIST_METHOD), AGG_METHOD) # cluster GENES (slow, due to the large number of them)

# ======= ***************************
# ======= ***************************
# ======= ***************************
nClustGroups <- min(N_GROUPS, MAX_ALLOWED_CLUSTER_GROUPS) # cut the tree at the number of experimental groups, but no more than 20!

#if (nClustGroups > 

COL_DENDRO_FILENAME             <- paste0("Out.", basefile, ".Heatmap_dendrogram_by_sample_highest_CV_genes", ".pdf")
COL_DENDRO_EVERYTHING_FILENAME  <- paste0("Out.", basefile, ".Heatmap_dendrogram_by_sample_all_genes", ".pdf")
ROW_DATA_FILENAME               <- paste0("Out.", basefile, ".Heatmap_data_in_dendrogram_order--", nClustGroups, "_clusters.txt")
FIG_FILENAME_PREFIX             <- paste0("Out.", basefile, ".Heatmap_top_", nrow(xp), "_genes_by_CV")

# ======= ***************************
# ======= ***************************
# ======= ***************************

row.cluster.vec <- cutree(rowClust, k=nClustGroups) # Actually generate the clusters by cutting the tree
stopifnot(nrow(xp) == length(row.cluster.vec))

# Write the output in the SAME ROW ORDER AS IN THE DENDROGRAM!!!
data_ordered_plus_cluster <- dfInClusterOrder(xp, rowClust, row.cluster.vec, includeClusterAssignmentColumn=TRUE)
#xpclust <- dfInClusterOrder(xp, rowClust)

write.table(data_ordered_plus_cluster, ROW_DATA_FILENAME, quote=F, sep="\t",row.names=T, col.names=T)
print0("Generated the output file <", ROW_DATA_FILENAME, ">...")

clusColors.vec <- rainbow(n=max(row.cluster.vec), s=1.0, v=0.8)[row.cluster.vec]

for (ttt in c("png", "pdf")) {
  if (ttt == "png") {
    out <- paste0(FIG_FILENAME_PREFIX,".",ttt)
    png(out, width=3000, height=3000, res=144)
  } else if (ttt == "pdf") {
    out <- paste0(FIG_FILENAME_PREFIX,".",ttt)
    pdf(out, width=20, height=20)
  }
  heatmap.2(xp, trace="none", Colv=NA, Rowv=as.dendrogram(rowClust), dendrogram="row"
            , col=ac, breaks=brr, margins=c(15,25), labRow="", cexRow=0.5
            , scale="none"
            , main=paste0("The ", nrow(xp), " highest-CV genes with mean(FPKM) between ", MIN_ALLOWED_AVG, " and ", MAX_ALLOWED_AVG, ".\nPlotted values are log2(1.0 + FPKM).")
            , denscol="#FF00FF"
            , RowSideColors=clusColors.vec
            , keysize=0.9)
  dev.off()
  print0("Generated the output file <", out, ">...")
}


pdf(COL_DENDRO_EVERYTHING_FILENAME, width=20, height=20)
# should probably "par(mar=c(A,B,C,D))" to make sure that the labels stay on the figure even if they're really long!
plot(as.dendrogram(colClustEverything), main=paste0("Dendrogram by log2(FPMK) of all ", nrow(x.everything.log), " genes, by sample\n", DIST_METHOD, " distance / ", AGG_METHOD, " agglomeration"))
dev.off()
print0("Generated the output dendrogram figure <", COL_DENDRO_EVERYTHING_FILENAME, ">...")


pdf(COL_DENDRO_FILENAME, width=20, height=20)
# should probably "par(mar=c(A,B,C,D))" to make sure that the labels stay on the figure even if they're really long!
plot(as.dendrogram(colClust), main=paste0("Dendrogram by log2(FPMK) of top ", nrow(xp), " highest-CV acceptable-FPKM genes\n", DIST_METHOD, " distance / ", AGG_METHOD, " agglomeration"))
dev.off()
print0("Generated the output dendrogram figure <", COL_DENDRO_FILENAME, ">...")


#plot(as.dendrogram(rowClust))

# ==================================================================================

# Sean's RAPID CLUSTERING (not yet integrated)
#
#Rapid Clustering!
#
#The output object v contains two variables
#
#v$od  contains the order of items in the heatmap
#v$allClust contains the meta cluster ID for each item
#
#so
#
#v$allClust[v$od]  would print  1 1 1 1 1 1 1 1 2 2 2 2 2  etc.
if (1 == 2) {
     ### start with a large matrix of data 'pks' that has up to a million rows and a couple hundred columns at most
     pks <- "???????" # some matrix of data
     # define colors
     bc <- c(0.1000,0.1000,0.3)
     wc <- c(1,1,1)
     yc <- c(1,1,0)
     xc <- 1:100/100
     cvWB <- rgb((bc[1]-wc[1])*xc+wc[1], (bc[2]-wc[2])*xc+wc[2], (bc[3]-wc[3])*xc+wc[3])

     # rapid hopach-based clustering
     getClusters <- function(m,nc) {
          k <- kmeans(m,nc,nstart=10,iter.max=100)
          require(hopach)
          disttype  <- "spearman"
          distfile  <- "distfile.Robject"
          clustfile <- "hopach.Robject"
          mat <- cor(t(k$centers), method = disttype)
          mat <- 1 - mat
          mat[is.na(mat)] <- -99999
          gene.dist <- as.hdist(mat)
          gene.dist@Call <- disttype
          gene.hobj = hopach(k$centers, dmat=gene.dist, clusters="best", verbose=TRUE)
          kclustn <- 1:length(gene.hobj$final$order)*0
          for (i in 1:length(gene.hobj$clustering$medoids)) {
               kclustn[which(gene.hobj$clustering$labels == gene.hobj$clustering$labels[gene.hobj$clustering$medoids[i]])] <- i
          }
          cod <- order(kclustn)
          rbc <- rainbow(max(kclustn))
          allClustN <- kclustn[k$cluster]
          od <- order(allClustN)
          tlist <- list("rbc" = rbc, "allClustN" = allClustN, "od" = od)
          return(tlist);
     }

     # set nClust to a manageable number of kmeans clusters that hopach can handle
     # don't set nClust too high - the function will fail if nClust is larger than the number of data points
     nClust <- 50

     # do the clustering and make a heatmap
     v <- getClusters(pks,nClust)
     heatmap(pks[v$od,], scale="none", col=cvWB, Rowv=NA, Colv=NA, RowSideColors=v$rbc[v$allClustN[v$od]], labRow=NA, labCol=NA)
}

