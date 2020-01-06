<<<<<<< HEAD
#!/bin/bash
/bin/env Rscript /wynton/group/gladstone/biocore/MonkeyPipeline/bin/51a.monocle.bak.R
=======
#!/wynton/group/gladstone/third_party/Rscript  --no-save  --no-restore  --no-site-file  --no-init-file




# Warning: single quotes are NOT ALLOWED on any of the PBS lines above! They trigger a "qsub: unmatched" error and exit code 1.

# Note: you can DEBUG this script by setting "export AGWDEBUG=1; R" and then loading the script from within R. It will use FAKE data instead of real data! Make sure to set AGWDEBUG=0 afterward though, unless you want to keep running the debug code!!

require(data.table)
require(gplots)
require(monocle) # http://www.bioconductor.org/packages/release/bioc/vignettes/monocle/inst/doc/monocle-vignette.pdf
require(reshape) # reshape for 'melt'
#source("/wynton/home/pollard/wmaguire/gladstone/biocore/references/agwUtil.R"); # Has "pairs.agw" in it
options(error=recover)
options(stringsAsFactors=FALSE)

# Note: if you get an X11 error when debugging with R and loading this file using "source", make sure to not just say source, but include echo=FALSE!
# Example: source("/data/work/Code/alexgw/monkey_agw/bin/51a.monocle.R", echo=FALSE)
# See this thread for details: http://stackoverflow.com/questions/6675066/ggplots-qplot-does-not-execute-on-sourcing

print(paste("The 'interactive' status for the R script being run this way is:", interactive()))

if (interactive() && Sys.getenv("DISPLAY") != "") {
     print("[ERROR]: This stuff totally fails to run, and pops up an X11 window annoyingly, if 'export DISPLAY='' ' has not been run beforehand to prevent X11 and Cairo from running. Note that you CANNOT use 'Sys.setenv('DISPLAY','')' in the code--- by then, it's too late to properly affect Cairo support, and it still tries to pop up an X11 window.")
     print("[ERROR]: Note also that the behaviors are DIFFERENT when using 'Rscript <script name>' and sourcing the script from within R. Solution: run export DISPLAY='' before you try this script.")
     stopifnot(Sys.getenv("DISPLAY") == "")
}

looks_true <- function(string) { return(!is.null(string) && !(toupper(string) %in% c("", "F","FALSE","UNDEF","0"))) }
print_stdout_and_stderr <- function(...) { write(paste(..., sep=''), stderr()); write(paste(..., sep=''), stdout()); }

VERBOSE                   <- looks_true(Sys.getenv("verbose"))
FORCE                     <- looks_true(Sys.getenv("force"))
DEBUG                     <- looks_true(Sys.getenv("debug"))

FPKM_MATRIX_FILE          <- Sys.getenv("edgeRNaiveFpkmFile")
OUTDIR                    <- Sys.getenv("monocleDir") # Should be the full path
MONOCLE_FINAL_OUTPUT_FILE <- file.path(OUTDIR, "z.done.monocle")

MIN_EXPR_THRESH                <- 0.50 # in fpkm
QTHRESH                        <- 0.10 # q-value (FDR) threshold: really should be like 0.1 or something
MIN_SAMP_FRACT                 <- 0.50 # Among all the samples, this gene must be expressed in at least X fraction of them (to start--if the result is 0, we will gradually reduce this)

SAMPLE_INFO_SHEET   <- "monocle_sample_sheet.auto_generated.tmp.txt" # We generate these files from the input fpkm file...
GENE_ANNOT          <- "monocle_gene_info.auto_generated.tmp.txt"    # We generate these files from the input fpkm file...

ourPdf <- function(filename, width=12, height=12) {
     par(pty='s') # square plotting region
     pdf(filename, width=width, height=height)
}

AGWDEBUG            <- Sys.getenv("AGWDEBUG") # To use the testing code, do the following in the terminal:  export AGWDEBUG=1 ; export DISPLAY='' ; R, then source("/data/work/Code/alexgw/monkey_agw/bin/51a.monocle.R", echo=FALSE)
if (nchar(AGWDEBUG) > 0 && as.numeric(AGWDEBUG)>0 && toupper(AGWDEBUG) != "FALSE" && toupper(AGWDEBUG) != "F") {
     # Make fake "FPKM" data for testing purposes!
     if (grepl("[/]monkey", getwd())) { print(">>>>>>>>> Hey you! Don't run this test code from the MONKEY directory, because it will pollute it with a bunch of test files! Make a test directory (that is ***not*** a subdirectory of the monkey directory, and doesn't start with 'monkey') and run this test code in there instead!!!!"); stopifnot(1==2); }
     FPKM_MATRIX_FILE = "tmp.test.fpkm.fake.txt"
     OUTDIR <- "./"
     NR   <- 250 # num fake genes
     FAKE_NAMES = sample( paste("ENSFAKE", sprintf("%06d", seq(NR)), sep='') ) # Note that we RANDOMLY re-order it so they aren't in alphabetical order!
     FAKE_ADJ = c(0, 99, runif(NR-2)*0.90 + 0.10)  # randomly adjust it. Note the 0 (low) and the 99999 (too high) expression items
     A_ADJ = 1.0 * FAKE_ADJ
     B_ADJ = 1.5 * FAKE_ADJ
     C_ADJ = 1.2 * FAKE_ADJ
     exp_scale_rate <- 0.01
     FAKE_EXPR_BASE = rexp(NR, rate=exp_scale_rate) # Fake data is exponentially distributed. Note: the numbers will be totally unrealistic!
     fake_rescale <- function(NUM) {  return (runif(NUM)*0.20+0.40)  }
     x <- data.table("GENE_ID"=FAKE_NAMES
                     , "aa.1"=fake_rescale(NR) * FAKE_EXPR_BASE * A_ADJ
                     , "aa.2"=fake_rescale(NR) * FAKE_EXPR_BASE * A_ADJ
                     , "aa.3"=fake_rescale(NR) * FAKE_EXPR_BASE * A_ADJ
                     , "bb.1"=fake_rescale(NR) * FAKE_EXPR_BASE * B_ADJ
                     , "bb.2"=fake_rescale(NR) * FAKE_EXPR_BASE * B_ADJ
                     , "bb.3"=fake_rescale(NR) * FAKE_EXPR_BASE * B_ADJ
                     , "cc.1"=fake_rescale(NR) * FAKE_EXPR_BASE * C_ADJ
                     , "cc.2"=fake_rescale(NR) * FAKE_EXPR_BASE * C_ADJ   )
     write.table(x, file=FPKM_MATRIX_FILE, quote=F, sep="\t", row.names=FALSE, col.names=TRUE)
     print("[NOTE]: WE ARE RUNNING THE FAKE TESTING CODE, NOT THE REAL CODE, SINCE 'AGWDBEUG' WAS SET AS AN ENVIRONMENT VARIABLE.")
}

stopifnot(nchar(OUTDIR) > 0)
stopifnot(nchar(FPKM_MATRIX_FILE) > 0)

if (!FORCE && file.exists(MONOCLE_FINAL_OUTPUT_FILE)) {
     print("[OK] Not re-running monocle, because the output appears to already exist.")
} else {
     print("Running Monocle now...")     
     dir.create(OUTDIR, recursive=TRUE)
     print_stdout_and_stderr("Changing directory to the (possibly newly-created) output directory ", OUTDIR, "...")
     setwd(OUTDIR) # Let's go to the output directory! Wherver it is.
     print("Loading the FPKM matrix file now...")
     everything <- read.table(FPKM_MATRIX_FILE, sep="\t", header=T, row.names=1, check.names=FALSE)
     HSMM_expr_matrix     <- data.matrix(everything)
     groupNamesPerSample  <- gsub("[.].*", "", colnames(HSMM_expr_matrix), perl=T) # Remove everything after the first "."
     sampleAnnot.df       <- data.frame("SAMPLE"=colnames(HSMM_expr_matrix), "SAMPLEGROUP_CHAR"=groupNamesPerSample, "SAMPLEGROUP"=as.integer(factor(groupNamesPerSample)), "Pseudotime"=groupNamesPerSample) # Do NOT change 'SAMPLEGROUP' here without also changing it in the expression~SAMPLEGROUP code

     # ========================== GENERATING TWO TEMPORARY FILES IN THE MONOCLE DIRECTORY ========================
     print("Generating some temporary files that we will need...")
     write.table(sampleAnnot.df, SAMPLE_INFO_SHEET, sep="\t", row.names=F, quote=F)
     HSMM_sample_sheet    <- read.delim(SAMPLE_INFO_SHEET, sep="\t", header=T, row.names=1) # <------- ============ WRITING A TEMP FILE
     stopifnot(all(rownames(HSMM_sample_sheet) == colnames(HSMM_expr_matrix)))
     geneAnnot.df <- data.frame("GENE"=rownames(HSMM_expr_matrix), "INFO"="No_info", "gene_short_name"=rownames(HSMM_expr_matrix))
     write.table(geneAnnot.df, GENE_ANNOT, sep="\t", row.names=F, quote=F) # <------- ============ WRITING A TEMP FILE
     HSMM_gene_annotation <- read.delim(GENE_ANNOT, sep="\t", header=T, row.names=1)
     stopifnot(all(rownames(HSMM_expr_matrix) == rownames(HSMM_gene_annotation)))
     # ========================== GENERATING TWO TEMPORARY FILES IN THE MONOCLE DIRECTORY ========================
     
     pd <- new("AnnotatedDataFrame", data=HSMM_sample_sheet)
     fd <- new("AnnotatedDataFrame", data=HSMM_gene_annotation)
     HSMM <- newCellDataSet(HSMM_expr_matrix, phenoData=pd, featureData=fd)
     HSMM <- detectGenes(HSMM, min_expr=MIN_EXPR_THRESH) # Look for EXPRESSED genes that meet the MIN_EXPR_THRESH cutoff!
     print(dim(HSMM)) # Print size of sample
     
     print("Expression looks like this:")
     print(head(exprs(HSMM)))
     print("----------------------------")

     print(head(fData(HSMM)))

     oft_exp_gn <- c() # often-expressed gene names
     MIN_SAMP <- NA
     while (length(oft_exp_gn) == 0 && MIN_SAMP_FRACT >= 0.00) {
          MIN_SAMP <- MIN_SAMP_FRACT * ncol(HSMM) # Minimum number of samples that this gene must be expressed in
          oft_exp_gn    <- row.names(subset(fData(HSMM), num_cells_expressed >= MIN_SAMP)) # "The vector oft_exp_gn now holds the identifiers for genes expressed in at least N cells of the data set."
          if (length(oft_exp_gn) == 0) { MIN_SAMP_FRACT = MIN_SAMP_FRACT - 0.05 }  # Increase the q-value threshold if we seriously get NO genes
          print("Adjusting MIN_SAMP_FRACT...")
     }
     stopifnot(length(oft_exp_gn) >= 1)

     # How to get info about one of these objects:
     print(head(pData(HSMM))) # let's get some other info...
     print(head(fData(HSMM))) # let's get some other info...

     # ======== ERROR CHECKING HERE =======================================
     min_in_matrix = min(exprs(HSMM))
     if (min_in_matrix < 0) {
          print("WARNING: What in the world! The minimum value in the input matrix, which is SUPPOSED to be linear-scaled FPKM data, is less than 0. This should be impossible.")
          system("touch ./WARNING_PROBABLY_MONOCLE_MESSED_UP--INPUT_DATA_HAD_A_VALUE_LESS_THAN_ZERO--IT_WAS_SUPPOSED_TO_BE_FPKM_DATA")
     }
     stopifnot(min_in_matrix >= 0)
     max_minus_min = max(exprs(HSMM)) - min_in_matrix
     if (max_minus_min < 100) { # We really expect our NON-LOG input data to have a linear difference of at least 100. If not, then probably we have log data, which will mess up the log transform below
          print("WARNING: Something is probably SERIOUSLY WRONG with your input data! The (max - min) is not at least 100, which means you probably ALREADY have log-transformed data. But here we're going to log-transform it again? This is PROBABLY an error.")
          system("touch ./WARNING_PROBABLY_MONOCLE_MESSED_UP--INPUT_DATA_DOES_NOT_APPEAR_TO_BE_LINEAR_SCALE_FPKM_DATA")
     }
     # ======== DONE WITH SOME BASIC ERROR CHECKING ========================
     
     print("Log-transforming (1 + expression value) to avoid problems.")
     print("Note: this assumes that the input data is FPKM.")
     L <- log2( 1.0 + exprs(HSMM[oft_exp_gn,])) # Log2-transform each value in the expression matrix, PLUS ONE!

     # ================================================================
     # Standardize each gene, so that they are all on the same scale,
     # Then melt the data with plyr so we can plot it easily"
     print("Melting...")
     melted_dens_df <- reshape::melt(t(scale(t(L))))
     # ================================================================
     
     # Note: GGplot has extremely specific plotting requirements: always use 'ggsave', not pdf(...)
     print("Generating a GGPlot 'qplot' of FPKM density. Note that we added 1 to all the values before log-transforming them...")
     ourGG <- ggplot2::qplot(value, geom="density", data=melted_dens_df) + stat_function(fun = dnorm, size=0.5, color='red') + xlab("Standardized log(FPKM)") + ylab("Density") # Plot the distribution
     ggplot2::ggsave("tmp.fig.monocle-plot-1-FPKM-density.pdf", plot=ourGG, width=12, height=12) # save the plot...
     # ================================================================

     marker_genes <- head(oft_exp_gn, n=25) #head(oft_exp_gn, n=25)
     stopifnot(length(marker_genes) >= 1)
     print(paste("Got a total of this many marker genes (i.e., genes that were expressed frequently enough in our cells to be 'interesting'): ", length(marker_genes), sep=''))
     
     print("Below: SAMPLEGROUP is the column header from the annotation file. Do not change this name without also changing the other occurrences of 'SAMPLEGROUP'. The word 'expression' is literal for whatever reason. Don't change it either!")

     print("Running monocle::differentialGeneTest now...")
     diff_test_res         <- monocle::differentialGeneTest(HSMM[marker_genes,], fullModelFormulaStr="expression~SAMPLEGROUP") # previous fomula str: expression~Media
     print("Done with monocle::differentialGeneTest...")
     ordered_diff_test_res <- diff_test_res[order(diff_test_res$pval),] # order by pval

     sig_gene_names <- c()
     while (length(sig_gene_names) <= 2 && QTHRESH <= 1.0) {
          sig_gene_names <- row.names(subset(diff_test_res, qval < QTHRESH))
          if (length(sig_gene_names) <= 2) { QTHRESH = QTHRESH + 0.05 }  # Increase the q-value threshold if we seriously get NO genes
          print("Adjusting Q-thresh...")
     }
     stopifnot(length(sig_gene_names) >= 2) # Must be at least TWO for this to be meaningful!
     sig_genes <- data.frame(subset(diff_test_res, qval < QTHRESH))
     sig_genes <- merge(fData(HSMM), sig_genes, by="row.names")

     print(paste("The QTHRESH we settled on was: ", QTHRESH, sep=''))
     print(paste("The MIN_SAMP_FRACT we settled on was: ", MIN_SAMP_FRACT, sep=''))
     print(paste("  * Which is equivalent to MIN_SAMP of: ", MIN_SAMP, sep=''))
     
     WANT_TO_PLOT_SOME_SPECIFIC_GENES = FALSE # CURRENTLY WE DO NOT HANDLE THIS IN ANY WAY
     if (WANT_TO_PLOT_SOME_SPECIFIC_GENES) {
          cool_gene_names <- head(rownames(HSMM), n=10)      #cool_gene_names <- HSMM[row.names(subset(fData(HSMM), gene_short_name %in% c("ENSMUSG00000000001", "ENSMUSG00000000088"))),]
          cool_genes      <- HSMM[row.names(subset(fData(HSMM), rownames(HSMM) %in% cool_gene_names)),]

          thresholdStr = paste("MinQVALUE=", QTHRESH, "__MinEXPR_FRACT=", MIN_SAMP_FRACT, "", sep='') # Info about the thesholds involved in selecting genes. Should probably go in each file as a title, or in the filename
          
          print("Generating the plot_genes_jitter figure with plot_trend...")
          ourPdf("tmp.fig.monocle.trend.1.pdf")
          print(monocle::plot_genes_jitter(cool_genes, grouping="SAMPLEGROUP", ncol=2, color_by="SAMPLEGROUP", plot_trend=TRUE))
          dev.off()

          ourPdf("tmp.take2.fig.monocle.trend.1.pdf")
          monocle::plot_genes_jitter(cool_genes, grouping="SAMPLEGROUP", ncol=2, color_by="SAMPLEGROUP", plot_trend=TRUE)
          dev.off()
          
          print("Generating the plot_genes_jitter figure WITHOUT plot_trend...")
          ourPdf("tmp.fig.monocle.no_trend.1.pdf")
          print(monocle::plot_genes_jitter(cool_genes, grouping="SAMPLEGROUP", ncol=2, color_by="SAMPLEGROUP", plot_trend=FALSE))
          dev.off()
     }
     #stopifnot(1==2)

     ordering_gene_names <- intersect(sig_gene_names, oft_exp_gn)    #Only use genes are detectably expressed in a sufficient number of cells
     stopifnot(length(ordering_gene_names) >= 2) # Gotta be 2+ for this to make any sense!
     
     print("Set ordering filter...")
     HSMM_FILT <- setOrderingFilter(HSMM, ordering_gene_names)
     
     print("Reduce dimension...")
     HSMM_RED <- reduceDimension(HSMM_FILT, use_irlba=FALSE)

     # ==================== ==================== ==================== ==================== ====================
     NUM_PATHS_FOR_SPANNING_TREE <- 2 # Turns out, this doesn't really matter, it just changes the colors!
     # See if we can find a valid number for "num_paths"...
     found_valid_num_paths = FALSE
     for (i in 5:1) { # Start at 5, this is basically arbitrary though. Crashes if this number is too... high? Unclear.
          NUM_PATHS_FOR_SPANNING_TREE = i
          result <- tryCatch({
               print(paste("[Monocle] Ordering cells with num_paths=", NUM_PATHS_FOR_SPANNING_TREE, "...", sep=''))
               HSMM2 <- monocle::orderCells(HSMM_RED, num_paths=NUM_PATHS_FOR_SPANNING_TREE, reverse=TRUE)
               print(paste("[Monocle] Ordering cells with num_paths=", NUM_PATHS_FOR_SPANNING_TREE, " was SUCCESSFUL.", sep=''))
               found_valid_num_paths = TRUE
          }, error=function(e) {
               print(paste("[Monocle] could not run 'orderCells' with num_paths=", NUM_PATHS_FOR_SPANNING_TREE, ". Error message is: ", e, sep=''))
          }, finally={
               # Nothing to do afterward
          })
          if (found_valid_num_paths) { break; } # Get out of this for loop!
     }

     #HSMM2 <- monocle::orderCells(HSMM_RED, num_paths=NUM_PATHS_FOR_SPANNING_TREE, reverse=TRUE)
     NAME_SIZE <- 3.0
     print("Generating the plot_spanning_tree figures...")
     ourPdf(paste("tmp.fig.monocle.tree.detailed.", "paths=", NUM_PATHS_FOR_SPANNING_TREE, ".pdf", sep=''))
     print(monocle::plot_spanning_tree(HSMM2, color_by="SAMPLEGROUP", show_cell_names=TRUE, cell_name_size=NAME_SIZE))
     dev.off()

     ourPdf(paste("tmp.fig.monocle.tree.detailed.nobackbone.", "paths=", NUM_PATHS_FOR_SPANNING_TREE, ".pdf", sep=''))
     print(monocle::plot_spanning_tree(HSMM2, color_by="SAMPLEGROUP", show_cell_names=TRUE, cell_name_size=NAME_SIZE, show_backbone=F))
     dev.off()

     ourPdf(paste("tmp.fig.monocle.tree.detailed.nobackbone.", "paths=", NUM_PATHS_FOR_SPANNING_TREE, ".nonames.pdf", sep=''))
     print(monocle::plot_spanning_tree(HSMM2, color_by="SAMPLEGROUP", show_backbone=F))
     dev.off()

     ourPdf(paste("tmp.fig.monocle.tree.state.", "paths=", NUM_PATHS_FOR_SPANNING_TREE, ".nonames.pdf", sep=''))
     print(monocle::plot_spanning_tree(HSMM2, color_by="State", show_cell_names=F))
     dev.off()

     ourPdf(paste("tmp.fig.monocle.tree.state.nobackbone.", "paths=", NUM_PATHS_FOR_SPANNING_TREE, ".nonames.pdf", sep=''))
     print(monocle::plot_spanning_tree(HSMM2, color_by="State", show_cell_names=F, show_backbone=F))
     dev.off()

     # ==================== ==================== ==================== ==================== ====================
     HSMM_filtered <- HSMM[oft_exp_gn, pData(HSMM)[["SAMPLEGROUP"]] != "useless_group_or_something_basically_just_include_everything"]

     stopifnot(ncol(HSMM_filtered) > 0) # Too few SAMPLES
     stopifnot(row(HSMM_filtered) > 0)  # Too few GENES
     
     my_gene_names <- tail(head(rownames(HSMM_filtered), n=10), n=5)
     #my_genes <- row.names(subset(fData(HSMM_filtered), rownames(HSMM_filtered) %in% my_gene_names))   #my_genes <- row.names(subset(fData(HSMM_filtered), gene_short_name %in% c("ENSMUSG00000000001", "ENSMUSG00000000088", "ENSMUSG00000000194")))
     cds_subset <- HSMM_filtered[my_gene_names,]
     print(dim(cds_subset@assayData$exprs)) # just FYI...

     figPstime = "tmp.fig.monocle.pseudo.pdf"
     result <- tryCatch({
          print("Generating the plot_genes_in_pseudotime figures...")
          ourPdf(figPstime)
          print(monocle::plot_genes_in_pseudotime(cds_subset, color_by="Hours")) # may need a column labeled Hours or Pseudo-time or something?
     }, error=function(e) {
          system(paste("/bin/rm ", figPstime, sep='')) # Delete the broken image / figure
          print(paste("[Monocle] was (as expected...) unable to plot a 'plot_genes_in_pseudotime' figure: Unable to generate plot_genes_in_pseudotime plot for some reason. Message is:", e))
     }, finally={
          dev.off() # Always close the output device afterward
     })
}


>>>>>>> b3752661da232c4b0bd72fdd3b3252226bd15be4
