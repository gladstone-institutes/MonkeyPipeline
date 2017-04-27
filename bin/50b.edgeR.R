require(data.table); require(edgeR); require(gplots)
require(methods) # <-- just in case we use RScript to run this
#source("/data/home/alexgw/TimeForScience/Lab_Code/R/AGW/agwUtil.R"); # Has "pairs.agw" in it
#options(error=traceback) # useful even non-interactively
#options(error=recover)
options(stringsAsFactors=FALSE)

SUBREAD_FILE <- Sys.getenv("SUBREAD_COUNTS_FILE") # <-- adds 'SUBREAD_COUNTS_FILE' as a new variable in the namespace. Requires the environment variable to exist!
# Note: subread file CAN be any of these: 1) uncompressed, 2) gzipped, 3) bzip2'd! All these are transparently handled by "gzfile(...)"
GROUP_LIST         <- Sys.getenv("GROUP_LIST") # <-- separated by MAJOR_DELIM, example:  ctrl:experiment:drugX:drugY
BAM_BY_GROUP       <- Sys.getenv("BAM_LIST_BY_GROUP") # <-- groups are separated by MAJOR_DELIM, files by MINOR_DELIM. e.g. ctrl1.bam~ctrl2.bam:exp1.bam~exp2.bam
MAJOR_DELIM        <- Sys.getenv("MAJOR_DELIM") # Separates GROUPS (e.g. "Group1:Group2"). Probably a ':' colon.
MINOR_DELIM        <- Sys.getenv("MINOR_DELIM") # Separates REPLICATES within a group (e.g. "Drug1,Drug2,Drug3"). Probably a comma.

ROUND_TO           <- 2 # Round CPM and FPKM decimal places to this amount.
FPKM_COL_TEXT      <- "Naive_FPKM_across_all_genes:"
CPM_NAIVE_COL_TEXT <- "Naive_CPM_across_all_genes:"
CPM_COL_TEXT       <- "Normalized_CPM:"              # <-- Only for FILTERED genes that are used in the Diff Expression calculation. Omits the LOW and HIGH genes.

FINAL_DE_FILE      <- "EDGER.all.table.txt" # <-- note: this is HARD-CODED here and is also used in 50a.edgeR.pl (to create the variable "$finalDiffExpressionFile". This value should be identical to the filename in that variable! Do NOT change this variable without changing the one in the perl script 50a.edgeR.pl!

MAX_DE_COMPARISONS_TO_DO_BEFORE_SKIPPING_THE_REST <- (12*11)+1 # if there are more than 12 groups, stop doing all the pair-wise comparisons. It's just too many!

# Maybe these should be acquired from Sys.getenv?
#NUM_GENES_REQ_TO_FILTER_OUT_HIGH_VALUES <- 1000 # If there are fewer than 1000 genes/features, then we are not looking at a whole-genome sample and should NOT filter out "too high expression" values
NUM_EXPR_GENES_REQ_TO_FILTER_HI <- 5000 # If there are fewer than 5000 genes/features with expression, then we are not looking at a whole-genome sample and should NOT filter out "too high expression" values. It's probably microRNA data then.
# Maybe this should instead be a cutoff like "if the HIGH threshold filters more than 1% of genes, then keep them"

RAW_COUNTS_MIN_CUTOFF <- 5 # require at least 2 samples have 5+ counts in order to keep that gene/feature as having ANY expression. NOT a CPM cutoff!
#CPM_MIN_CUTOFF <- 0.5   # Conservative! EdgeR recommends 100. CPM = COUNTS PER MILLION. So 5 reads per 10,000,000 reads.
CPM_MAX_CUTOFF <- 20000  # Should filter out basically only mitochondrial genes and possibly a couple others at this point.  # Remove any genes where there aren't at least 2 samples with less than fraction of total reads for a sample. Only applies for whole-genome experiments where the number of genes is at least NUM_GENES_REQ_TO_FILTER_OUT_HIGH_VALUES. A gene with 10,000 CPM = 1% of total reads for that sample. 100,000 CPM = 10% of total reads, etc.

HEATNAME <- "edger_heatmap"

AGWDEBUG = Sys.getenv("AGWDEBUG")
if (AGWDEBUG != "" && as.numeric(AGWDEBUG)>0) {
     # ********************* HOW TO RUN THIS DEBUGGING CODE ************************
     # From the terminal, run this:     $ export AGWDEBUG=1 && R
     # Now, inside R,     run this:     R>  source("/work/monkey_development/bin/50b.edgeR.R")
     # *****************************************************************************
     # This is DEBUGGING CODE for testing without a subread file.
     if (grepl("[/]monkey", getwd())) { print(">>>>>>>>> Hey you! Don't run this test code from the MONKEY directory, because it will pollute it with a bunch of test files! Make a test directory (that is ***not*** a subdirectory of the monkey directory, and doesn't start with 'monkey') and run this test code in there instead!!!!"); stopifnot(1==2); }

     GROUP_LIST = "A:B:C"
     BAM_BY_GROUP = "aa.1.bam~aa.2.bam~aa.3.bam:bb.1.bam~bb.2.bam~bb.3.bam:cc.1_sample.bam~cc.2.human_or_something.bam"
     MAJOR_DELIM = ":"
     MINOR_DELIM = "~"
     SUBREAD_FILE = "test.subread"
     NR   <- 250 # num fake rows. Needs to be at least 100 or so or it fails. Always fails for N = 10!
     set.seed(123); # set rand seed
     SPOS <- floor(runif(NR)*10000)
     LEN  <- 100+floor(runif(NR)*1000)
     FAKE_EXPR_BASE = floor(rexp(NR)*100)
     FAKE_NAMES = sample( paste0("ENSFAKE", sprintf("%06d", seq(NR))) ) # Note that we RANDOMLY re-order it so they aren't in alphabetical order!
     FAKE_ADJ = c(0, 99, runif(NR-2)*0.90 + 0.10)  # randomly adjust it. Note the 0 (low) and the 99999 (too high) expression items
     A_ADJ = 1.0 * FAKE_ADJ
     B_ADJ = 1.8 * FAKE_ADJ
     C_ADJ = 1.5 * FAKE_ADJ
     x <- data.table("Geneid"=FAKE_NAMES, "Chr"=rep("chrFake", NR), "Start"=SPOS, "End"=SPOS+LEN, "Strand"=rep("+",NR), "Length"=LEN
                     , "aa.1.bam"=floor((runif(NR)*0.20+0.40) * FAKE_EXPR_BASE * A_ADJ)
                     , "aa.2.bam"=floor((runif(NR)*0.20+0.40) * FAKE_EXPR_BASE * A_ADJ)
                     , "aa.3.bam"=floor((runif(NR)*0.20+0.40) * FAKE_EXPR_BASE * A_ADJ)
                     , "bb.1.bam"=floor((runif(NR)*0.20+0.40) * FAKE_EXPR_BASE * B_ADJ)
                     , "bb.2.bam"=floor((runif(NR)*0.20+0.40) * FAKE_EXPR_BASE * B_ADJ)
                     , "bb.3.bam"=floor((runif(NR)*0.20+0.40) * FAKE_EXPR_BASE * B_ADJ)
                     , "cc.1_sample.bam"            =floor((runif(NR)*0.20+0.40) * FAKE_EXPR_BASE * C_ADJ)
                     , "cc.2.human_or_something.bam"=floor((runif(NR)*0.20+0.40) * FAKE_EXPR_BASE * C_ADJ)   )
     write.table(x, file=SUBREAD_FILE, quote=F, sep="\t", row.names=FALSE, col.names=TRUE)

     # Manually load a subread file from a different project
     #GROUP_LIST = "N_MG132:NH_Neuron_noHS:A_DMSO:A3_Astro_3hrHS:A_MG132:N_DMSO:N3_Neuron_3hrHS:AH_Astro_noHS"
     #BAM_BY_GROUP = 'N_MG132.1_mm9.chr_q30.bam~N_MG132.4_mm9.chr_q30.bam~N_MG132.3_mm9.chr_q30.bam~N_MG132.2_mm9.chr_q30.bam:NH_Neuron_noHS.1_mm9.chr_q30.bam~NH_Neuron_noHS.3_mm9.chr_q30.bam~NH_Neuron_noHS.2_mm9.chr_q30.bam:A_DMSO.4_mm9.chr_q30.bam~A_DMSO.1_mm9.chr_q30.bam~A_DMSO.3_mm9.chr_q30.bam~A_DMSO.2_mm9.chr_q30.bam:A3_Astro_3hrHS.1_mm9.chr_q30.bam~A3_Astro_3hrHS.3_mm9.chr_q30.bam~A3_Astro_3hrHS.2_mm9.chr_q30.bam:A_MG132.1_mm9.chr_q30.bam~A_MG132.4_mm9.chr_q30.bam~A_MG132.3_mm9.chr_q30.bam~A_MG132.2_mm9.chr_q30.bam:N_DMSO.4_mm9.chr_q30.bam~N_DMSO.1_mm9.chr_q30.bam~N_DMSO.3_mm9.chr_q30.bam~N_DMSO.2_mm9.chr_q30.bam:N3_Neuron_3hrHS.1_mm9.chr_q30.bam~N3_Neuron_3hrHS.3_mm9.chr_q30.bam~N3_Neuron_3hrHS.2_mm9.chr_q30.bam:AH_Astro_noHS.1_mm9.chr_q30.bam~AH_Astro_noHS.3_mm9.chr_q30.bam~AH_Astro_noHS.2_mm9.chr_q30.bam'
     #MAJOR_DELIM = ":"
     #MINOR_DELIM = "~"
     #SUBREAD_FILE = "a150.txt"
}

OUT_PREFIX             <- "EDGER.OUT."
NORM_COUNT_FILENAME    <- paste0(OUT_PREFIX, "000.normalized_counts_per_million.txt")
FPKM_ONLY_FILENAME     <- "fpkm_naive.txt" # <-- Do NOT change this without also changing the monocle code! This literal filename will be searched for!
RAW_COUNTS_FILE        <- "raw_counts.txt" # <-- Do NOT change this without also changing the 51b.summary.R code! This literal filename will be searched for!
NORM_STAT_SUMMARY_FILE <- paste0("EDGER.STATS.normalization_summary_stats.txt")

stopifnot(     nchar(SUBREAD_FILE) >= 1 && file.exists(SUBREAD_FILE))
stopifnot(       nchar(GROUP_LIST) >= 1)
stopifnot(nchar(BAM_BY_GROUP) >= 1)
stopifnot(      nchar(MAJOR_DELIM) >= 1)
stopifnot(      nchar(MINOR_DELIM) >= 1)

pgsub <- function(...) { (gsub(..., perl=T, ignore.case=T)) } # shorthand for perl=T and ignore.case=T

rescue_us_if_there_are_no_replicates <- function(yyy) {
     if (all(is.na(yyy$trended.dispersion))) {
          print("WARNING: Looks like some samples have no replicates, so we have to do a workaround...")
          print("WARNING: ******** This is probably only ok if you: 1) Actually expect some group to have zero replicates")
          print("WARNING: ********                               or 2) Are running monkey on --test2 or other test data (those don't have replicates either).")
          system("echo 'WARNING_50b.edgeR.R found some conditions with no replicates, and set the dispersion MANUALLY. This may or may not be OK.' > WARNING_SOME_SAMPLES_HAD_ZERO_REPLICATES_AND_DISPERSION_WAS_SET_MANUALLY.txt")
          # Handle the "no replicates" condition gracefully:
          # From the EdgeR manual:  http://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
          # As a last resort (don't trust these P-values too much!)...
          # ""Simply pick a reasonable dispersion value, based on your experience with similar data,""
          # ""and use that for exactTest or glmFit. Typical values for the common BCV (squareroot-dispersion)""
          # ""for datasets arising from well-controlled experiments are 0.4 for human""
          # ""data, 0.1 for data on genetically identical model organisms or 0.01 for technical replicates.""
          # ""Here is a toy example with simulated data:""
          placeholder_no_replicates_dispersion <- 0.40 # Hard-coded to be decent-ish for human and conservative for other species.
          yyy$trended.dispersion[ is.na(yyy$trended.dispersion) ] <- (placeholder_no_replicates_dispersion**2)
     }
     return(yyy)
}

groupNames.vec <- unlist(strsplit(GROUP_LIST, MAJOR_DELIM, fixed=TRUE)) # fixes things like "74", which really needs to be a name like 'X74' for R's purposes
bamByGroup.vec <- unlist(strsplit(BAM_BY_GROUP, MAJOR_DELIM, fixed=TRUE)); stopifnot(length(groupNames.vec) == length(bamByGroup.vec))
nGroups        <- length(groupNames.vec)
#bamsByGroup.list = list() # List: name = group name, value = vector of bam filenames in this group
bamGrpUnordered.list = list() # List: name = BAM, value = GROUP name
for (i in seq_along(groupNames.vec)) {
     grName             <- groupNames.vec[i]
     bamsInThisGroupStr <- bamByGroup.vec[i]
     #if (grName != make.names(grName)) { stop(paste("ERROR: one of your group names probably started with a number or something, or is otherwise unsuitable for edgeR's naming requirements. Here is the offending name: <", grName, ">, and this is what it looked like when it was passed through the 'make.names' function: <", make.names(grName), ">. It should have been unchanged, but as you can see, it was changed! Fix this.", sep='')) }
     temp.bams.vec <- unlist(strsplit(bamsInThisGroupStr, MINOR_DELIM, fixed=TRUE)) # <-- minor delim, NOT the major delim!
     temp.bams.vec <- pgsub(".*[/]", "", temp.bams.vec); # <-- Stripping full paths from bam names...
     stopifnot(all(grepl(".bam$", temp.bams.vec, ignore.case=TRUE, perl=TRUE))) # make sure all the bam filenames end in '.bam'!
     stopifnot(length(temp.bams.vec) >= 1) # At least one bam per group, please!!
     #bamsByGroup.list[[thisGroupName]] = temp.bams.vec # Save all the individual bams!
     for (bam in temp.bams.vec) {
          bamGrpUnordered.list[[bam]] <- grName # save each BAM's group
     }
}

# bamGrpUnordered.list is the data structure that you MOST LIKELY want to use from this point onward.
# It has the BAM filename (NOT the full path!) as the key, and the group (as a string) as a value

if (TRUE) { #!exists("dfilt.counts.mat")) {
     ddd <- read.table(gzfile(SUBREAD_FILE), sep="\t", header=TRUE
                       , check.names=F, stringsAsFactors=FALSE, fill=FALSE,
                       , row.names=1, quote='', comment.char='#') # <-- comment.char is important: the first line of a subread file is a '#' comment (can also use nskip = 1)
     #print("Note: we treat lines starting with '#' as comments (because subread generates a comment line)!")
     colnames(ddd) <- basename(colnames(ddd)) # Strip the full paths---only keep the filename
     all.dt  <- data.table::data.table(ids=rownames(ddd), original.data=ddd) # all.dt includes the ANNOTATION ROWS
     firstDataColIdx <- (1 + which(toupper(colnames(ddd)) == "LENGTH")); print("Assuming that the last column of annotation is named 'Length' (case-insensitive)...")
     all.counts.mat   <- data.matrix(ddd[, firstDataColIdx:ncol(ddd), drop=F]) # all.counts.mat = only DATA
     all_gene_lengths <- as.numeric(ddd$Length) ; names(all_gene_lengths) <- rownames(ddd) # Vector
     # All.Counts.Mat does NOT have the ANNOTATION ROWS that subread provides in cols1 1-5

     # ============ THIS IS THE PART WHERE WE FIGURE OUT WHICH COLUMN HEADERS IN THE 'SUBREAD' COUNTS FILE
     # ============ ACTUALLY MATCH UP TO WHICH BAM FILES AND GROUPS
     #print("Unordered list:")
     #print(names(bamGrpUnordered.list))
     #print("Then if it is ordered:")
     #print(sort(names(bamGrpUnordered.list)))
     #print("Length:")
     #print(length(bamGrpUnordered.list))
     bamGrpOrdered.vec <- vector("character", length(bamGrpUnordered.list))
     if (ncol(all.counts.mat) != length(bamGrpUnordered.list)) {
          stop(paste0("ERROR: Serious problem in the input files to edgeR---it looks like your 'all counts' matrix (probably generated by subread featureCounts had ", ncol(all.counts.mat), " columns of data, but the number of input BAM files was ", length(bamGrpUnordered.list), "---this can happen if you are trying to MANUALLY generate the subread file due to some issues in the input data. This is NOT a recoverable problem---you must fix the input file!"))
     }
     for (ib in seq_along(colnames(all.counts.mat))) {
          thisBamName <- colnames(all.counts.mat)[ib]
          bamRE       <- paste0("^", gsub("\\.","[.]", thisBamName), "$") # The gsub turns all '.' (which means "any character" in perl-regexp land into LITERAL [.] only-match-a-period). We also use ^ and $ to make sure that the FULL STRING matches
          bamIsAt     <- grep(bamRE, x=names(bamGrpUnordered.list), ignore.case=TRUE, perl=TRUE)
          if (length(bamIsAt) > 1) {       stop(paste0("50b.edgeR.R: ERROR: this bam file (named <", thisBamName, ">) was found at MULTIPLE LOCATIONS in the 'unordered' list! This should be impossible, even if one name is a strictly longer version of another... debug this!")) }
          else if (length(bamIsAt) == 0) { stop(paste("50b.edgeR.R: ERROR: this bam file (named <", thisBamName, ">) was found at ZERO LOCATIONS in the 'unordered' list! This should be impossible! The unordered list contained the following (commad-separated) names: ", paste0(names(bamGrpUnordered.list), collapse=','), "  .")) }
          #print(paste("Checking the bam indexes for bam named ", thisBamName, " in these groups: ", names(bamGrpUnordered.list), sep=''))
          #print(paste0("This bam name: ", thisBamName, " assigned to ib at ", ib))
          bamGrpOrdered.vec[ib] <- bamGrpUnordered.list[[thisBamName]]
     }
     if (!all(bamGrpOrdered.vec != "")) {
          print(paste0("50b.edgeR.R: ERROR: FAILURE in edgeR: Uh oh, some of your 'bamGrpOrdered.vec' was (somehow) blank. Here is the FULL list, which consists of ", length(bamGrpOrdered.vec), " items, of which a total of ", sum(bamGrpOrdered.vec == ""), " were apparently blank. The list: ")); print(bamGrpOrdered.vec)
          stopifnot(all(bamGrpOrdered.vec != "")) # "NO BLANK ENTRIES ALLOWED PLEASE"
     }
     # ============ DONE MATCHING UP BAM FILES AND GROUPS =================
     Grps      <- factor(bamGrpOrdered.vec)
     dc.table  <- data.table::data.table(ids=rownames(all.counts.mat), original.counts=all.counts.mat)
     #print("Dimensions of the input data are:"); print(dim(all.counts.mat))
     low.ok.vec <- rowSums(all.counts.mat >= RAW_COUNTS_MIN_CUTOFF) >= 2 # at least 2+ samples had raw counts above min cutoff
     
     we_should_filter_high_cpm <- (sum(low.ok.vec) >= NUM_EXPR_GENES_REQ_TO_FILTER_HI) # Are there enough actually-expressed genes for us to think this is probably a whole-genome experiment where we should also filter out highly-expressed genes?
     if (!we_should_filter_high_cpm) {
          hi.ok.vec <- rep(TRUE, times=nrow(all.counts.mat)) # Just keep ALL of them.
          print(paste0("We apparently are not doing a whole-genome mRNA analysis (there are < ", NUM_EXPR_GENES_REQ_TO_FILTER_HI, " genes/features being tested), so we should NOT filter out genes with high expression."))
     } else {
          hi.ok.vec  <- rowSums(edgeR::cpm(all.counts.mat) <= CPM_MAX_CUTOFF) >= 2 # See whether at least 2+ samples weren't "too high"
          if (sum(!hi.ok.vec) > 25) {
               hi.ok.vec <- rep(TRUE, times=nrow(all.counts.mat))
               we_should_filter_high_cpm <- FALSE
               print("There were over 25 (!!!) genes with 'too high' expression, which probably means we aren't actually dealing with regular mRNA whole-genome data. So we should definitely NOT be disqualifying genes for having 'too high' expression in this case. Skipping the 'too high' disqualification. Normally we expect to see 3 or 4 mitochondrial genes disqualified here AT MOST. CONCLUSION: we will not disqualify ANY gene for having 'too much' expression. (We can still disqualify them for insufficient expression, however.).")
          }
     }
     
     keep.all.vec <- low.ok.vec & hi.ok.vec
     print(paste0("Keeping a total of ", sum(keep.all.vec), " genes out of "
                 , length(keep.all.vec), ". ", sum(!low.ok.vec), " had insufficient expression, and "
                 , sum(!hi.ok.vec), " had too much."))
     dfilt.counts.mat <- all.counts.mat[keep.all.vec, , drop=FALSE]
     #print("'dfilt.counts.mat' is FILTERED 'all.counts.mat' with only 'reasonable CPM' genes.")
}

if (TRUE) { #!exists("reord.counts.and.cpm")) {
     #print("EdgeR DGEList processing now...")
     y.all              <- edgeR::DGEList(counts=all.counts.mat, group=Grps) # y.all has ALL GENES, not just the "has a reasonable amount of expression" genes.
     stopifnot(all(rownames(y.all$counts) == names(all_gene_lengths)))
     rpkm.all.mat       <- round(edgeR::rpkm(y.all, gene.length=all_gene_lengths, normalized.lib.sizes=TRUE, log=FALSE), ROUND_TO)
     cpm.naive.all.mat  <- round(edgeR::cpm( y.all, normalized.lib.sizes=TRUE, log=FALSE), ROUND_TO)
     
     # ============= WRITE AN OUTPUT FILE FOR USE BY THE 'MONOCLE' (Cole Trapnell's project) SINGLE CELL ANALYSIS PACKAGE
     rpkm_write.dt   <- data.table::data.table(ids=rownames(rpkm.all.mat), rpkm.all.mat) # This is specially formatted for use by the "Monocle" package later!!
     #                                          Note about the regexp below: SOMETHING.### <-- note: everything after the number gets deleted. See example on line below's comment.
     data.table::setnames(rpkm_write.dt, old=names(rpkm_write.dt), new=gsub("([^.]*[.]\\d+).*", "\\1", names(rpkm_write.dt), perl=T, ignore.case=T) ) # <-- Replaces something like "CONTROL.24_mouse_mm9.something.bam with just "CONTROL.24"
     write.table(rpkm_write.dt, file=FPKM_ONLY_FILENAME, row.names=F, quote=F, sep="\t")
     # ============ DONE WRITING THE OUTPUT FILE FOR THE 'MONOCLE' PACKAGE

     # ============ NOW JUST WRITE THE COUNTS OUTPUT  -- actually, people should just use the subread file for this
     #counts_write.dt   <- data.table::data.table(ids=rownames(all.counts.mat), all.counts.mat) # This will be used elsewhere
     ##                                          Note about the regexp below: SOMETHING.### <-- note: everything after the number gets deleted. See example on line below's comment.
     #data.table::setnames(counts_write.dt, old=names(counts_write.dt), new=gsub("([^.]*[.]\\d+).*", "\\1", names(counts_write.dt), perl=T, ignore.case=T) ) # <-- Replaces something like "CONTROL.24_mouse_mm9.something.bam with just "CONTROL.24"
     #write.table(counts_write.dt, file=RAW_COUNTS_FILE, row.names=F, quote=F, sep="\t")
     # ============ DONE WRITING THE COUNTS
     
     # AFTER copying it to "rpkm_formatted_to_write_dt", Fix the names for the rpkm.dt
     rpkm.dt            <- data.table::data.table(ids=rownames(rpkm.all.mat), rpkm.all.mat)
     data.table::setnames(rpkm.dt, names(rpkm.dt)[2:ncol(rpkm.dt)], paste(FPKM_COL_TEXT, names(rpkm.dt)[2:ncol(rpkm.dt)], sep=" ") )

     cpm.naive.all.dt   <- data.table::data.table(ids=rownames(cpm.naive.all.mat), cpm.naive.all.mat)
     data.table::setnames(cpm.naive.all.dt, names(cpm.naive.all.dt)[2:ncol(cpm.naive.all.dt)], paste(CPM_NAIVE_COL_TEXT, names(cpm.naive.all.dt)[2:ncol(cpm.naive.all.dt)], sep=" ") )
     
     y.filt       <- edgeR::DGEList(counts=dfilt.counts.mat, group=Grps) # y.filt only has the 'keep' genes from above
     cpm.filt.mat <- round(edgeR::cpm(y.filt, normalized.lib.sizes=TRUE, log=FALSE), ROUND_TO) # rounded values
     
     #print("Note that cpm.filt.mat contains ONLY the non-NA (OK CPM) genes!")
     # Above: a really non-obvious way of computing cpm for a 'DGEList'.
     #        Can look VERY weird if the norm.factors in y.filt$samples is far from 0,
     #        but if you plot a histogram you'll see that this number actually
     #        does a better normalization
     #print("EdgeR calcNormFactors...")
     y.filt <- edgeR::calcNormFactors(y.filt)
     # Above: "normalizes for RNA composition by finding a set of scaling factors
     #        for the library sizes that minimize the log-fold changes between
     #        the samples for most genes."
     #print("EdgeR GLM...")
     
     y.filt <- edgeR::estimateGLMCommonDisp(y.filt)
     #print("EdgeR GLM common is done...")
     y.filt <- edgeR::estimateGLMTrendedDisp(y.filt)
     #print("EdgeR GLM trended is done...")
     y.filt <- rescue_us_if_there_are_no_replicates(y.filt) # Has to come right AFTER estimateGLMTrendedDisp
     y.filt <- edgeR::estimateGLMTagwiseDisp(y.filt)
     #print("EdgeR GLM tagwise disp is done...")
     #print("Making cpm data table and printing it to a file...")
     cpm.CHAR.mat <- matrix(as.character(cpm.filt.mat),nrow=nrow(cpm.filt.mat), dimnames=list(rownames(cpm.filt.mat), colnames(cpm.filt.mat)))
     cpm.dt <- data.table::data.table(ids=rownames(cpm.filt.mat), cpm.CHAR.mat) # Yes, this has to become a character matrix, because we will use things like "NA_LOW_EXPRESSION" and "NA_HIGH_EXPRESSION" below, which don't work in a numeric matrix
     data.table::setnames(cpm.dt, names(cpm.dt)[2:ncol(cpm.dt)], paste(CPM_COL_TEXT, names(cpm.dt)[2:ncol(cpm.dt)], sep=" ") )

     reord.all <- cbind(all.dt, "NORM_PROCEDURE_PLACEHOLDER"="") # <-- insert a placeholder where we'll put descriptive text below
     reord.all <- merge(reord.all, cpm.dt, by="ids", all.x=T)            # Reorder so that we write the data in a predictable order! Also, include ALL the ids (from all.dt) --- cpm has FEWER rows than all.dt!
     z <- reord.all
     cpm_col_indexes.vec <- grep(paste("^", CPM_COL_TEXT, sep=''), colnames(reord.all), perl=T, ignore.case=T)
     if (sum(!low.ok.vec) > 0) {
          too_low_names = names(low.ok.vec[low.ok.vec == FALSE])
          reord.all[too_low_names, cpm_col_indexes.vec := "NA--low" , with=F] # <- ":=" <--This is the SUPER WEIRD way to assign values in data.table to specific elements only
     }

     if (sum(!hi.ok.vec) > 0) {
          too_high_names = names(hi.ok.vec[hi.ok.vec == FALSE])
          reord.all[too_high_names , cpm_col_indexes.vec := "NA--high" , with=F] # <- ":=" <-- This is the SUPER WEIRD way to assign values in data.table to specific elements only
          # NB: You can't reassign these unless the column classes are already non-numeric
     }

     stopifnot(0 == sum(is.na(cpm.naive.all.dt)))
     stopifnot(0 == sum(is.na(rpkm.dt)))

     reord.all <- cbind(reord.all, "NOTE: the 'naive' CPM and FPKM values to the right are calculated using ALL genes---however, only a subset of genes within a certain expression range are actually used to calculate differential expression. See the 'Normalized_CPM' columns at left. Note that the normalized values are calculated WITHOUT the 'NA--high' and 'NA--low' genes, so the two CPM values are expected to be different."="") # <-- insert a placeholder where we'll put descriptive text below
     reord.all <- merge(reord.all, cpm.naive.all.dt, by="ids", all.x=T)
     reord.all <- merge(reord.all, rpkm.dt, by="ids", all.x=T)

     # Set some descriptive text
     if (we_should_filter_high_cpm) {
          describe_cutoffs <- paste(" there are not both 2+ samples (in any experimental groups) with ", RAW_COUNTS_MIN_CUTOFF, "+ raw reads AND ALSO 2+ samples with a naive CPM value of less than ", CPM_MAX_CUTOFF, " (this generally excludes 2-5 highly-expressed mitchondrial genes).", sep='')
     } else {
          describe_cutoffs <- paste(" there are not 2+ samples (in any experimental groups) with ", RAW_COUNTS_MIN_CUTOFF, "+ raw reads.", sep='')
     }
     
     setnames(reord.all, old="NORM_PROCEDURE_PLACEHOLDER", new=paste("Note: Certain genes are omitted from the normalized CPM (counts per million). In this project, we omitted a total of ", sum(!low.ok.vec), " genes/features with insufficient expression and ", sum(!hi.ok.vec), " genes/features with extremely high expression. The omission procedure is: 1) CPM is calculated by edgeR. 2) Any genes/features are set to 'NA' values if ", describe_cutoffs, " 3) Then the CPM is recalculated with those genes removed. 4) As a result, all normalized CPM columns should sum to exactly 1 million.", sep=''))
     setnames(reord.all, old=colnames(reord.all), new=pgsub("^ids$", "ENSEMBL IDs", colnames(reord.all)))  # Make the column names nicer for the final output
     setnames(reord.all, old=colnames(reord.all), new=pgsub("^original[.]data[.](Chr|Start|End|Strand|Length)$", "\\1", colnames(reord.all)))
     setnames(reord.all, old=colnames(reord.all), new=pgsub("^Chr$", "Chromosome", colnames(reord.all)))
     setnames(reord.all, old=colnames(reord.all), new=pgsub("^(Chromosome|Start|End|Strand)$", "\\1 (for all exons in annotation)", colnames(reord.all)))
     setnames(reord.all, old=colnames(reord.all), new=pgsub("^Length$", "Feature length (base pairs)", colnames(reord.all)))
     setnames(reord.all, old=colnames(reord.all), new=pgsub("^original[-._]data[-._](.*[.]bam)$", "Raw counts (non-normalized) for \\1", colnames(reord.all))) # Make the column names nicer for the final output
     setnames(reord.all, old=colnames(reord.all), new=pgsub("^Normalized[-._]CPM[-._](.*[.]bam)$", "Normalized CPM (Counts Per Million fragments) for \\1", colnames(reord.all)))

     options(scipen=999) # <-- Never output SCIENTIFIC NOTATION!

     # reord.all[, grep("^(Naive|Normalized)", colnames(reord.all)), with=F]
     
     write.table(reord.all, file=NORM_COUNT_FILENAME, row.names=F, quote=F, sep="\t")
     stopifnot(nrow(all.dt) == nrow(rpkm.dt))
     stopifnot(nrow(all.dt) == nrow(reord.all))
}

write.table(y.filt$samples, file=NORM_STAT_SUMMARY_FILE, col.names=NA, row.names=T, quote=F, sep="\t")
#print(y.filt$samples); # Has counts and stuff, useful to check this out
#print(summary(y.filt$samples$norm.factors)); # These numbers should be around 1-ish!

require(matrixStats) # for "rowMaxs(...)" / "rowMins(...)"

# ========= GENERATE A HEATMAP OF THE 1000 'TOP' GENES ============
REQ_CPM_DIFF <- 10 # Show in heatmap: genes with at least this ABSOLUTE CPM diff. between some two samples, to omit the non-expressed genes
REQ_CPM_RATIO <- 2 # Show in heatmap: genes with this ratio between two samples
HEATMAP_TOP_N <- 1000 # Show this many genes at most, ordered by fold change
OK.vec <- c()
while (sum(OK.vec) < 2) {
     maxDiff.vec  <- abs(matrixStats::rowMaxs(cpm.filt.mat) - matrixStats::rowMins(cpm.filt.mat))
     maxRatio.vec <- (1+matrixStats::rowMaxs(cpm.filt.mat)) / (1+matrixStats::rowMins(cpm.filt.mat))
     diffOK.vec  <- maxDiff.vec  >= REQ_CPM_DIFF
     ratioOK.vec <- maxRatio.vec >= REQ_CPM_RATIO
     OK.vec <- diffOK.vec & ratioOK.vec

     if (sum(OK.vec) >= 2 || REQ_CPM_RATIO <= 1.05) {
          break; # ok, great, this worked! otherwise, we'll have to relax the criteria for differential expression
     } else {
          REQ_CPM_RATIO = REQ_CPM_RATIO * 0.8 # relax the stringency...
     }
}

# From EdgeR docs:
#The function plotMDS draws a multi-dimensional scaling plot of the RNA samples in which
#distances correspond to leading log-fold-changes between each pair of RNA samples. The
#leading log-fold-change is the average (root-mean-square) of the largest absolute log-foldchanges
#between each pair of samples

# Probably should plot the heatmap with PSEUDOCOUNTS (currently we are not doing this). See line below for details.
# fixed.cpm.with.pseudocount <- cpm(d, prior.count=2, log=TRUE) # Recommended way to get values for a heatmap!

if (sum(OK.vec) < 2) {
     print("Ah, a terrible situation, we were not able to find ANY genes for edgeR to make a nice heatmap with, due to lack of differential expression and/or quantity of genes!")
     png(paste(HEATNAME, ".png", sep=''), width=2000, height=2000, res=144)
     plot(c(1), main="No genes met the required threshold in edgeR, so we did not generate a heatmap.")
     dev.off()
     pdf(paste(HEATNAME, ".pdf", sep=''), width=14, height=14)
     plot(c(1), main="No genes met the required threshold in edgeR, so we did not generate a heatmap.")
     dev.off()
} else {
     reord.cpm.filt.mat <- cpm.filt.mat[OK.vec,]
     reord.cpm.filt.mat <- reord.cpm.filt.mat[ order(rowMeans(reord.cpm.filt.mat), decreasing=TRUE), ]
     log2re.cpm.filt.mat <- log2(1+head( reord.cpm.filt.mat, HEATMAP_TOP_N)) # Top 1000 rows only!

     NAME_FOR_YES <- paste("Top ", nrow(log2re.cpm.filt.mat), " genes\n(selected by fold change)", sep='')
     NAME_FOR_NO  <- paste("All ", nrow(log2re.cpm.filt.mat), " genes", sep='')
     topNText = ifelse(nrow(reord.cpm.filt.mat) > HEATMAP_TOP_N, yes=NAME_FOR_YES, no=NAME_FOR_NO)
     labRowLabel <- paste("Y-axis:\n", topNText, "\nwith\n1) at least ", REQ_CPM_DIFF, " C.P.M. between\n    low and high\nand\n2) ", REQ_CPM_RATIO, "x fold change\n    i.e. ((max+1) / (min+1) >= ", REQ_CPM_RATIO, ")", sep='')

     png(paste(HEATNAME, ".png", sep=''), width=2000, height=2000, res=144)
     heatmap.2(log2re.cpm.filt.mat, trace="none", margin=c(20,15), density.info="histogram", denscol="black", labRow=labRowLabel, cexRow=0.8)
     dev.off()

     pdf(paste(HEATNAME, ".pdf", sep=''), width=14, height=14)
     heatmap.2(log2re.cpm.filt.mat, trace="none", margin=c(20,15), density.info="histogram", denscol="black", labRow=labRowLabel, cexRow=0.8)
     dev.off()
}
# ========= DONE GENERATING A HEATMAP OF THE 1000 'TOP' GENES ============

# If you want EVERY PAIRWISE COMPARISON in 'Grps' (which is what we want in this case), here's one way to do it!
interest.list = list()
for (i in seq_along(levels(Grps))) {
     for (j in seq_along(levels(Grps))) {
          if (i >= j) { next; }
          theName <- paste(levels(Grps)[i], "__vs__",levels(Grps)[j], sep='')
          theComp <- paste(levels(Grps)[i], " - " ,levels(Grps)[j], sep='')
          interest.list[[theName]] <- list(group=Grps, comparison=theComp)
     }
}

# NOT GENERATING ANOVA STUFF NOW
if (0 == 1) {
     # ============ ANOVA STUFF ==============
     aovP.vec <- rep(1.0, times=nrow(cpm.filt.mat))
     for (iii in seq_along(rownames(cpm.filt.mat))) {
          vv <- data.frame("GROUP"=Grps, "GENE"=cpm.filt.mat[iii,,drop=T])
          theAOV <- aov(vv$GENE ~ vv$GROUP)
          aovP.single <- summary(theAOV)[[1]][[1,"Pr(>F)"]] # This is seriously how you get a P-value from AOV.
          #print(summary.aov(theAOV))
          #print("Can print the TukeyHSD if you want to break up the AOV by group:")
          #TukeyHSD(theAOV)
          aovP.vec[iii] <- aovP.single
          if (iii > 100) { break; }
     }
     # Sean's "anova"-like edgeR code
     ## x <- read.delim("fileofcounts.txt",row.names="Symbol")
     ## group <- factor(c(1,1,2,2,3,3,4,4))
     ## y <- DGEList(counts=x,group=group)
     ## design <- model.matrix(~group)
     ## y <- estimateGLMCommonDisp(y,design)
     ## y <- estimateGLMTrendedDisp(y,design)
     ## y <- estimateGLMTagwiseDisp(y,design)
     ## fit <- glmFit(y,design)
     ## lrt <- glmLRT(fit,coef=2:4)
     ## w <- which(lrt$table$PValue <= 0.05)
}


for (i in seq_along(interest.list)) {
     if (i > MAX_DE_COMPARISONS_TO_DO_BEFORE_SKIPPING_THE_REST) {
          print(paste0("WARNING: edgeR wanted to do a HUGE number of pair-wise comparisons (in fact, it wanted to perform ", length(interest.list), ") but we have a hard-coded stopping number of only ", MAX_DE_COMPARISONS_TO_DO_BEFORE_SKIPPING_THE_REST, ", so we did the first ones and then skipped the remainder."))
          system("touch 0000_WARNING_edgeR_did_not_perform_all_pairwise_comparisons--there_were_too_many.touch.txt")
          break;
     }
     thename    <- names(interest.list)[[i]]
     comparison <- interest.list[[i]]$comparison
     thegroup   <- factor(interest.list[[i]]$group) # Should already be a factor, but just in case.
     print(paste("EdgeR: ", i, ": running this comparison: ", deparse(comparison), "", sep=''))

     design           <- stats::model.matrix(~0+factor(thegroup)) # or factor(thegroup)
     colnames(design) <- levels(factor(thegroup))
     contr <- limma::makeContrasts(CCC=comparison, levels=design)
     
     #yyy <- edgeR::DGEList(counts=dfilt.counts.mat[, OFINTEREST], group=thegroup[OFINTEREST])
     # Don't subset here, if you do, edgeR will fixate on outliers and not be able
     # to estimate the "true" variance for each gene.
     # You'll get signif-p-value genes where only ONE SAMPLE in a group was aberrant.
     yyy <- edgeR::DGEList(counts=dfilt.counts.mat, group=thegroup)

     yyy <- edgeR::calcNormFactors(yyy)        # "normalizes for RNA composition by finding a set of scaling factors for the lib sizes that minimize the log-fold changes between the samples for most genes."
     yyy <- edgeR::estimateGLMCommonDisp(yyy)  # Estimate COMMON dispersion. May be redunant when Trended is used.
     yyy <- edgeR::estimateGLMTrendedDisp(yyy) # <-- apparently this makes the COMMON and TAG-WISE disperson obsolete? In other words, you can just run this (and not estimateGLMCommonDisp or estimateGLMTagwiseDisp, and you get the same P-values, as far as I can tell)
     yyy <- rescue_us_if_there_are_no_replicates(yyy) # Has to come right AFTER estimateGLMTrendedDisp
     yyy <- edgeR::estimateGLMTagwiseDisp(yyy) # Estimate TAG-WISE dispersion. May be redundant with Trended

     #yyy <- edgeR::estimateDisp(yyy, design=design) # <-- this fails to handle the N = 1 case with NO REPLICATES
     
     fit   <- edgeR::glmFit(yyy, design)
     lrt   <- edgeR::glmLRT(glmfit=fit, coef=ncol(fit$design), contrast=contr)
     # Note: only one comparison at a time!!
     fin   <- data.table(ids=rownames(dfilt.counts.mat), lrt$table, "Adj_BH_FDR"=p.adjust(lrt$table$PValue, method="BH"))
     # NOTE: we should probably remove the "LR" and "CPM" here because those are
     #       misleading -- it's the average of ALL SAMPLES!
     #       Not anything to do with the specific p-value or anything.
     data.table::setnames(fin, names(fin)[2:ncol(fin)]
                          , paste(thename, names(fin)[2:ncol(fin)], sep="."))

     # Below (f2): Just get the IDS from the 'master' table, not all the other metadata
     #             (counts, etc). Those were already written to a file above.
     f2.everything <- merge(dc.table[,"ids",with=F], fin, by="ids", all.x=T)
     f2 <- f2.everything[,grepl("[.](logFC|PValue|Adj_BH_FDR)", colnames(f2.everything), perl=T),with=F] # Only get the fold change, P-value, and FDR columns
     f2      <- cbind("PLACEHOLDER"="", f2)
     setnames(f2, old="PLACEHOLDER", new=paste("Differential expression: ", thename, sep='')) # Spacer column
     setnames(f2, old=colnames(f2), new=pgsub("(.*)(__vs__)(.*)[-._]logFC$", "\\1\\2\\3: log2 fold change (values > 0 indicate that '\\1' is higher than '\\3')", colnames(f2)))
     setnames(f2, old=colnames(f2), new=pgsub("[.]PValue$", ": P-value (raw, not accounting for multiple testing)", colnames(f2)))
     setnames(f2, old=colnames(f2), new=pgsub("[.]Adj_BH_FDR$", ": FDR (false discovery rate, Benjamini-Hochberg method)", colnames(f2)))
     #statFilename <- paste("EDGER.FULL.", sprintf("%03d",i),".",thename,".txt", sep='')
     #write.table(f2, file=statFilename, col.names=NA, row.names=T, quote=F, sep="\t")
     subset.filename <- paste0(OUT_PREFIX, sprintf("%03d",i),".",thename,".txt")
     write.table(f2, file=subset.filename, col.names=TRUE, row.names=FALSE, quote=F, sep="\t")
}

print("Pasting all the various D.E. output files together into one larger file")
FILES_TO_PASTE <- paste0(OUT_PREFIX,"*") # <-- Example: "EDGER.OUT.*"

TEMP_COLLECTION_FILE <- "temp.paste.edgeR.everything.rename.me--if-you-see-this-file-on-the-filesystem-it-means-something-got-interrupted-OR-pasting-is-in-progress"

pasteCmd       <- paste("paste", FILES_TO_PASTE, ">", TEMP_COLLECTION_FILE, sep=" ")
system(pasteCmd) # Paste the files to a temporary location, just so that if it gets interrupted, the "final" file will not be present in malformed state
moveToFinalLocCmd <- paste("mv", TEMP_COLLECTION_FILE, "", FINAL_DE_FILE, sep=" ")
system(moveToFinalLocCmd) # Now "atomically" move the temporary file to the final file

print("[DONE]")
print("Next you should probably join this output file with some human-readable gene annotation file...")
print("Example for MOUSE: paste EDGER.OUT.*.txt > all.txt ; join.pl -o 'NO_ANNOT' all.txt /work/Common/Data/Annotation/mouse/2011_Archive/Mouse_Ensembl_mm9_Biomart_2012.txt  > all_annot_mm9.txt ")
print("Example for HUMAN: paste EDGER.OUT.*.txt > all.txt ; join.pl -o 'NO_ANNOT' all.txt /work/Common/Data/Annotation/human/2011_Archive/Human_Ensembl_hg19_Biomart_2011.txt > all_annot_hg19.txt ")
