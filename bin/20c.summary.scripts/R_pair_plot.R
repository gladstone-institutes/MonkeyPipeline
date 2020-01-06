require(edgeR); require(gplots)
source("/wynton/home/pollard/wmaguire/gladstone/biocore/references/agwUtil.R"); # Has "pairs.agw" in it

SUBREAD_FILE <- "subread.featureCounts.counts.txt"
OUT_PREFIX      <- "subread.norm.cpm.pairs"
# Looks for a file IN THE SAME DIRECTORY named "subread.featureCounts.counts.txt" (SUBREAD_FILE)--a "subread featureCounts"-format file
# and makes a "Pairs" plot PDF file named "subread.norm.cpm.pairs.pdf" (OUT_PDF)

ddd <- read.table(SUBREAD_FILE, sep="\t", header=TRUE, check.names=F, stringsAsFactors=FALSE, fill=FALSE, row.names=1, quote='', comment.char='', skip=1) # Skip the first line, with '#'
ddd[ddd == '']  <- 0
ddd[is.na(ddd)] <- 0
gene.len.vec <- ddd[["Length"]]
ds <- ddd[,6:ncol(ddd)]
dcomplete <- data.matrix(ds) # only the ACTUAL DATA, not the various header rows
#dc.table  <- data.table::data.table(ids=rownames(dcomplete), length=gene.len.vec, counts=dcomplete)

COUNT_PER_MILLION_MIN_CUTOFF <- 0.5     # Very conservative! EdgeR recommends 100.
COUNT_PER_MILLION_MAX_CUTOFF <- 5000  # Should filter out mitochondrial genes only at this point
keepLo <- rowSums(edgeR::cpm(dcomplete) >= COUNT_PER_MILLION_MIN_CUTOFF) >= 2 # At least TWO samples had a count-per-million greater than the cutoff
keepHi <- rowSums(edgeR::cpm(dcomplete) <= COUNT_PER_MILLION_MAX_CUTOFF) >= 2 # At least TWO samples were under the "too many counts" threshold.
keep <- keepLo & keepHi
print(paste("Keeping a total of ", sum(keep), " genes out of a population of ", length(keep), ". ", sum(!keepLo), " didn't have enough expression, and ", sum(!keepHi), " had too much.", sep=''))
dfilt <- dcomplete
dfilt <- dcomplete[keep, ] # Delete rows with too little/much expression

y   <- edgeR::DGEList(counts=dfilt, group=colnames(dfilt))
cpm <- edgeR::cpm(y, normalized.lib.sizes=TRUE, log=FALSE)

mmm <- log2(1+cpm)
colnames(mmm) <- gsub(".*[/]", "", colnames(mmm), perl=T)
colnames(mmm) <- gsub("[.](sam|bam|sam.gz)$", "", colnames(mmm), perl=T)
colnames(mmm) <- gsub("[_.](mm9|mm10|hg18|hg19).*", "", colnames(mmm), perl=T)

# Make a PDF
pdf(paste(OUT_PREFIX, ".pdf", sep=''), width=10+0.75*ncol(mmm), height=10+0.75*ncol(mmm)) # Increase the plot size by 3/4 inch per sample. Start at 10 inches.
par(pty='s')
pairs.agw(mmm, main=paste("log2(Normalized CPM). 'r' is Pearson's R. Points are excluded if CPM in all samples is over ", COUNT_PER_MILLION_MAX_CUTOFF, " or under ", COUNT_PER_MILLION_MIN_CUTOFF, sep=''))
dev.off()

# Make a PNG image. Not sure this is actually a good resolution though!
png(paste(OUT_PREFIX, "-144dpi.png", sep=''), res=144, width=2000+125*ncol(mmm), height=2000+125*ncol(mmm)) # Increase the plot size by 3/4 inch per sample. Start at 10 inches.
par(pty='s')
pairs.agw(mmm, main=paste("log2(Normalized CPM). 'r' is Pearson's R. Points are excluded if CPM in all samples is over ", COUNT_PER_MILLION_MAX_CUTOFF, " or under ", COUNT_PER_MILLION_MIN_CUTOFF, sep=''))
dev.off()

# Make a PNG image. Not sure this is actually a good resolution though!
png(paste(OUT_PREFIX, "-72dpi.png", sep=''), res=72, width=2000+125*ncol(mmm), height=2000+125*ncol(mmm)) # Increase the plot size by 3/4 inch per sample. Start at 10 inches.
par(pty='s')
pairs.agw(mmm, main=paste("log2(Normalized CPM). 'r' is Pearson's R. Points are excluded if CPM in all samples is over ", COUNT_PER_MILLION_MAX_CUTOFF, " or under ", COUNT_PER_MILLION_MIN_CUTOFF, sep=''))
dev.off()

# Make a PNG image. Not sure this is actually a good resolution though!
png(paste(OUT_PREFIX, "-300dpi.png", sep=''), res=300, width=2000+125*ncol(mmm), height=2000+125*ncol(mmm)) # Increase the plot size by 3/4 inch per sample. Start at 10 inches.
par(pty='s')
pairs.agw(mmm, main=paste("log2(Normalized CPM). 'r' is Pearson's R. Points are excluded if CPM in all samples is over ", COUNT_PER_MILLION_MAX_CUTOFF, " or under ", COUNT_PER_MILLION_MIN_CUTOFF, sep=''))
dev.off()
