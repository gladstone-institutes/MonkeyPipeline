#!/usr/bin/perl

# Max walltime is 336 hours (14 days)! Jobs requesting more are auto-rejected.

# ABOUT WALLTIME: Make sure "walltime" is high for running the RSEQC jobs, or else jobs will be terminated while they're still running. RSEQC is VERY SLOW!
# RSEQ POTENTIALLY TAKES A VERY LONG TIME TO RUN!
# BE AWARE OF THAT!

# Runs "RSEQC" quality control metrics.
# Rseqc requires "R"!
# RSEQC is basically like a more inconvenient version of FastQC---but it has a few DIFFERENT metrics!
# RNASeq QC program.
#      Docs: http://dldcc-web.brc.bcm.edu/lilab/liguow/CGI/rseqc/_build/html/#download-rseqc
# You'll also need these annotation files first:
# ## Human hg19:         wget http://dldcc-web.brc.bcm.edu/lilab/liguow/RSeQC/dat/hg19_Ensembl.bed.gz
# ## Mouse mm9:          wget http://dldcc-web.brc.bcm.edu/lilab/liguow/RSeQC/dat/mm9_NCBI37_Ensembl.bed.gz
# ##Mouse mm10           wget http://dldcc-web.brc.bcm.edu/lilab/liguow/RSeQC/dat/GRCm38_mm10_Ensembl.bed.gz
# #Fly (D. melanogaster) (BDGP R5/dm3)   wget http://dldcc-web.brc.bcm.edu/lilab/liguow/RSeQC/dat/fly_dm3_EnsemblGene.bed.gz

use strict; use warnings; use Carp; # Carp = "confess" function
use File::Spec::Functions qw(catfile);
use Data::Dumper;
use File::Basename;
use Digest::MD5; # Let's us make a "digest" (a hash function) of a long full-path string to an MD5 one
use lib "/data/work/Code/alexgw/monkey_agw"; # To make sure we find bananas_agw
use bananas_agw;

bananas_agw::requireEnvOrDie('outputDir','inputBAM','bedAnnot','python2_7_exe','scriptParentDir'); # <-- these vars MUST be defined in the %ENV hash

my $force              = bananas_agw::envLooksTrue("force");
my $verbose            = bananas_agw::envLooksTrue("verbose");
my $debug              = bananas_agw::envLooksTrue("debug"); # extra-verbose

#my $PERL_PBS_JOBNAME       = $ENV{'PBS_JOBNAME'};  ($PERL_PBS_JOBNAME ne '') or confess "Could not get the 'PBS_JOBNAME' field from the ENV.";
#my $PERL_PBS_JOBID         = $ENV{'PBS_JOBID'};    ($PERL_PBS_JOBID   ne '') or confess "Could not get the 'PBS_JOBID' field from the ENV.";
my $outd                   = $ENV{'outputDir'};
my $inputBAMFullPath       = $ENV{'inputBAM'};
my $refBedFile             = $ENV{'bedAnnot'}; # <-- you can download these from the RSEQC web site if you don't have them for whatever reason (as described in the error message below)
my $python2_7_exe          = $ENV{'python2_7_exe'};   # A version of python2 that is >= 2.7. 2.6 won't work.
my $scriptParentDir        = $ENV{'scriptParentDir'}; # The location of RSEQC python files
my $samtools               = $ENV{'samtools'};
my $shouldUseAllReads      = $ENV{'shouldUseAllReads'}; # 0 = just use a subset of reads. 1 = use ALL reads.
my $N_BAM_LINES_TO_EXAMINE = $ENV{'nBamLines'}; # 250000; # If we aren't '$shouldUseAllReads'-ing, we will only examine this many bam lines for the (slow) RSEQC operations

(-e $inputBAMFullPath) or confess "The BAM input file <$inputBAMFullPath> (which should probably be a FULL PATH to a BAM file) could not be found! Check to make sure that this file exists and is a full path (not a relative path).";
(-e $refBedFile)       or confess "The input reference BED file <$refBedFile> (which should probably be a FULL PATH to a BED file) could not be found! Double check that it exists. You can also download these reference BED files DIRECTLY from the RSEQC web site, at this URL: http://dldcc-web.brc.bcm.edu/lilab/liguow/CGI/rseqc/_build/html/#download-rseqc .";
(-x $python2_7_exe)    or confess "The input python 2.7 (or later) python executable could NOT be run (it was not executable)! Check this path to make sure it is valid: ${python2_7_exe}";
(-d $scriptParentDir)  or confess "The rseqc parent directory ($scriptParentDir) (where we expected to find files like 'geneBody_coverage.py' and the like) was not found! Check to make sure that this is a valid directory with rseqc python files in it: ${scriptParentDir}";

# if RSEQC is installed, we also want to know gene body coverage, junction annotation, etc etc...
# 'refBed' is a reference bed file downloaded from http://dldcc-web.brc.bcm.edu/lilab/liguow/CGI/rseqc/_build/html/#download-rseqc
#my $refBed="/data/info/genome/rseqc_bed_files/mm9_NCBI37_Ensembl.bed";
my $bamBase = basename($inputBAMFullPath);
$bamBase  =~ s/[.][sb]am$//i;     # <-- remove the bam/sam suffix
my $outpre = catfile($outd, $bamBase);

# ================== EXIT EARLY IF THE OUTPUT ALREADY EXISTS!!!! =================
my $successFileFullPath = "${outpre}_stats.txt"; # Just some output file that we expect to exist.
if (!$force and (-e "${successFileFullPath}" or -e "${successFileFullPath}.gz" or -e "${successFileFullPath}.bz2")) {
    ($verbose) && print STDERR "[OK] [Skipping RSEQC re-run] since the 'successfully run' sentinel file <$successFileFullPath> already exists.\n";
    exit(0);
}

bananas_agw::mkdirOrDie($outd);

my $inb       = 'ERROR_THIS_MUST_BE_SET_BELOW';
my $tmpHead   = `/bin/mktemp "${outd}/rseqc_tmp_header_XXXXXXX" --suffix ".sam"`; chomp($tmpHead);                                   ###catfile($outd, "tmp_header_${PERL_PBS_JOBNAME}.${PERL_PBS_JOBID}.out.sam");
my $tmpSubBam = `/bin/mktemp "${outd}/rseqc_tmp_only_${N_BAM_LINES_TO_EXAMINE}_reads_XXXXXXX" --suffix=".out.bam"`; chomp($tmpSubBam); #catfile($outd, "tmp_only_${N_BAM_LINES_TO_EXAMINE}_reads_${PERL_PBS_JOBNAME}.${PERL_PBS_JOBID}.out.bam");
my $tmpSubBai = $tmpSubBam . ".bai"; chomp($tmpSubBai); # the bam index file that is auto-generated

if ($shouldUseAllReads) {
	$inb = $inputBAMFullPath; # We are using ALL the reads, which will take a SUPER DUPER long time for a normal-sized file.
	# 6+ hours for the geneBody_coverage and RPKM_saturation parts.
	# Generally there is no reason to use ALL the reads. A subset should be sufficient.
} else {
	# Grab a SUBSET of reads!
	($N_BAM_LINES_TO_EXAMINE >= 1) or confess "The N_BAM_LINES_TO_EXAMINE variable must be at least 1! Really it should be at least 10,000... you set it to '$N_BAM_LINES_TO_EXAMINE' however, which is not valid!";
	bananas_agw::systemAndLog("${samtools} view -H ${inputBAMFullPath} > ${tmpHead}", $verbose);
	(-e $tmpHead) or confess "ERROR: Failed to generate the temporary header file '$tmpHead'...";
	bananas_agw::systemAndLog("${samtools} view ${inputBAMFullPath} | head -n ${N_BAM_LINES_TO_EXAMINE} | cat ${tmpHead} - | ${samtools} view -bS - > ${tmpSubBam}", $verbose);
	bananas_agw::systemAndLog("${samtools} index ${tmpSubBam}", $verbose); # '.bai' file is also mandatory. This will generate the file "$tmpSubBai", even though we don't explicitly specify that name here.
	(-e $tmpSubBam) or confess "ERROR: Failed to generate the subset bam file '$tmpSubBam'...";
	(-e $tmpSubBai) or confess "ERROR: we were not able to properly generate the subset .bai BAM index file <$tmpSubBai>...";
	$inb = $tmpSubBam;
}

# Note: some of these commands are QUITE SLOW! This can take 4+ hours per BAM file!!!
(-e "${inb}.bai") or print STDERR "Warning: geneBody_coverage.py is, for some reason, super finnicky and now requires a .bai file, which did not exist for the bam file ${inb}. The generation of that figure may fail.\n";
bananas_agw::systemAndLog("$python2_7_exe ${scriptParentDir}/inner_distance.py      -r ${refBedFile} -i ${inb} -o ${outpre}", $verbose); # <-- Makes PDF histogram: for paired end files, shows the distribution of inner distance between pairs
bananas_agw::systemAndLog("$python2_7_exe ${scriptParentDir}/junction_annotation.py -r ${refBedFile} -i ${inb} -o ${outpre}", $verbose); # <-- Makes 2 PDF pie charts: show whether novel/known junctions were found
bananas_agw::systemAndLog("$python2_7_exe ${scriptParentDir}/junction_saturation.py -r ${refBedFile} -i ${inb} -o ${outpre}", $verbose); # <-- Makes 1 PDF line plot: shows junctions that are discovered at various sequencing depths
bananas_agw::systemAndLog("$python2_7_exe ${scriptParentDir}/read_duplication.py                     -i ${inb} -o ${outpre}_duplication", $verbose); # Makes a PDF scatterplot of reads vs frequency: Not sure how to interpet it...
#bananas_agw::systemAndLog("$python2_7_exe ${scriptParentDir}/read_NVC.py               --nx      -i ${inb} -o ${outpre}_nvc", $verbose); # Nucleotide frequency. Redundant with FASTQC's (better) version of this same plot.
bananas_agw::systemAndLog("$python2_7_exe ${scriptParentDir}/infer_experiment.py    -r ${refBedFile} -i ${inb} &> ${outpre}_inferred_experiment.txt", $verbose); # <-- Generates .TXT file showing whether the bam file is paired end and/or has strand-specific reads
bananas_agw::systemAndLog("$python2_7_exe ${scriptParentDir}/read_distribution.py   -r ${refBedFile} -i ${inb} &> ${outpre}_read_dist.txt", $verbose); # <-- Generates .TXT file with number of reads found in Introns/UTR/Exons/etc

# bam_stats.py is VERY SLOW (4+ hours for a run)
bananas_agw::systemAndLog("$python2_7_exe ${scriptParentDir}/bam_stat.py                   --mapq=30 -i ${inb} &> ${outpre}_stats.txt", $verbose); # <-- Generates .TXT files with data like "Splice reads: 10736655 / Mapped in proper pairs: 56207918"

# RPKM_saturation.py is SUPER DUPER slow (6+ hours for a run)
bananas_agw::systemAndLog("$python2_7_exe ${scriptParentDir}/RPKM_saturation.py     -r ${refBedFile} -i ${inb} -o ${outpre}_saturation", $verbose); # <-- PDF: Some quartile plot that I don't understand. 4 box-and-whisker graphs in a square

# geneBody_coverage.py is SUPER DUPER slow (6+ hours for a run)
bananas_agw::systemAndLog("$python2_7_exe ${scriptParentDir}/geneBody_coverage.py   -r ${refBedFile} -i ${inb} -o ${outpre}", $verbose); # <-- Makes PDF histogram: shows read distribution from 5' to 3' across all genes in $refBedFile

# Bzip2 all the (not immediately useful to us) text files, R files, and XLS files so they don't take up a ton of space.
bananas_agw::systemAndLog("bzip2 ${outpre}*.txt ${outpre}*.[rR] ${outpre}*.xls", $verbose); # <-- GZIP the (potentially large!) outputs. In my testing, reduces the output from over 1 GB to ~300 MB.
#bananas_agw::systemAndLog(qq{echo "[DONE with RSEQC for $baseBam]" > $successFileFullPath }, $verbose); # <-- We are done! Do not re-run this (slow) set of operations!

# If the temporary files were created above, then delete them now.
(-e $tmpHead)   and bananas_agw::systemAndLog("/bin/rm $tmpHead", $verbose);   # Delete a temp header file...
(-e $tmpSubBam) and bananas_agw::systemAndLog("/bin/rm $tmpSubBam", $verbose); # Delete a temp "bam subset" file...
(-e $tmpSubBai) and bananas_agw::systemAndLog("/bin/rm $tmpSubBai", $verbose); # Delete a temp "bai subset" file...

