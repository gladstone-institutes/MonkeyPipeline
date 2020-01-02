#!/bin/perl





# Warning: if you have a single quote on a PBS comment line above (e.g., a word like "don't"), it will trigger a "qsub: unmatched" error and exit code 1.
my $syslog = "X.46.qc.multibam.syslog.txt";
# Runs "RSEQC" quality control metrics, but only the ones that run on ALL input bam files, not the ones that require single bam files only.
# See "47.rseq.qc.pl" for the SINGLE BAM FILE QC metrics.
# Rseqc requires "R" to be installed to run! Docs: http://dldcc-web.brc.bcm.edu/lilab/liguow/CGI/rseqc/_build/html/#download-rseqc

use strict; use warnings; use Carp; # Carp = "confess" function
use File::Spec::Functions qw(catfile);   use Data::Dumper;   use File::Basename;
use lib "/wynton/group/gladstone/biocore/MonkeyPipeline"; use bananas_agw; # To make sure we find bananas_agw

bananas_agw::requireEnvOrDie('inputBamDir','outputDir','scriptParentDir','bedAnnot','python2_7_exe'); # <-- these vars MUST be defined in the %ENV hash
my $force              = bananas_agw::envLooksTrue("force");
my $verbose            = bananas_agw::envLooksTrue("verbose");
my $debug              = bananas_agw::envLooksTrue("debug");
my $inputBamDir        = $ENV{'inputBamDir'};
my $scriptParentDir    = $ENV{'scriptParentDir'}; # The location of RSEQC python files
my $refBedFile         = $ENV{'bedAnnot'}; # <-- you can download thes from the RSEQC web site if you don't have them for whatever reason (as described in the error message below)
my $python2_7_exe      = $ENV{'python2_7_exe'};     # A version of python that

(-d $inputBamDir)      or confess "The BAM input directory <$inputBamDir> (which should probably be a FULL PATH to a directory full of BAM files) could not be found! Check to make sure that this directory exists and is a full path (not a relative path).";
(-e $refBedFile)       or confess "The input reference BED file <$refBedFile> (which should probably be a FULL PATH to a BED file) could not be found! Double check that it exists. You can also download these reference BED files DIRECTLY from the RSEQC web site, at this URL: http://dldcc-web.brc.bcm.edu/lilab/liguow/CGI/rseqc/_build/html/#download-rseqc .";
(-x $python2_7_exe)    or confess "The input python 2.7 (or later) python executable could NOT be run (it was not executable)! Check this path to make sure it is valid: ${python2_7_exe}";
(-d $scriptParentDir)  or confess "The rseqc parent directory ($scriptParentDir) (where we expected to find files like 'geneBody_coverage.py' and the like) was not found! Check to make sure that this is a valid directory with rseqc python files in it: ${scriptParentDir}";

my $outPrefix = catfile($ENV{'outputDir'}, "0_all_bam_files"); # Gbc = genebody coverage

# DO NOT MANUALLY EDIT THE $outPdf name below! This is a filename that is auto-generated by RSEQC and we can't specify it (we could rename it, though).
my $outPdf    = "${outPrefix}.geneBodyCoverage.curves.pdf"; # <-- the ".geneBodyCoverage.curves.pdf" part is auto-added by the RSEQC script and is not specifyable by us

# ================== EXIT EARLY IF THE OUTPUT ALREADY EXISTS!!!!
if (!$force and -e $outPdf) { # Output file exists, AND we aren't forcing a re-run!
    ($verbose) && print STDERR "[OK] [Skipping RSEQC multi-bam re-run] since the output file <$outPdf> already exists.\n";
    exit(0);
}

bananas_agw::mkdirOrDie($ENV{'outputDir'});     # Guess we need to make the output files!
# Note that geneBodyCoverage is super picky now and requires .bai index files as well as .bam files! Annoying!

# Possible options for '-i': directory containing one or more bam files. Note that it appears to traverse subdirectories!
my $cmd = "$python2_7_exe ${scriptParentDir}/geneBody_coverage.py   -r ${refBedFile} -i ${inputBamDir} -o $outPrefix";
my $exitCode = bananas_agw::systemAndLog($cmd, $verbose, $syslog); # <-- Makes PDF histogram: shows read distribution from 5' to 3' across all genes in $refBedFile
if ($exitCode != 0) {
	print STDERR "[WARNING] Exit code of the geneBody_coverage.py command was not equal to zero!\n";
	print STDERR "               Here is the command that may have failed: $cmd\n";
} else {
	($verbose) && print STDERR "[OK] Looks like the all-samples gene body coverage finished and exited with code 0. Hopefully it was successful!\n";
}
