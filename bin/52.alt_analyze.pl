#!/bin/perl -w

use strict; use warnings; use Carp; # Carp = "confess" function
use File::Spec::Functions qw(catfile);   use Data::Dumper;   use File::Basename;
use lib "/wynton/group/gladstone/biocore/MonkeyPipeline"; use bananas_agw; # To make sure we find bananas_agw

# Note about AltAnalyze: all paths MUST be full paths at all times, relative paths never work.

# REQUIRES the "numpy" package to be installed
# REQUIRES that the "requests" package be installed
#      sudo pip install --upgrade requests
#      sudo pip install --upgrade numpy
#      sudo pip install --upgrade matplotlib
#      pip  install --ignore-installed --install-option="--prefix=/data/applications/2015_06/libpython2.7/"

# How to test it:

# First you definitely need the required databases:
# python `which AltAnalyze.py` --species Mm --update Official --version EnsMart65 --additional all ; python `which AltAnalyze.py` --species Hs --update Official --version EnsMart65 --additional all ; python `which AltAnalyze.py` --species Rn --update Official --version EnsMart77 --additional all

# 1. mkdir ~/a
# 2. cd ~/a
# 3. monkey_agw --test1
# 4. cat ~/a/20.subread_counts/subread.featureCounts.counts.txt |  featureCounts_to_plain_matrix.sh > simple.txt
# 5. mkdir -p ~/a/ALTTEST/
# 6. python2.7 /data/applications/bin/AltAnalyze.py --input ~/a/20.subread_counts/simple.txt --runLineageProfiler yes --platform RNASeq --species Hs --output ~/a/20.subread_counts/ALTTEST/

# Note that ALTANALYZE actually likes to run in python2.6 for some reason. So maybe:
# grep -r python2.6 /path/to/altanalyze/

# Make sure those "shebang" invocations have been changed to "/usr/bin/env python2"

# Optional, but recommended:
#      sudo yum install tkinter python-imaging

bananas_agw::requireEnvOrDie('python2_exe', 'alt_analyze_py', 'featureCountsToPlainMatrix', 'subreadCountsFilename', 'species', 'outputDir'); # <-- MUST already be defined in the %ENV hash

sub speciesToAltAnalyzeFormat($) {
	# Input: species name in some (string) format). Output: the AltAnalyze-style format (like 'Hs' for human)
	my ($in) = @_;
	if ($in =~ m/^(human|hg\d|Hs)/i or $in =~ m/(\b|[._-])(human|hg19|hg38)([\b|[._-])/i) { return "Hs"; } # Hs = human (homo sapiens)
	if ($in =~ m/^(mouse|mm\d|Mm)/i or $in =~ m/(\b|[._-])(mouse|mm9|mm10)([\b|[._-])/i)  { return "Mm"; } # Mm = mouse (mus musculus)
	if ($in =~ m/^(rat)/i or $in =~ m/(\b|[._-])(rat)([\b|[._-])/i)  { return "Rn"; } # Rn = rat
	return undef; # huh, guess it wasn't recognized
}

my ($force, $verbose, $dbug) = bananas_agw::envLooksTrue("force", "verbose", "debug");
my $python2_exe                = $ENV{python2_exe};
my $altAnalyzePy               = $ENV{alt_analyze_py};
my $featureCountsToPlainMatrix = $ENV{featureCountsToPlainMatrix}; # <-- this is a really basic script that just turns a 'subread featureCounts' matrix into a plain data table. Normally it is named "featureCounts_to_plain_matrix.sh"
my $subreadInputMatrix         = $ENV{subreadCountsFullPath};
my $speciesAltAnalyzeFormat    = speciesToAltAnalyzeFormat($ENV{species});

my $outDir                     = $ENV{outputDir};
my $altOutFile                 = catfile($outDir, "altanalyze.run.touch");

if (!$force and -e $altOutFile) {
    ($verbose) && print STDERR "[OK -- SKIPPING]: Looks like the output AltAnalyze data already existed, so we are not re-running it.\n";
    exit(0);
}


if (!defined($speciesAltAnalyzeFormat)) {
	warn("[WARNING]: Skipping AltAnalyze lineage profiler analysis, because the input species (\"$ENV{species}\") seemed not to be human or mouse. We do not support LineageProfiler on any other species at the moment.");
} else {

	bananas_agw::dieIfFileAccessFails($python2_exe);
	bananas_agw::dieIfFileAccessFails($altAnalyzePy);
	bananas_agw::dieIfFileAccessFails($featureCountsToPlainMatrix);
	bananas_agw::dieIfFileAccessFails($subreadInputMatrix);

	bananas_agw::dieIfMissingPythonLib($python2_exe, "matplotlib"); # altanalyze REQUIRES matplotlib, and fails later if it is missing.
	bananas_agw::dieIfMissingPythonLib($python2_exe, "requests"); # altanalyze REQUIRES matplotlib, and fails later if it is missing.
	bananas_agw::dieIfMissingPythonLib($python2_exe, "Tkinter"); # altanalyze REQUIRES matplotlib, and fails later if it is missing.

	bananas_agw::mkdirOrDie($outDir);
	my $simplifiedMatrix  = catfile($outDir, "simplified.subread.featureCounts.matrix.txt");
	my $simplifyCmd       = qq{$featureCountsToPlainMatrix $subreadInputMatrix --base > $simplifiedMatrix}; # run a command that baiscally just cuts off the first few columns (removing them)
	bananas_agw::systemAndLog($simplifyCmd, $verbose, undef, "DIE IF NONZERO");
	(-e $simplifiedMatrix) or confess "[ERROR]: Somehow failed to generate the simplified matrix '$simplifiedMatrix'. Probably was a problem with the command '$simplifyCmd'---most likely, the script '$featureCountsToPlainMatrix' had some issue.";
	my $altCmd = qq{$python2_exe $altAnalyzePy }
	  . qq{ --input              $simplifiedMatrix }
	  . qq{ --runLineageProfiler yes }
	  . qq{ --platform           RNASeq }
	  . qq{ --species            $speciesAltAnalyzeFormat }
	  . qq{ --output             "$outDir" };

	print STDERR "GENERAL NOTE: If you need to update the databases (mouse), you can run:   python2 AltAnalyze.py --update Official --species Mm --version EnsMart65 --additional all  \n";
	print STDERR "GENERAL NOTE: If you need to update the databases (human), you can run:   python2 AltAnalyze.py --update Official --species Hs --version EnsMart65 --additional all  \n";
	print STDERR "GENERAL NOTE: Note: if you got the message 'Please note: LineageProfiler not currently supported for this species...' above, that is (probably) not actually true! It really means you need to download the databases. Supported species are: Hs (human), Mm (mouse), Rn (rat). Here are some update commands you could try running to get the required databases:\n * python `which AltAnalyze.py` --species Mm --update Official --version EnsMart65 --additional all\n * python `which AltAnalyze.py` --species Hs --update Official --version EnsMart65 --additional all\n * python `which AltAnalyze.py` --species Rn --update Official --version EnsMart77 --additional all\n";
	bananas_agw::systemAndLog($altCmd, $verbose, undef, "DIE IF NONZERO");
	# If we could capture the STDERR output, we could look for this to detect missing databases
	# >>> Please note: LineageProfiler not currently supported for this species...

	# AltAnalyze REQUIRES full paths to all files at all times!
	# python2.7 /data/applications/bin/AltAnalyze.py --input /data/home/alexgw/a/20.subread_counts/simple.txt --runLineageProfiler yes --platform RNASeq --species Hs --output "/data/home/alexgw/a/20.subread_counts/ALTTEST"
}

bananas_agw::systemAndLog(qq{touch "$altOutFile"}); # indicate that we finished successfully. Even do this if we skipped the analysis due to it being a non-applicable species, since it was still a "success"

# From Nathan:
#One of the most useful outputs that Alex Zambon pushed me to get integrated is the heatmap output with the enriched pathways shown on the left. To run this from the command-line:
#
#python AltAnalyze.py --image hierarchical --platform RNASeq --species Hs --display True --color_gradient red_black_sky  --normalization median --column_method hopach --row_method hopach --column_metric cosine --row_metric correlation --input /Users/saljh8/Desktop/Stanford/EB_d17-SingleCell.txt --clusterGOElite BioMarkers
#
#  or --clusterGOElite GeneOntology
#
#  or --clusterGOElite WikiPathways
#
#  or --clusterGOElite KEGG
#
#  or --clusterGOElite TFTarget
#
#  or --clusterGOElite PathwayCommons
#
#  or --clusterGOElite CTDOntology
#

#OUTDIR=~/Desktop/ALT_OUT ; mkdir -p $OUTDIR ; python2.7 AltAnalyze.py --output="$OUTDIR" --input "/Users/alexgw/Desktop/ELIJ/EM_460_COUNTS.txt" --runLineageProfiler yes --vendor RNASeq --platform RNASeq  --species Mm
  
#bananas_agw::systemAndLog("chmod a+r $trackInfoLocalPath $localOutDirPath $localOutDirPath/Browser*", $verbose); # Make sure everyone can READ the browser files


