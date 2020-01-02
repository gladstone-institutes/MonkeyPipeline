#!/bin/perl





use strict; use warnings; use Carp; # Carp = "confess" function
use File::Spec::Functions qw(catfile);
#use Data::Dumper;
use lib "/wynton/group/gladstone/biocore/MonkeyPipeline"; use bananas_agw; # To make sure we find bananas_agw

bananas_agw::requireEnvOrDie('R_EXE','pairScript','subreadDir','subreadCountsFilename','outputDir'); # <-- these variables MUST be defined in the %ENV hash
my $force   = bananas_agw::envLooksTrue("force");
my $verbose = bananas_agw::envLooksTrue("verbose");
my $debug   = bananas_agw::envLooksTrue("debug"); # extra-verbose

#my $SUBREAD_COUNT_FILE_BASENAME = "subread.featureCounts.counts.txt";
my $OUTPUT_SUBREAD_COUNT_PDF    = "subread.norm.cpm.pairs.pdf"; # <-- this is the expected output file name! This is important, since we want to detect whether it was made properly or not

my $R_EXE        = $ENV{'R_EXE'};
my $pairScript   = $ENV{'pairScript'};
my $subreadDir   = $ENV{'subreadDir'};
my $outDir       = $ENV{'outputDir'};
($ENV{'subreadCountsFilename'} =~ m/[-_.A-Za-z.]/) or confess "[ERROR]: The output subread filename contained some weird characters besides A-Z, 0-9, ., _, and -. This is NOT ALLOWED. The offending filename was: <$ENV{'subreadCountsFilename'}>.";

my $subreadCountsFileFullPath  = catfile(${outDir}, $ENV{'subreadCountsFilename'}); # Single and paired together. DO NOT CHANGE THIS FILENAME without changing the R scripts in "20b.scripts"--they use HARD CODED filenames! Sorry about that.

(-d $subreadDir)                or confess "[ERROR]: The input directory <$subreadDir> needs to exist. We could not find it!";
(-e $subreadCountsFileFullPath) or confess "[ERROR]: Could not find the expected mandatory input subread counts file <$subreadCountsFileFullPath>.";
(-x $R_EXE)                     or confess "[ERROR]: Could not find the specified R executable in <$R_EXE> (or possibly it could not be executed by the current user)! Double check that it exists and has the right permissions.";
(-e $pairScript)                or confess "[ERROR]: Could not find the specified 'pair plot' R script in <$pairScript>! Double check that it exists.";

my $outPdfFullPath = catfile($outDir, $OUTPUT_SUBREAD_COUNT_PDF);
if (!$force and (-e $outPdfFullPath)) { # Output files exist already, AND we aren't forcing a re-run!
    ($verbose) && print STDERR "[OK -- SKIPPING] Skipping the re-run of the Pairs Plot, because the output PDF plot named <$outPdfFullPath> exists already.\n";
    exit(0);
}

bananas_agw::mkdirOrDie($outDir);
bananas_agw::systemAndLog("cd $outDir && ${R_EXE} --vanilla < $pairScript", $verbose); # 'cd' to Run it from the output directory! It makes the output file from wherever the script was run from.


