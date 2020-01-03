#!/bin/perl
eval 'exec /bin/perl -w -S $0 ${1+"$@"}'
if 0; # not running under some shell

use strict; use warnings; use Carp; # Carp = "confess" function
use File::Spec::Functions qw(catfile);
#use Data::Dumper;
use lib "/wynton/group/gladstone/biocore/MonkeyPipeline"; use bananas_agw; # To make sure we find bananas_agw

bananas_agw::requireEnvOrDie('R_EXE','script','cuffdiffDir','outputDir'); # <-- these variables MUST be defined in the %ENV hash
my $force   = bananas_agw::envLooksTrue("force");
my $verbose = bananas_agw::envLooksTrue("verbose");
my $debug   = bananas_agw::envLooksTrue("debug"); # extra-verbose

my $R_EXE        = $ENV{'R_EXE'};
my $script       = $ENV{'script'};
my $cuffdiffDir  = $ENV{'cuffdiffDir'};
my $outd         = $ENV{'outputDir'};

my $OUTPUT_SUCCESS_FILE_FULL_PATH = catfile($outd, "z.success.cummerbund"); # Do NOT change this without changing the corresponding call in "23b.cummerbund.R". This file will be generated at the end of "23b.cummerbund.R", assuming it runs successfully!
my $cuffIn                       = catfile($cuffdiffDir, "gene_exp.diff"); # Just a file we expect to exist as a prerequisite, assuming that Cuffdiff was run properly.

(-e $R_EXE)           or confess "[ERROR]: Could not find the specified R executable in <$R_EXE>! Double check that it exists.";
(-e $script)          or confess "[ERROR]: Could not find the specified 'cummerbund plot' R script in <$script>! Double check that it exists.";
(-d $cuffdiffDir)     or confess "[ERROR]: The input directory <$cuffdiffDir> needs to exist. We could not find it!";
(-r $cuffIn)          or confess "[ERROR]: Could not find certain expected mandatory input file in the $cuffdiffDir!";

if (!$force and -e $OUTPUT_SUCCESS_FILE_FULL_PATH) { # Output files exist already, AND we aren't forcing a re-run!
    ($verbose) && print STDERR "[OK] [SKIPPING RE-RUN]: Skipping re-run of the cummeRbund cuffdiff graphical plot generation, because the output file ($OUTPUT_SUCCESS_FILE_FULL_PATH) already exists.\n";
} else {
    bananas_agw::mkdirOrDie($outd);
    bananas_agw::systemAndLog("cd $outd && ln -f -s ${cuffdiffDir}/* ./", $verbose); # Link all the cuffdiff input files in this folder
    bananas_agw::systemAndLog("cd $outd && ${R_EXE} --vanilla < $script", $verbose); # 'cd' to Run it from the output directory! It makes the output file from wherever the script was run from.
}

