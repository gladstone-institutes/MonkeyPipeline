#!/usr/bin/perl




# Warning: if you have a single quote on a PBS comment line above (e.g., a word like "don't"), it will trigger a "qsub: unmatched" error and exit code 1.

use strict; use warnings;
use File::Basename;
use File::Spec::Functions qw(catfile);
use Data::Dumper;
use Carp; # backtrace on errors. Has the "confess" function. Use this instead of "die" if you want useful information!
use lib "/data/work/Code/alexgw/monkey_agw"; # To make sure we find bananas_agw
use bananas_agw;

bananas_agw::requireEnvOrDie('RSCRIPT_EXE', 'edgeRScript', 'subreadCountsFilePath', 'outputDir', 'groupList', 'bamsPerGroup', 'majorDelim', 'minorDelim'); # <-- these vars MUST be defined in the %ENV hash

my ($force, $verbose, $dbug) = bananas_agw::envLooksTrue("force", "verbose", "debug");
my $RSCRIPT_EXE     = $ENV{'RSCRIPT_EXE'};             (-x $RSCRIPT_EXE)    or confess "[ERROR]: Could not find/execute the R executable in <$RSCRIPT_EXE>.";
my $edgeRscript   = $ENV{'edgeRScript'};       (-e $edgeRscript)  or confess "[ERROR]: Could not find mandatory edgeR script file <$edgeRscript>.";
my $inCounts  = $ENV{'subreadCountsFilePath'}; (-e $inCounts) or confess "[ERROR]: Could not find mandatory 'subread featureCounts' input counts file at the following *FULL PATH*: <$inCounts>.";
my $glist     = $ENV{'groupList'}; #(-e $inCounts) or confess "[ERROR]: Could not find mandatory 'subread featureCounts' input counts file at the following *FULL PATH*: <$inCounts>.";
my $bams      = $ENV{'bamsPerGroup'}; #(-e $inCounts) or confess "[ERROR]: Could not find mandatory 'subread featureCounts' input counts file at the following *FULL PATH*: <$inCounts>.";
my $major     = $ENV{'majorDelim'};
my $minor     = $ENV{'minorDelim'};
my $outputDir = $ENV{'outputDir'};

my $finalDiffExpressionFile = catfile($outputDir, "EDGER.all.table.txt"); # <-- Note: HARD-CODED "EDGER.all.table.txt" -- If you change this, also change it in 50b.edgeR.R. Do NOT change this one without also changing the "FINAL_DE_FILE" variable in 50b.edgeR.pl!

if (!$force and -e $finalDiffExpressionFile) {
    ($verbose) && print STDERR "[OK -- SKIPPING]: Looks like the output edgeR differnetial expression file '$finalDiffExpressionFile' already exists, so we are not re-running edgeR.\n";
    exit(0);
}

bananas_agw::mkdirOrDie($outputDir);
bananas_agw::systemAndLog(qq{cd "$outputDir" } # 'cd' to run it from the output directory! It makes the output file from wherever the script was run from.
			  . qq{ && export SUBREAD_COUNTS_FILE="$inCounts" } # <-- R reads 'export'-ed environment variables in (like command line arguments)
			  . qq{                    GROUP_LIST="$glist" }       # <-- R reads 'export'-ed environment variables in (like command line arguments)
			  . qq{                   MAJOR_DELIM="$major" }
			  . qq{                   MINOR_DELIM="$minor" }
			  . qq{             BAM_LIST_BY_GROUP="$bams" } # <-- R reads 'export'-ed environment variables in (like command line arguments)
#			  . qq{ &&  ${RSCRIPT_EXE} --vanilla < "$edgeRscript" }, $verbose);
			  . qq{ &&  ${RSCRIPT_EXE} --vanilla "$edgeRscript" }, $verbose);

#bananas_agw::systemAndLog("chmod a+r $trackInfoLocalPath $localOutDirPath $localOutDirPath/Browser*", $verbose); # Make sure everyone can READ the browser files


