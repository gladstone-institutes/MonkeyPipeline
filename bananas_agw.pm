#!/bin/perl
eval 'exec /bin/perl -w -S $0 ${1+"$@"}'
if 0; # not running under some shell

# bananas.pm version 2 - parallelized for rigel
# Authors: Sean Thomas and Alex Williams, 2013
# This program is designed to automate many of the tasks that accompany next generation sequencing tasks. 

# MERGER of CHIP/EXO and RNA pipelines completed Dec 15, 2014

# Alex's changed package name so as to not inadvertently interfere with Sean's.
# Thing to note:
#        * when jobs are queued up, they are added IN ALPHABETICAL ORDER
#        * a job that begins alphabetically earlier therefore CANNOT depend on a later job---example: job "08b" can depend on "08a" but not vice versa.
#        * thus, jobs must be numbered in such a way that alphabetization is related to prerequisites.


# To do:
#      * All of this assumes that you're running the executable from the same filesystem as the cluster.
#               So when we check for binaries, we are checking on the CLIENT (head/login node),
#               but what we REALLY need to do is submit a 'check dependencies' cluster job to check everything
#               on the cluster machine.     -- AlexW, July 18 2016

package bananas_agw;

use strict; use warnings; use Carp; # Carp has "confess"
use Data::Dumper;
use File::Basename;
use Scalar::Util qw(looks_like_number);
use File::Spec::Functions qw(catfile); # File::Spec has "catfile"
use IO::Compress::Gzip qw(gzip $GzipError) ;
use List::Util qw(max);
use Exporter;
use Sys::Hostname; # for the 'hostname()' command used in getUserPriv
use POSIX; # for 'getgroups()' in getUserPriv

### 0.i Export functions
our (@ISA, @EXPORT, @EXPORT_OK);
@ISA       = qw(Exporter);
@EXPORT_OK = qw(GLOBAL_VERBOSE GLOBAL_BIN_DIR GLOBAL_DEBUG GLOBAL_DRY_RUN GLOBAL_COLOR FORCE_RERUN OUTPUT_DIR_OVERRIDE GLOBAL_WARNINGS RUN_DIRECTLY_WITHOUT_TORQUE NA_VALUE);
@EXPORT    = qw(loadConfig checkConfig buildJobSubmissionScript runJobSubmissionScript noExtensions printInstallInstructions isNA);

(defined($ENV{BINFBINROOT}) and length($ENV{BINFBINROOT}) >= 1) or confess "-----------\n[ERROR]: You need to define the UNIX system-wide environment variable 'BINFBINROOT' (for 'bioinformatics binary root') in your ~/.bashrc file before you can run Monkey.\nHere is an example of something you can add to the ~/.bashrc:  export BINFBINROOT=/data/applications/2015_06/bin \n(Make sure to log out and log back in after you set that, or type 'source ~/.bashrc' to reload the bash config file). If you add that and it doesn't work, try checking the value of 'echo \$BINFBINROOT' and make sure that actually prints the proper path.\n------------ ";

our ($GLOBAL_BIN_DIR) = $ENV{BINFBINROOT}; # "/data/bin/2015_06"; # <== note: hard-coded, and NOT THE "monkey" /bin dir with the perl scripts in it! (That's the 'monkeyPoo' directory)

if (not -d $ENV{BINFBINROOT}) {
	warn("[WARNING]: Your environment variable 'BINFBINROOT' (for 'bioinformatics binary root') was not set to a valid directory that already exists! You will need to properly set the 'BINFBINROOT' in your ~/.bashrc file. You can check to see what this variable is set to by typing 'echo \$BINFBINROOT' on the command line. It should be a directory with a bunch of programs in it. The current (probably incorrect) one is set to this value: <$GLOBAL_BIN_DIR> (which does not seem to be a valid directory)");
	confess "Failure to find a valid directory with binaries in it, which we expect to be specified in the 'BINFBINROOT' environment variable.\nSee text above for details.\nIn other words---your environment variables in your ~/.bashrc are NOT set up correctly.";
}

## These are settings SPECIFIC to our server and queue system, NOT general settings!!!!! #####
my $ADMIN_GROUP_NAME  = "admin";
my $NORMAL_GROUP_NAME = "normal";

my %QSETTINGS = (
#			"dest"                  => {"$ADMIN_GROUP_NAME"   => "-q Bio"
#					     , "$NORMAL_GROUP_NAME" => "-q General" }
#		  , "grouplist"           => {"$ADMIN_GROUP_NAME"   => "-A bioqueue"
#					     , "$NORMAL_GROUP_NAME" => "-A genqueue" }
		    "ncpus_for_alignment" => {"$ADMIN_GROUP_NAME"   => "6"  # Privileged users = 4 threads per job, other users = 1 thread
					     , "$NORMAL_GROUP_NAME" => "2" }
		  , "tophat_needs"        => {"pbs_mem" => "12g" }
		  , "star_aligner_needs"  => {"pbs_mem" => "33g" } # needs a lot of ram! Default is 31 GB limit for just the aligner
		  , "gem"                 => {"pbs_mem" => "25g", "pbs_ncpus"=>1 } # needs a lot of ram! Default is 31 GB limit for just the aligner
		  , "default_ncpus"       => {"$ADMIN_GROUP_NAME"    => "1"
					      , "$NORMAL_GROUP_NAME" => "1" }
		  , "default_mem"         => {"$ADMIN_GROUP_NAME"    => "4g"
					    , "$NORMAL_GROUP_NAME"   => "4g" }
		  , "default_walltime"    => {"$ADMIN_GROUP_NAME"    => "259200"
					     , "$NORMAL_GROUP_NAME"  => "259200" }
		  );
my $UNIX_GRP_WORK       = "1004"; # <== Anybody in this group gets the be 'privileged' for job submission
my $UNIX_GRP_CAVE       = "3050"; # <== Anybody in this group gets the be 'privileged' for job submission
my $UNIX_GRP_B2B        = "2050"; # <== Anybody in this group gets the be 'privileged' for job submission
my $UNIX_GRP_BIOQUEUE   = "35098"; # <== Anybody in this group gets the be 'privileged' for job submission
my @UNIX_BIO_GRP_ID_ARRAY = ($UNIX_GRP_WORK, $UNIX_GRP_CAVE, $UNIX_GRP_B2B, $UNIX_GRP_BIOQUEUE); # HARD CODED right now for rigel specifically. Check /etc/groups for a list by ID.

my $SHOULD_ATTEMPT_TO_RUN_MONOCLE = 0; # Never did really get this working right...

## ABove: these are settings SPECIFIC to our server and queue system, NOT general settings!!!!! #####

our ($GLOBAL_VERBOSE)      = 0; # 0 = quiet, 1 = verbose
our ($GLOBAL_DEBUG)        = 0; # 0 = normal, 1 = ultra-verbose diagnostic messages
our ($GLOBAL_DEV_MODE)     = 0; # 0 = normal, 1 = ultra-verbose diagnostic messages
our ($GLOBAL_DRY_RUN)      = 0; # 0 = normal, 1 = dry run (dry = "don't actually do anything")
our ($GLOBAL_COLOR)        = 0; # 0 = only black and white messages. 1 = COLOR text to console
our ($FORCE_RERUN)         = 0; # 0 = normal, 1 = re-run all jobs, even completed ones
our ($OUTPUT_DIR_OVERRIDE) = undef; # If defined, this means the user has specified their own resultDir on the command line, rather than in the config file.
our ($GLOBAL_WARNINGS)     = ""; # Accumulate all the WARNING messages that we might want to tell the user both 1) when they occur and 2) at the end of the job submission
our ($RUN_DIRECTLY_WITHOUT_TORQUE) = 0; # By default, we will run with the torque job scheduler. But if this is specified, we will instead just run ONE job at a time, very slowly. Works even if Torque is not installed
our ($NA_VALUE)            = "NA"; # Don't change this unless you ALSO change the 'NA' value that is used to specifiy *intentionally* missing files in the input files.

my $poo_data_subdir = "data"; # subdirectory in "poo" that has certain data files

my $ZLM_LICENSE_HARD_CODED = "/home/sthomas/tools/vmatch-2.2.4-Linux_x86_64-64bit/vmatch.lic"; # <== warning: do not change the variable name: "ZLM_LICENSE" is some ENV variable that is needed by peakMotifs. See 10.motifDiscovery.pl for the exact command

# For the 'gem' peaks
my $GEM_READ_DIST_DEFAULT  = "Read_Distribution_default.txt"; # This file should be in the "monkeyPoo" directory
my $GEM_READ_DIST_CHIP_EXO = "Read_Distribution_ChIP-exo.txt"; # This file should be in the "monkeyPoo" directory

my $HOLD_UNTIL_ALL_JOBS_SUBMITTED = 1; # <== Should be 1. This is VERY IMPORTANT for jobs with large numbers of files (--test6 for an example). Set it to 0 to remove Alex's workaround for PBS failing when jobs depend on other already-completed jobs.
my $SIGNAL_TO_START_VARNAME = "SIGNAL_TO_START_ALL_JOBS";

my $TRACKFILE_SUFFIX = ".tracks.txt";

# List recognized "key" names from the config file. Anything in the config file but NOT in here will raise a warning.
my @OK_KEY_NAMES         = qw(studyName sampleDir resultDir minMapQ libraryAdapterFile bowtie2Index bwaIndex starIndexDir genomicBins genomicBinsID chromSizesFile geneWindows exonFile transcriptFile symbolXref genomeFasta genomeFastaPerChr genomeFastaDir repeatMask monkeyPoo aligner gtfFile tracksURL tracksRsync bedAnnotForRSEQC species version );
my @OK_BOOLEAN_KEY_NAMES = qw(forceRerun tracksShowEveryRead skipFastQC skipQC doQC skipBrowser browserBins browserWigs skipWindows doWindows skipTagsAndDensity doDensity shouldFilterAdapters doAltAnalyze );

my @OK_ALIGNERS  = qw(bowtie bowtie2 tophat tophat2 tophat2_1_1 star hisat hisat2_0_1 bwa);

my $SAFE_JOB_NAME_REGEXP = '^' . "[_A-Za-z0-9]+" . '$'; # Job names must be in this format. Note that hyphens are not allowed here.
my $SAFE_FILENAME_REGEXP = '^' . "[-+~._A-Za-z0-9]+" . '$'; # A filename, NOT a path, so no '/' are allowed. Dots, hyphens, and even the 'plus' sign are OK here.
my $SAFE_SPECIES_REGEXP = $SAFE_FILENAME_REGEXP; # species must either be undefined OR look something like this. Anything that would be OK in a filename (no spaces!) is ok here.

my $VERBOSE_DIAGNOSTIC_OK_COLOR = "cyan on_black";
my $VERBOSE_SKIP_STEP_OK_COLOR  = "magenta on_black";

# Currently running job names from qstat:   qstat -f| grep 'Job_Name' | perl -pe 's/^\s*Job_Name\s*=\s*//'

my $BANANAS_VERSION = "1.1.2017.01";

my $LOGDIR = "X.logs"; # ALL logs go into here -- note that this location is HARD CODED in each of the scripts in "monkey/bin/*.pl", so DO NOT CHANGE THIS unless you also change it in all those scripts too!

my $globalNumJobsSubmitted = 0; # counter

# =================== "REMEMBER" stores various file relationships (e.g., which BAM files were generated) which are not stored anywhere else ===================
my $REM_ALIGN_BAM           = "alignbam"; # where we remember which bam files we aligned
my $REM_EXPERIMENT_CATEGORY = "category"; # remember which experimental categories exist
my $REM_SCRIPTS_CHECKED     = "scripts-checked-for-ok-syntax";
my $REM_ALL_JOBS            = "alljobs";
my $REM_DEPENDENCIES_STR_COLON_DELIM     = "dep-colon-delim";
my $REM_SAMPLE_ORDER_ARR    = "sample-ids-in-order"; # Remember the original order that the files were specified in the config file. Important for grouping samples in the 'subread' pairs plot figure. This is an ARRAY (@) of the sample IDs (not filenames)!
my %remember = ($REM_ALIGN_BAM            =>{} # new empty hash ref
		, $REM_EXPERIMENT_CATEGORY=>{} # new empty hash ref
		, $REM_SCRIPTS_CHECKED    =>{} # new empty hash ref
		, $REM_ALL_JOBS           =>{} # new empty hash ref
		, $REM_DEPENDENCIES_STR_COLON_DELIM    =>{} # new empty hash ref
    ); # A hash that stores important file relationships that we need to remember. Somewhat ad hoc at the moment.
@{$remember{$REM_SAMPLE_ORDER_ARR}} = (); # New empty array
# ==============================================================================

my $STEP_FILTER               = "s02_filter";       # <== These need to appear ALPHABETICALLY,
my $STEP_MAPPING              = "s04a_map";         #     as that is how jobs are queued up and run.
my $STEP_TAGS                 = "s05a_tag";         # In other words: jobs can only depend on other ones that come ALPHABETICALLY EARLIER.
my $STEP_TAG_STATS            = "s05b_tag_stats";
my $STEP_DENSITY              = "s06_density";
my $STEP_PEAKS                = "s07a1_peaks";
my $STEP_PEAK_SUMMARY         = "s07a2_peaks_summary";
my $STEP_GEM_PEAKS            = "s07a3_gemPeaks";
my $STEP_GEMPEAKS_SUMMARY     = "s07a4_gemPeaks_summary";
my $STEP_BCP_PEAKS            = "s07b1_bcpPeaks";
my $STEP_BCPPEAKS_SUMMARY     = "s07b2_bcpPeaks_summary";
my $STEP_NUCLEOATAC           = "s07c1_nucleoATAC";          # ATAC-seq
my $STEP_BROWSER_BINS         = "s08a_browser_bins";         # Note: Different from "s48_browser_wiggle"
my $STEP_BROWSER_BIN_SUMMARY  = "s08b_browser_sum";
my $STEP_WINDOWS              = "s09_windows";
my $STEP_MOTIFS               = "s10a_motif";
my $STEP_MOTIF_SUMMARY        = "s10b_motif_sum";
my $STEP_SUBREAD              = "s20a_subread_count";          # RNA
my $STEP_SUBREAD_SUMMARY_FIGS = "s20b_subread_figs";           # RNA
my $STEP_CUFFDIFF             = "s22_cuffdiff";                # RNA
my $STEP_CUMMERBUND           = "s23_cummerbund";              # RNA
my $STEP_RSEQC_SINGLE         = "s47a_rseqc";
my $STEP_RSEQC_TOGETHER       = "s46_rseqc_all";
my $STEP_BROWSER_WIGGLE       = "s48_browser_wig";
my $STEP_FASTQC_UNMAPPED      = "s49a_fastqc_unaligned";
my $STEP_FASTQC_ALIGNED       = "s49c_fastqc_aligned";
my $STEP_EDGER_DIFF_EXPR      = "s50a_edgeR";                   # RNA
my $STEP_MONOCLE              = "s51a_monocle";                 # RNA... single-cell really, though
my $STEP_SUMMARY_AFTER_EDGER  = "s51b_sum";
my $STEP_CLUSTER_AFTER_EDGER  = "s51c_clust";
my $STEP_ALT_ANALYZE          = "s52_alt_analyze";




# job times:
# tags: 8 hours
# alignment: 99 hours
# filering: 8 hours
# 05.xSummary.pl: 2 hours

#   07a.peaksGem.pl:#PBS -l walltime=8:00:00
#   07a.xSummaryGem.pl:#PBS -l walltime=1:00:00
#   07b.peaksBCP.pl:#PBS -l walltime=4:00:00
#   07b.xSummaryBCP.pl:#PBS -l walltime=1:00:00
#   07c.NucleoAtac.pl:#PBS -l walltime=300:00:00
#   07.peaks.pl:#PBS -l walltime=4:00:00
#   07.xSummary.pl:#PBS -l walltime=1:00:00
#   08.browser.pl:#PBS -l walltime=24:00:00
#   08.xSummary.pl:#PBS -l walltime=1:00:00
#   09.windows.pl:#PBS -l walltime=4:00:00
#   10.motifDiscovery.pl:#PBS -l walltime=4:00:00
#   10.xSummary.pl:#PBS -l walltime=24:00:00
#   20a.count.subread.pl:#PBS -l walltime=72:00:00
#   20b.summary.figures.pl:#PBS -l walltime=8:00:00
#   grep: data: Is a directory
#   22.cuffdiff.pl:#PBS -l walltime=300:00:00
#   23a.cummerbund.pl:#PBS -l walltime=8:00:00
#   46.qc.multiple.bam.pl:#PBS -l walltime=72:00:00
#   47.rse.qc.pl:#PBS -l walltime=100:00:00
#   48.browser_agw.pl:#PBS -l walltime=72:00:00
#   49.qc.pl:#PBS -l walltime=2:00:00
#   50a.edgeR.pl:#PBS -l walltime=24:00:00
#   51a.monocle.R:#PBS -l walltime=24:00:00
#   51b.summary.R:#PBS -l walltime=24:00:00
#   51c.cluster.R:#PBS -l walltime=24:00:00
#   52.altanalyze.pl:#PBS -l walltime=24:00:00

#########################################################
# Code for the sub-scripts to (for example) evaluate arguments

sub getMostExtensivePythonPath() {
	my $current = (defined($ENV{PYTHONPATH}) ? $ENV{PYTHONPATH} : "");
	return join(':'
		    , ($current
			, "/wynton/group/gladstone/third_party/libpython2.7/lib64/python2.7/site-packages"
			, "/wynton/group/gladstone/third_party/libpython2.7/lib/python2.7/site-packages")
		   );
}


# envLooksTrue is also used by the sub-scripts in the bin directory--do not delete it!
sub envLooksTrue {
	# Works fine for multiple values: e.g. my ($t1, $t2) = envLooksTrue("t_one", "t_two");
	# Returns whether $ENV{'the_input_thing'} is a TRUE perl value, and NOT "false" or "undef" (case-insensitive).
	# Example correct usage:   my $boolean = envLooksTrue("varname");      # Example wrong usage:     my $wrongAnswer = envLooksTrue($ENV{"varname"} <- WRONG ) <== WRONG! Do not pass "$ENV" in here!
	# Looks like it can also return an array, if you want to check multiple items at once
	if (scalar(@_) == 0) { confess "ERR: Can't call envLooksTrue without at least one (string) argument."; }
	my @listOfBools = map { defined($_) && $ENV{$_} && (uc($ENV{$_}) ne "FALSE") && (uc($ENV{$_}) ne "UNDEF") } @_;
	if (scalar(@_) == 1) { return $listOfBools[0]; } # Return just the FIRST argument. Allows this:  my $var = envLooksTrue("something") to work properly. Otherwise that will fail (it will evaluate the return value of the list in a scalar context, thus returning the length of the list (1) instead of the value. This is confusing because Perl is weird!
	else {                  return (@listOfBools); } # Return an ARRARY if multiple arguments are specified. Example:   my ($t1, $t2) = envLooksTrue("t_one", "t_two");
}

sub isNA($) { return (defined($_[0]) && (uc($_[0]) eq uc($NA_VALUE))); } # Checks for the literal "NA" text value. This is how we distinguish our own "intentionally not set" values from Perl undef values, which generally indicate a PROGRAMMING ERROR.
sub sampleHasMatchingInput($) { return(defined($_[0]->{'input'}) and !isNA($_[0]->{'input'}) and $_[0]->{'input'} ne ''); }
sub getMatchingInput($) { return (sampleHasMatchingInput($_[0])) ? $_[0]->{'input'} : $NA_VALUE; }

sub catRequiredFile { # same as catfile, but REQUIRES that the file exists, unless it's a dry run.
	my $filePath = catfile(@_);
	($GLOBAL_DRY_RUN or -e $filePath) or confess "[ERROR] We expected the following file to exist, but it did not: $filePath";
	return($filePath);
}
sub mkdirOrDie($) { # This is SURPRISINGLY complicated because of race conditions! 'mkdirOrDie' is also used by the sub-scripts in the bin directory--do not delete it!
	my ($dir) = @_;
	defined($dir)               or confess "[ERROR] What in the world! 'mkdirOrDie' got an undefined value passed to it! This is probably a programming error.";
	($dir !~ m/\s/)             or ourWarn("[WARNING] Creating a directory with WHITESPACE in the name. The offending name is: $dir");
	($dir =~ m/^[-+._~\w\/]+$/) or ourWarn("[WARNING] Creating a directory with NON-BASIC-ALPHANUMERIC (and [-+._~]) characters in it. The offending name is: $dir");
	if ($GLOBAL_DRY_RUN) {
		print STDERR "Dry run: we are not actually creating the directory named '$dir'.\n";
		return 1;
	}
	if (-d $dir) { return 1; } # ok, it exists already # race condition info: note: file could be created here BETWEEN checking for existence and making the directory! 
	if (not mkdir($dir) and ($! =~ m/exists/i)) { return 1; } # Actually MAKE the directory here. The failure message was "File exists" -- so this means that some OTHER process must have made it! So it's ok. 	# Note: mkdir returns 0 (FAILURE!!) if the directory already exists at this point, which is possible if some OTHER process made it
	if (-d $dir) { return 1; } # last chance to succeed!
	confess qq{[ERROR] Failed in our attempt to make (or find) the directory that we expected to be at the following location: "$dir". The error message was as follows: '$!'};
}


sub pythonHasLib($$) {
	my ($pyExe, $pyLibName) = @_;
	my $exitCode = system(("$pyExe", "-c", "import $pyLibName"));
	return ($exitCode == 0); # 0 == got it!
}
sub dieIfMissingPythonLib($$) {
	my ($pyExe, $pyLibName) = @_;
	(pythonHasLib($pyExe, $pyLibName)) or confess "[ERROR]: The required python library '$pyLibName' was NOT FOUND on the system! This is for the python executable '$pyExe'---double check to make sure that this is the version of python that you expect to be running! You may have multiple versions of python on your system (2.6, 2.7, 3.x, etc...). Try running 'import ___library_name_here___' on the command line, and seeing if you get the same problem. Also remember that the machine where you SUBMIT a job may not be the same machine that actually RUNS the job, if you are running on a cluster! You may be able to install the missing library using 'pip'.";
}

sub dieIfCannotWriteTo($) { # Argument: a file path. Tells you if you could theoretically write there or not.
	my ($path) = @_;
	($path =~ /^\//) or ourWarn("[WARNING] You appear to have given a non-absolute path to the 'dieIfCannotWriteTo()' function---with Monkey, everything should be an ABSOLUTE path. Continuing anyway..."); # check to make sure it starts with a '/' (and is therefore probably a full path)
	if (-d $path) { -r $path and -w $path or confess "[ERROR] Weird, the directory <$path> was not both readable and writable by the current user! Disaster! Fix this."; }
	if (-e $path) { -r $path and -w $path or confess "[ERROR] Weird, the file at <$path> was not both readable and writable by the current user! Catastrophe! Fix this."; }
	if (not -e $path) { # Nothing ALREADY exists at this location, so let's see if we could (hypothetically) create something here by looking to see if the parent directory is writable.
		my $pathParentDir = dirname($path); # always the parent dir, even if this path is ITSELF a directory
		if (not -e $pathParentDir) {
			ourWarn("[WARNING] This is weird... the parent directory of <$path> also doesn't exist. This may be ok, though, so we will continue... but it may be a problem, too.\n");
		} else { # see if the parent "directory" is actually a directory and is writeable
			if (-f $pathParentDir) { confess "[ERROR] The parent directory of <$path> cannot be an (existing) plain file, this is too bizarre."; }
			if (-d $pathParentDir) { -r $path and -w $path or confess "Hey, we need to be able to write to the parent directory <$pathParentDir> of path <$path>! But right now we cannot. Fix the permissions!"; }
		}
	}
	return 1;
}

sub fileIsReadableAndNonzero($) { my ($f)=@_; return (($GLOBAL_DRY_RUN) or (defined($f) and -e $f and -r $f and -s $f)); } # -e: exists. -r: is readable. -s: is nonzero

sub dieIfFileAccessFails($;$) {
	# Input: a filename. This function will do nothing UNLESS that file is missing or unreadable, in which case it should give you an informative message about the specific inaccessible file.
	if ($GLOBAL_DRY_RUN) { return 1; }
	my ($filename, $optionalMessage) = @_;
	if (!defined($optionalMessage)) { $optionalMessage = "FILE HANDLING FAILED"; }
	defined($filename)      or confess "[ERROR] $optionalMessage: For whatever reason, a required filename was passed in as 'undefined'. This probably means that you did not specify a required parameter/filename in the monkey config file. See if other messages give any additional indication of what happened!";
	if (-d $filename) { return 1; } # ok, it was a directory, so I guess we can just return early without checking for anything else!
	(-e $filename)          or confess "[ERROR] $optionalMessage: Cannot open the seemingly-nonexistent file '$filename'! Quitting now.";
	(-r $filename)          or confess "[ERROR] $optionalMessage: this user cannot READ the file '$filename'! Quitting now.";
	return 1; # Looks like it was OK
}


sub numLinesInFile { # pass in a FILENAME
    my ($filename) = @_;
    my $lineCount = 0;
    open(FFF,$filename) or die "Cannot open $filename\n"; {
	    while(<FFF>) { $lineCount++; }
    } close(FFF);
    return($lineCount);
}

sub bash_system { my @args = ( "bash", "-c", shift ); return(system(@args)); } # Runs system calls through BASH instead of regular SH, allowing for things like subshells.

# systemAndLog is also used by the sub-scripts in the bin directory--do not delete it!i
sub systemAndLog($;$$$) {
    # Takes a command and (optionally) a boolean (should log to STDERR?) and a filename (if defined, append to this log file?)
    # Returns the exit code of the command, just like plain 'system'
    my ($cmd, $shouldLogToStderr, $logToThisFile, $dieIfNonzero) = @_;
    if (!defined($shouldLogToStderr)) { $shouldLogToStderr = 0; }
    if (!defined($dieIfNonzero)) { $dieIfNonzero = 0; }
    my $date = `date`; chomp($date);
    my $msg      = ($GLOBAL_DRY_RUN ? "DRY RUN of this command: " : "[${date}] Running this command:") . "\n     ${cmd}\n";
    my $exitCode = ($GLOBAL_DRY_RUN) ? 0 : bash_system($cmd); # Only run the command if this is NOT a dry run.
    $msg .= "     (Returned exit code $exitCode)\n";
    if ($exitCode != 0) { $msg .= "[WARNING]: Exit code was nonzero!\n"; }
    if ($shouldLogToStderr) { print STDERR $msg; }
#    if (defined($logToThisFile)) {
#	open(LOG,">>$logToThisFile") or confess "ERROR: Failed to open log file $logToThisFile. Probably this is a programming error.";
#	if ($logToThisFile) { print LOG $msg; }
#	close(LOG) or confess "ERROR: Could not close the file $logToThisFile.";
#    }
    if ($dieIfNonzero && $exitCode != 0) { confess "[ERROR / FATAL SYSTEM CALL FAILURE]! Non-zero exit code:\n$msg"; }
    return($exitCode);
}

sub systemAndLogDieNonzero($;$$) { my ($cmd, $toerr, $f) = @_; return(systemAndLog($cmd, $toerr, $f, "die if nonzero exit")); }

sub setIfUndefined($$$) {
	my ($config, $varname, $defaultValue) = @_; # Sets a variable in the config only if that value was NOT defined before!
	if (not defined($config->{$varname})) { $config->{$varname} = $defaultValue; }
}

# requireEnvOrDie: This function is used by the SUB scripts (the ones in "poo/bin") in order to validate that %ENV variables were passed in correctly.
sub requireEnvOrDie { # Takes any number of string inputs
    # Example usage:
    #      bananas_agw::requireEnvOrDie('input', "secondary_thing", "outfilename"); # If any of those are undefined, 'confess'es with an error. Otherwise returns true.
    # This function can accept any number of (string) arguments.
    # It checks to make sure they are all DEFINED. If any one of them is NOT, then it exits with a FATAL ERROR.
    # Variables ARE allowed to be blank. But they must be *defined*.
    # Note: "%ENV" is automatically populated by the 'qsub' queue-handling command. It is not a hash that we create explicitly.
    my (@varList) = @_;
    #defined(%ENV) or confess "The global 'ENV' hash variable (which should have been passed in by the 'qsub' queue submission program) was *not* defined in script submission! This indicates a fatal error somewhere in the queue submission process. Possibly this job was run by some means OTHER than qsub?";
    for my $requiredVarName (@varList) {
	defined($ENV{$requiredVarName}) or confess "[ERROR] The required variable '$requiredVarName' was NOT passed into this script by bananas.pm! This is a bug in bananas--fix it! Also check your capitalization--variable names are CASE SENSITIVE.";
    }
    return 1; # <== If we get here, then all the required variables were defined.
}

sub removeTrailingSlash($) { # Removes (at most) ONE single trailing slash from a path (string)
    my ($s) = @_; $s =~ s/\/$//;
    return($s);
}

sub getNumExperimentalGroups() { return scalar(keys(%{$remember{$REM_EXPERIMENT_CATEGORY}})); } # number of unique GROUPS

sub getGroupLabelsAndBamStringByGroup($$$) { # Gets the bams **by group** from the "remember" hash.
    my ($labelDelim, $betweenGroupDelim, $withinGroupDelim) = @_;
    # Example:
    # Returns TWO STRINGS, something like:
    # labelStr:    WT,CTRL,EXP2 (<== delimited by $labelDelim, a comma in this case)
    # groupStr:    wt1.bam,wt2.bam,wt3.bam ctl1.bam,ctl2.bam exp2.1.bam,exp2.2.bam (<== delimited by spaces BETWEEN groups and commas WITHIN group here)
    my $bamsByGroupStr = ""; # <== Will eventually look something like "exp1.bam,exp2.bam  ctrl1.bam,ctrl2.bam,ctrl3.bam"
    foreach my $category ( keys( %{$remember{$REM_EXPERIMENT_CATEGORY}} ) ) {
	my @bamArray = ();
	foreach my $replicate (keys(%{$remember{$REM_EXPERIMENT_CATEGORY}{$category}})) {
	    my $bamPath = $remember{$REM_EXPERIMENT_CATEGORY}{$category}{$replicate}{"bam"};
	    push(@bamArray, $bamPath);
	}
	$bamsByGroupStr .= ((length($bamsByGroupStr) == 0) ? '' : $betweenGroupDelim) . join($withinGroupDelim, @bamArray); # <== add a space between categories for every category EXCEPT the very first one!. The FINAL result will look something like "exp1.bam,exp2.bam  ctrl1.bam,ctrl2.bam,ctrl3.bam"
    }
    my $labelStr = join($labelDelim, keys(%{$remember{$REM_EXPERIMENT_CATEGORY}})); # Should look something like: "wildtype<COMMA>control<COMMA>drug<COMMA>whatever". Note that we can't use a literal comma (,), because qsub chews it up. So we use "<COMMA>" instead. The cuffdiff script knows how to deal with this.
    return ($labelStr, $bamsByGroupStr); # TWO strings!!
}

#########################################################
sub safeColor($;$) { # one required and one optional argument
    ## Returns colored text, but only if $SHOULD_USE_COLORS is set.
    ## Allows you to totally disable colored printing by just changing $SHOULD_USE_COLORS to 0 at the top of this file
    # Colorstring is OPTIONAL, and can be something like "red on_blue" or "red" or "magenta on_green"
    # Example usage:
    #    *    print STDERR safeColor("This warning message is red on yellow", "red on_yellow");
    my ($message, $color) = @_;
    if ($GLOBAL_COLOR) { use Term::ANSIColor; }
    return (($GLOBAL_COLOR && defined($color)) ? (Term::ANSIColor::colored($message, $color) . Term::ANSIColor::color("reset")) : $message);
}

sub printColorErr($;$) {    # prints color to STDERR *UNLESS* it is re-directed to a file, in which case NO COLOR IS PRINTED.
    my ($msg, $col) = @_; # Only prints in color if STDERR is to a terminal, NOT if it is redirected to an output file!
    if (! -t STDERR) { $col = undef; } # no coloration if this isn't to a terminal
    print STDERR safeColor($msg, $col);
}

sub printColorOut($;$) {    # prints color to STDOUT *UNLESS* it is re-directed to a file, in which case NO COLOR IS PRINTED.
    my ($msg, $col) = @_; # Only prints in color if STDOUT is to a terminal, NOT if it is redirected to an output file!
    if (! -t STDOUT) { $col = undef; } # no coloration if this isn't to a terminal
    print STDOUT safeColor($msg, $col);
}

# aid with checking configuration file key values
sub addKeyValuePairToConfig($$$$$) {
    my ($cfgPointer, $key, $value, $lineNum, $cfgFilename) = @_;

    # Check for deprecated variable names.
    ($key ne 'tracksSCP')   or confess "[ERROR] 'tracksSCP' is now a DEPRECATED CONFIG VARIABLE NAME. Replace it with 'tracksRsync', which is superior for various reasons.";
    ($key ne 'skipBrowser') or confess "[ERROR] 'skipBrowser' is now a DEPRECATED CONFIG VARIABLE NAME. Instead, specify 'browserWigs = FALSE' and 'browserBins = FALSE' in the config file.";

    (!exists(${$cfgPointer}->{$key})) or confess "\n[ERROR] DUPLICATE ASSIGNMENT OF '$key' DETECTED (filename: \"$cfgFilename\")\n>>>> The '$key' was defined a SECOND time on line $lineNum. You cannot have a key defined twice in the same file! ";
    (grep(/^$key$/, (@OK_KEY_NAMES, @OK_BOOLEAN_KEY_NAMES))) or warn "\n>>>> KEY WARNING: the key '$key' (on line $lineNum in the config file $cfgFilename) was NOT FOUND in the list of recognized keys.";
    ${$cfgPointer}->{$key} = $value;
}

sub getBamFullPath($$) { # get the final path to a bam file in the "04 mapping" directory only
	my ($cfgPointer, $sampName) = @_;
	my $file = ${sampName} . "_" . $cfgPointer->{genomeName} . "_" . 'q' . $cfgPointer->{minMapQ} . ".bam"; # something like MUTANT1_hg19_q30.bam or WT2_mm9_q3.bam
	return catfile($cfgPointer->{'mappingResultsDir'}, $file);
}


# Same as regular "print" in most regards, except it ONLY actually prints if the $GLOBAL_VERBOSE flag is set. Only accepts one string argument. No arrays!
sub verbosePrint($;$) { # If it is a dry run, then it prepends [DRY RUN] to the statement.
    my ($text, $color) = @_;
    ($GLOBAL_VERBOSE || $GLOBAL_DEBUG) && printColorErr($text, $color);  # Prints to STDERR always!
}

sub verboseOkPrint($)   { verbosePrint("[OK] $_[0]"       , $VERBOSE_DIAGNOSTIC_OK_COLOR); } # Print that something is OK and doesn't really require the user to do anything.
sub verboseSkipPrint($) { verbosePrint("[OK] [SKIP] $_[0]", $VERBOSE_SKIP_STEP_OK_COLOR); } # Print that we are skipping a step. This is OK also.

sub warnPrint($) {
	# Prints a warning regardless of the "verbose" settings
	my ($text) = @_;
	printColorErr($text, "red");
	#printColorOut($text, "red");
	warn $text;
}

sub debugPrint($;$) {
    my ($text, $color) = @_;
    ($GLOBAL_DEBUG) && printColorErr($text, (defined($color)?$color:"yellow on_red"));  # Prints yellow on_red if no color was specified
}

sub isLiteralStringTrue($) {    # Returns "TRUE" if input is both defined and the string "TRUE" or just "T", in any capitalization, or 0 otherwise.
	my ($inputString) = @_;
	return((defined($inputString) && ($inputString =~ m/^(TRUE|T)$/i)) ? "TRUE" : 0);
}

sub isBooleanStringOrMissing($) { 
	my ($inputString) = @_;     # Returns "TRUE" if and only if the input is TRUE/true, FALSE/false, T/F or undefined. Used for checking validity of inputs.
	return((!defined($inputString) || ($inputString =~ m/^(TRUE|T|FALSE|F)$/i)) ? "TRUE" : 0);
}

sub getArrayOfAlignedBamFilenames($) {  # (uses a global variable: the %remember hash). Returns an ORDERED array of full FILE PATHS in the same order they were in the config file.
    # Example usage: my @a = getArrayOfAlignedBamFilenames("paired_end"); print join(", ", @a);
    # Returns the filenames IN THE SAME ORDER that they appear in the original config file!
    my ($filterFor) = @_; # Tells us if we want ALL files, or just single end / paired end. There are three options: "all" or "paired_end" or "single_end"
    ($filterFor =~ m/^(all|paired_end|single_end)$/) or confess "[ERROR] The 'filterFor' parameter must be either 'single', 'paired', or 'both'";
    my @bamArr = (); # Return array
    foreach my $id (@{$remember{$REM_SAMPLE_ORDER_ARR}}) {
	my $isPaired = $remember{$REM_ALIGN_BAM}{$id}{'isPaired'};
	if (($filterFor eq "all") or ($isPaired and $filterFor eq "paired_end") or (!$isPaired and $filterFor eq "single_end")) {
	    push(@bamArr, $remember{$REM_ALIGN_BAM}{$id}{'bam'});
	}
    }
    return (@bamArr);
}

sub getHashPointerOfAllBamJobs() { # (uses a global variable: the %remember hash)
    # Gets a valid "dependencies" hash that requires that all the BAM jobs have completed
    # This is used for any job that requires the BAM files be aligned, but doesn't have any other requirements.
    my $hashref = {}; # {} (hash reference)
    foreach my $sample (keys(%{$remember{$REM_ALIGN_BAM}})) {
	my $prereqJob = $remember{$REM_ALIGN_BAM}{$sample}{"job"}; defined($prereqJob) or confess "[ERROR] Failed to locate the required pre-requisite job variable for sample <$sample>!";
	$hashref->{$prereqJob} = $prereqJob; # save this as a prerequisite!
    }
    return $hashref; # hash reference to all BAM-generating jobs
}

sub noExtensions($) {  # Returns the filename but WITHOUT any compression extensions (including the "period" character)
    my ($filename) = @_;
    if (!defined($filename)) { return undef; }
    $filename =~ s/[.](bz2|bzip2|gzip|gz|zip)$//i; # Remove common compression extensions from filename/path.
    return($filename);
}

sub validateScriptOrDie($$) { # Validate that a script exists. Intentionally does not check that it is executable (-x), since we might run it like "perl ./script.pl" and thus the script itself need only be READABLE.
	my ($scriptPath, $cfgPtr) = @_;
	
	dieIfFileAccessFails($scriptPath, "This function called for the script <$scriptPath>, but there did not appear to a (readable) file actually there! Double check this! ");
	if ($scriptPath =~ /.pl$/i and (!$remember{$REM_SCRIPTS_CHECKED}{$scriptPath})) { # Is this a perl script that we have NOT checked to make sure its syntax is OK? If so, check it.
		(0 == system(("perl", "-c", $scriptPath))) or confess("[ERROR in syntax checking while setting up qsub jobs]: Looks like perl was not able to parse the syntax of the script <$scriptPath> (in other words, perl -c $scriptPath FAILED to return exit code zero. That script probably has an error that needs to be fixed");
		$remember{$REM_SCRIPTS_CHECKED}{$scriptPath} = 1; # GLOBAL variable. Remember that this script was checked---don't bother checking it again
	}

	if ($scriptPath =~ m/51b.summary.R$/i) { # hard-coded check for this one file to skip
		print STDERR "[OK - NOTE] svTools::lint is going to skip checking the '51b.summary.R' file. There is actually a bug in svLint that causes it to report the incorrect failure message: 'in grep(regex, attr(exprs, 'srcref')[[j]]) ... 'Missing ')''";
		return;
	}

	# ================ CHECK R SYNTAX. This is VERY ANNOYING in comparison to perl ================
	if ($scriptPath =~ /[.](R|Rscript)$/i and (!$remember{$REM_SCRIPTS_CHECKED}{$scriptPath}))      {
		print STDERR "Skipping R checks";
		return;
	}
#	{ # Is this an R script that we have NOT checked to make sure its syntax is OK? If so, check it.
#		# Another thing that should be checked is: capabilities()["X11"] <== if this is FALSE, then no plots will works
#		my $R_CHECKER_CODE = <<"R_STUFF";
## Below: this is R code in a perl HEREDOC
#if (!require("svTools")) {
#        print("You need to run install.packages('svTools') in R to allow syntax checking.");
#        quit(save="no", 2); # Exit code 2 will be used to mean 'svTools' is not installed
#} else {
#	lintOut  = svTools::lint(file="${scriptPath}"); # <== Note the PERL variable here, which gets expanded to the proper file name!
#	errFound = (is.data.frame(lintOut) && nrow(lintOut) > 0); # <== This is R code
#        if (errFound) { print("Problems found!"); print(lintOut); }
#        else          { quit(save="no", 0); } # exit code 0 = everything is OK
#}
#quit(save="no", 1); # <== exit with code #1 to indicate that improper syntax was detected by svTools
#R_STUFF
#		my $R_CHECKER_TEMPFILE = "tmp.R_checker_code.tmp";
#		open(RC, ">$R_CHECKER_TEMPFILE") or confess("[ERROR] Couldn't make a temporary R file in this directory! Tried to write (unsuccessfully) to the following file: $R_CHECKER_TEMPFILE");
#		print RC $R_CHECKER_CODE;
#		close(RC);
#		my $rverb    = ($GLOBAL_VERBOSE) ? " --quiet --vanilla " : " --quiet --vanilla --slave "; # --slave is for background tasks and means "don't echo R commands" (even less verbose than 'quiet')
#		my $exitCode = verboseSystem(getBinPath($cfgPtr, "R_EXE") . qq{ $rverb < $R_CHECKER_TEMPFILE});
#		($exitCode == 1) and confess("[ERROR] The R file '${scriptPath}' had a syntax error as detected by svTools::lint. Source the file yourself to find the error.");
#		($exitCode == 2) and confess("[ERROR] The R syntax code checker requires that 'svTools' is installed. This probably means you do not have your R_LIBS set up properly. If you are POSITIVE you have it set up properly, you should run install.packages('svTools') (possibly as root) from within R in order to install the useful syntax checker package.");
#		($exitCode != 0) and confess("[ERROR] The R syntax code checker found a syntax error in the script '${scriptPath}' Exit code was $exitCode.");
#		unlink($R_CHECKER_TEMPFILE) or warn("[WARNING] Weirdly, we could not delete the R temp file ($R_CHECKER_TEMPFILE) that we made...");
#	}
	# ================ DONE CHECKING R SYNTAX. ================
}

# Runs a system command, and prints it out if we are using "verbose" mode (i.e., monkey.pl --verbose)
# DRY RUN behavior: Does NOT actually run the system call if this is a dry run!
sub verboseSystem($) {
    my ($cmd) = @_;
    if ($GLOBAL_DRY_RUN) { verbosePrint(">>>>> DRY RUN system call >>>>> $cmd\n", "cyan" ); }
    else {                 verbosePrint(">>>>> System call >>>>> $cmd\n"        , "green"); }
    my $exitCode = ($GLOBAL_DRY_RUN) ? 0 : system($cmd); # Only run the actual system command if it is NOT a dry run!
    verbosePrint("      (Returned exit code $exitCode)\n", ($GLOBAL_DRY_RUN ? "cyan" : "green"));
    return($exitCode);    # <== Mandatory that we RETURN the system call result!
}

sub ourWarn($) { # Prints a warning to the console, AND ALSO updates the 'GLOBAL_WARNINGS' variables that gets printed at the end.
    my ($msg) = @_;
    chomp($msg); # remove any newline so that the warning message shows the line number
    warn($msg); # Print the warning to the screen when it occurs.
    $GLOBAL_WARNINGS .= "$msg\n"; # Add the message to the GLOBAL warnings, plus a newline. We will print the GLOBAL warnings again when the script is done running.
}

sub printGlobalWarnings() {
	if (length($GLOBAL_WARNINGS) > 0) {
		print STDERR "We encountered the following warning messages:\n";
		print STDERR $GLOBAL_WARNINGS . "\n";
	}
}

sub fatalConfigErr($$$) { my ($lineNum, $configFilename, $msg) = @_; confess "[ERROR ON LINE $lineNum OF $configFilename]: $msg"; }
sub configWarning($$$)  { my ($lineNum, $configFilename, $msg) = @_; ourWarn("[WARNING ON LINE $lineNum OF $configFilename]: $msg"); }

sub isSplicedAligner($) {
    my ($alignerName) = @_;
    if ($alignerName =~ m/^(tophat|star)/i) { return 1; } # Tophat, STAR = SPLICE-AWARE
    if ($alignerName =~ m/^(bowtie|bwa)/i) { return 0; } # Bowtie, BWA = NO SPLICING
    confess "[Error] 'isSplicedAligner' was passed an UNRECOGNIZED aligner executable name: '$alignerName'...";
}

sub loadConfig { # <== loads the configuration file information into the "$cfg" hash pointer
	my ($file) = @_;
	my $cfg;		# <== hash (hash pointer?)
	my %dupeCheckHash = (); # used to check whether an ID is a dupe
	my $lineNum       = 0;
	$cfg->{filename}  = $file; # Just save the config filename in case we need it later!
	open(CF,$file) or confess "[ERROR] Cannot open the configuration file \"$file\" It may not exist, or not be readable by this user!";
	while (my $line = <CF>) {
		chomp($line);
		$lineNum++;
		$line =~ s/\s+[#].*//; # <== Remove non-full-line comments! Note that comments need a space before them! Example: turns "mydata something    # comment here" into "mydata something"
		$line =~ s/[ \t\s]+/\t/g; # <== Replace RUNS of spaces/tabs with a single tab. Note that ALL SPACES ARE TURNED INTO TABS here! That means you can't have errant spaces in (say) filenames and whatnot.
		$line = trimWhitespace($line); # <== Remove any other leading / trailing whitespace.
		next if ($line =~ /^[#]/); # <== Skip any lines that START with a '#' (comment character
		next if ($line =~ /^(\s)*$/); # <== Skip any ENTIRELY whitespace lines, in other words, there's nothing left after our modifications above.
		my @t = split(/\t/,$line);
		if ($line =~ /^(sample|input)\s/i) { #(scalar @t == 5) or (scalar(@t) == 6)) { # should probably just check to see if the line starts with "^sample\t" or "^input\t"
			(scalar(@t) == 5 or scalar(@t) == 6) or fatalConfigErr($lineNum, $file, "The line began with 'sample' or 'input', but it did NOT have 5 or 6 elements! All lines beginning with 'sample' or 'input' must have exactly 5 (for single end) or 6 (for paired end) values. The text of the offending line was:\n$line");
			my $sampleOrInput  = lc($t[0]); # either "sample" or "input" -- the literal text!
			my $idName         = $t[1];	# like "myDrugTest.1"
			my $dataType       = lc($t[2]); # rna, chip, or exo. Lower-case it.
			my $inputField     = $t[3]; # if it's CHIP-SEQ, then the "input" is a corresponding file that matches this sample. Or it can be "NA"
			my $firstPairFile  = $t[4];
			my $secondPairFile = (scalar(@t) >= 6) ? $t[5] : undef; # check to see if it's paired-end...
			
			($dataType =~ m/^(chip|exo|rna|atac|other)$/i) or fatalConfigErr($lineNum, $file, "The sample *type* that you specified was not recognized. You specified \"$dataType\", but the type must be one of these: {'chip', 'exo', 'rna','atac', or other}). Fix it!");

			if ((!defined($secondPairFile) or isNA($secondPairFile)) and ($firstPairFile =~ m/[-._](R1|pair1)[-._]/)) {
				ourWarn("[WARNING] Although there was no second pair specified for the single-end (?) file named '$firstPairFile', the filename contained a hint (like '.R1' or '.pair1') that the file might actually be paired end. Double check that this file is not missing its mate pair!");
				sleep(1); # give the user some time to ponder the warning before it swooshes by
			}
			if ($firstPairFile =~ m/[-._](R2|pair2)[-._]/) {
				ourWarn("[WARNING] The filename '$firstPairFile' had text in it that made it look like a second-of-a-pair paired-end file (like '.R1' or '.pair1'), but it was listed as either single-end of as the first file in a pair. Double check that this file is really NOT a second-mate-pair paired-end file!");
				sleep(1); # give the user some time to ponder the warning before it swooshes by
			}

			($sampleOrInput =~ m/^(sample|input)$/) or fatalConfigErr($lineNum, $file, "The first item must be either the literal text 'sample' or 'input'. Instead it was: <$sampleOrInput>!");

			(!exists($dupeCheckHash{uc($idName)})) or fatalConfigErr($lineNum, $file, "The ID name for each sample must be UNIQUE, but the ID \"$idName\" on line $lineNum had already been used on line $dupeCheckHash{$idName}. Fix this!");
			$dupeCheckHash{uc($idName)} = $lineNum; # <== Save the line number that this ID occurred on. Used to report any erroneous duplicates. Save the UPPER CASE version of each string to avoid allowing only-case-differences in filenames.
			($idName =~ m/${SAFE_FILENAME_REGEXP}/) or fatalConfigErr($lineNum, $file, "The sample ID CANNOT have any special characters in it. The allowable characters are the ones that match this regular expression: $SAFE_FILENAME_REGEXP. The invalid ID was this: \"$idName\"!");
			# sample.ID0414.fas
			($idName !~ m/^[0-9]/) or fatalConfigErr($lineNum, $file, "The sample name CANNOT *start* with a number (it confuses certain R commands, especially in edgeR), but this one (incorrectly) does: <$idName>. Recommendation: add a letter or underscore to the beginning.");
			
			if ($idName =~ m/^[^.]+[.]id/i) { # looking for "some stuff that is not a dot, then '.id' (case-insensitive). E.g. DRUGX.idWHATEVER or DRUG_A.id.4234.24.24.2
				# There is a special form for replicate names to be non-numeric, if they have a name like:
				#     group_whatever.id_429424.424.23423 . In this case, the magic part is ".id"--anything is OK as an identifier if that appears, even more decimal points!
				#     There is no additional error checking for an id here; anything (unique) is ok--even multiple '.', except spaces and special characters.
			} else {
				# This is a NORMAL style of sample description: GROUPNAME.replicate_number
				($idName !~ m/[.].*[.]/  )         or fatalConfigErr($lineNum, $file, "The sample ID CANNOT have more than one decimal point in it, but this one does: \"$idName\"!");
				($idName =~ m/[.][0-9]+$/)         or fatalConfigErr($lineNum, $file, "The name \"$idName\" DID NOT end in a period and then a (replicate) number. Valid names should look like this: \"WT.1\" or \"DRUG.13\" or \"Ctrl.18\". Numbers do *not* necessarily need to be sequential. Note that there is also a special form allowed for arbitrary names, you just have to put the word 'id' after the dot. So for example: DRUG.id_999.88_77 is a valid replicate name: it's group DRUG, with replicate id_999.88_77. Even multiple dots are allowed in this form.\n");
			}
			my @id_replicate_split = split(/[.]/, $idName, 2); # Split into TWO parts only! CRITICAL to have the 2 here!
			$cfg->{details}->{$idName}->{name}      = $id_replicate_split[0]; # example: "WILDTYPE" or "CONTROL"
			$cfg->{details}->{$idName}->{replicate} = $id_replicate_split[1]; # example: "1" (for first replicate of WILDTYPE)
			($cfg->{details}->{$idName}->{name} =~ m/$SAFE_JOB_NAME_REGEXP/) or fatalConfigErr($lineNum, $file, "The sample name <$cfg->{details}->{$idName}->{name}> is NOT acceptable, because it must only consist of standard alphanumeric characters plus the '-' and '_'. No special characters are allowed, as they may result in problematic filenames downstream. Please change this name to something without special characters, spaces, or punctuation");
			$cfg->{details}->{$idName}->{type} = $dataType; # Note that '$dataType' is already lower-cased

			$cfg->{details}->{$idName}->{needsPeakCall} = ($dataType =~ m/^(---|chip|exo|atac|other)$/i);  # peak calling for everything but RNA
			$cfg->{details}->{$idName}->{needsGemPeaks} = ($dataType =~ m/^(---|chip|exo|----|-----)$/i);  # chip or exo only for some reason
			$cfg->{details}->{$idName}->{needsATAC}     = ($dataType =~ m/^(---|----|---|atac|-----)$/i);  # chip or exo only for some reason
			$cfg->{details}->{$idName}->{needsBcpPeaks} = ($dataType =~ m/^(---|chip|---|----|-----)$/i);  # chip only for BCP
			$cfg->{details}->{$idName}->{needsDE}       = ($dataType =~ m/^(rna|----|---|----|-----)$/i);  # only RNA gets expression value analysis
			$cfg->{details}->{$idName}->{needsExo}      = ($dataType =~ m/^(---|----|exo|----|-----)$/i);  # only EXO gets exo-seq value analysis

			# Do peaks is GLOBAL
			$cfg->{globalDoPeaks} = (($cfg->{details}->{$idName}->{needsPeakCall}) || (defined($cfg->{globalDoPeaks}) && $cfg->{globalDoPeaks}));

			if (!isNA($inputField) and ($inputField ne '') and ($inputField =~ /\S/)) { # If this field isn't blank or NA or has at least one non-space character (\S)
				$cfg->{details}->{$idName}->{input} = $t[3]; # some samples (e.g. for CHIPSEQ) have a corresponding "INPUT" sample name associated with them.
			}
			$cfg->{$sampleOrInput}->{$idName}->{$firstPairFile} = 1; # indicate that we have obtained a sample with this specific name
			push(@{$remember{$REM_SAMPLE_ORDER_ARR}}, $idName); # "remember" is a GLOBAL variable that stores the sample IDs in the same order they were seen.
			if (defined($secondPairFile) && !isNA($secondPairFile)) {
				$cfg->{$sampleOrInput}->{$idName}->{$secondPairFile} = 1;
			}  # indicate that we have obtained a paired end file!!!
		} elsif ($line =~ m/^([^=]+)=(.*)$/) { # look for something with an EQUALS SIGN between it, e.g. "genome = myGenome.fa"
			my $key   = trimWhitespace($1);
			my $value = trimWhitespace($2); (length($value) > 0) or fatalConfigErr($lineNum, $file, "The value for key '$key' was BLANK or otherwise invalid. Note that every assignment MUST have a value, even if it is a placeholder like 'NA'. The offending line was:\n   $line\nPlease double-check the config file");
			#verboseOkPrint("Found a key / value $key / $value for species '$cfg->{species}').\n");
			# ==================== Handle the species-specific prefixes ================
			if ($key =~ m/^\[              # <== leading bracket. Example:  [hg19 ]gtfFile = something
				      \s*([^\]]*)\s*   # <== ($1): whatever is inside the brackets and is NOT a space is the SPECIES (e.g., hg19)
				      \]               # <== closing bracket
				      \s*              # <== optional whitespace, which is NOT part of the key
				      (.*)             # <== ($2): everything else is part of the key
				     /x) {
				my $keyIsForThisSpecies = trimWhitespace($1);
				(length($keyIsForThisSpecies) > 0) or fatalConfigErr($lineNum, $file, qq{Seems like you were using the species-specific bracket notation for key '$key', but the species name (\"$keyIsForThisSpecies\") was blank (???). Double check this offending line in the config file:\n   $line\n});
				if (defined($cfg->{species}) and not isNA($cfg->{species})) {
					if ($keyIsForThisSpecies eq $cfg->{species}) {
						#verboseOkPrint("FOUND a key for species '$keyIsForThisSpecies' '$key' = '$value', since that is the current species (which is '$cfg->{species}').\n");
						$key = $2;
					} else {
						verboseOkPrint("Skipping a key for species '$keyIsForThisSpecies' '$value', since that is not the current species (which is '$cfg->{species}').\n");
						next; # Skip it, and do NOT add this key to the config hash!
					}
				}
			}
			# ==================== Done handling the species-specific prefixes ================
			($key !~ m/\s/)      or fatalConfigErr($lineNum, $file, "The key '$key' had WHITESPACE in the key name, which is not valid! The offending line was:\n   $line\nPlease double-check the config file");
			(length($key) > 0)   or fatalConfigErr($lineNum, $file, "The key was BLANK or otherwise invalid. The offending line was:\n   $line\nPlease double-check the config file");
			addKeyValuePairToConfig(\$cfg, $key, $value, $lineNum, $file);
		} elsif (scalar(@t) == 2) { # look for any line with TWO ELEMENTS on it (i.e. tab-delimited key and value)
			my $key = $t[0]; my $value = $t[1];
			addKeyValuePairToConfig(\$cfg, $key, $value, $lineNum, $file);
		} else {
			fatalConfigErr($lineNum, $file, "A line in the config file could not be interpreted. Perhaps it is incorrectly formatted, or has a 'weird' character in it somewhere. The offending line was:\n   $line\nPlease double-check the config file");
		}
	}
	close(CF);
	return($cfg);
}

sub trimWhitespace { my $str=shift; $str =~ s/^\s+|\s+$//g; return($str) }; # trims leading and trailing whitespace

sub definedInConfigOrDie($$) { my ($cfg, $itemName) = @_; defined($cfg->{$itemName}) or confess "ERROR in config: Configuration file does not contain a '$itemName' setting--please specify this setting in the config file! Also, double check that it is properly capitalized!"; return 1; }

sub dirExistsInConfigOrDie($$) {
    my ($theCfg, $v) = @_;
    definedInConfigOrDie($theCfg, $v) and (-d $theCfg->{$v}) or confess "ERROR in config: The specified '$v' directory (" . $theCfg->{$v} . ") does not exist. Double check that the directory exists and is properly capitalized! ";
    return 1;
}

sub fileExistsInConfigOrDie($$) {
    my ($theCfg, $v) = @_;
    definedInConfigOrDie($theCfg, $v) and (-f $theCfg->{$v}) or confess "ERROR in config: The specified '$v' file (expected to be at \"" . $theCfg->{$v} . "\") does not exist. Double check that the file exists and is properly capitalized! ";
    return 1;
}

sub printInstallInstructions() {
	my $tempCfg;
	setupBin($tempCfg);
	print "ok\n";
}

sub addExe($$$$$$) {
	# The intent of this function, which is not actually used yet, is to list all the required executables
	# ...AND ALSO the command used to install them. The idea is that we could then easily install Monkey's dependencies
	# ...on a new machine by just iterating through the '$installCmd' variable here.
	my ($cfgPtr, $key, $dir, $exe, $desc, $installCmd) = @_;
	(!defined($cfgPtr->{bin}->{$key})) or confess "[ERROR] Attempt to DOUBLE DEFINE the key '$key'.\n";
	(defined($cfgPtr)) and $cfgPtr->{bin}->{$key} = catfile($dir, $exe); # Save the path to the executable into this specific $key
	# Note: right now we aren't saving the exe, dir, desc(ription) or installation cmd
	#addExe($cfg, "bowtie2", $binDir, "bowtie2", "Bowtie 2 non-splicing aligner", "brew tap homebrew/science; brew update; brew install bowtie2");
	#addExe(undef, undef, undef, "linuxbrew", "linuxbrew software installation system (similar to apt-get)", "brew tap homebrew/science; brew update;");
}

my %VERIFIED_EXE_HASH = (); # hash of executables that we verified do exist and are runnable
sub getBinPath($$) {
	# This verifies that binaries exist, but only when they are ***actually about to be USED***.
	# That way, we don't have to have every program installed on every system, if those programs aren't used for anything.
	my ($cfg, $exeName) = @_;
	exists($cfg->{bin}->{$exeName}) or confess "[ERROR] 'getBinPath' failed for the binary named '$exeName'--this is a programming error and indicates that we haven't added '$exeName' to the 'bin' hash yet in bananas_agw.pm. Fix this!";
	if (!defined($VERIFIED_EXE_HASH{$exeName})) {
		my $exePath = $cfg->{bin}->{$exeName};
		(-e $exePath) or confess "[ERROR in bananas_agw.pm]: While checking all the EXEs in cfg->{bin}, we discovered that the hard-coded path '$exePath' (for executable '$exeName') did not seem to have a valid executable! Check to make sure this file actually exists. If that is the wrong location for this executable, you may have to change the expected location in the code in 'bananas.pm'. Note that the problem may also be with your 'BINFBINROOT' enviroment variable (which should be a directory full of binaries). That directory is currently specified as <$GLOBAL_BIN_DIR>.";
		if ($exePath =~ m/[.](pl|sh|py)/) {
			# no need to check these for executable-ness, because they will be run by another program (perl, bash, python) anyway
		} else {
			(-x $exePath) or confess "[ERROR in bananas_agw.pm]: We discovered that the hard-coded path (for executable '$exeName') was present, but NOT EXECUTABLE by the current user! Use 'chmod' to investigate this! ";
		}
		$VERIFIED_EXE_HASH{$exeName} = 1; # Verified to exist / be an executable file
	}
	return($cfg->{bin}->{$exeName}); # returns the valid path
}

sub setupBin($) {
	my ($cfg) = @_;
	(scalar(@_) == 1) or confess "Wrong number of args to setupBin. The only argument must be the '$cfg' pointer to the hash where we store everything.";
	# Anything in $cfg->{bin} gets checked to make sure it both EXISTS and also CAN BE EXECUTED by the current user.
	# Do NOT put non-executables in $cfg->{bin}, or the automatic checking will fail!
	defined($cfg->{monkeyPoo}) or confess "[ERROR in code]: 'monkeyPoo' directory must have already been set up in the 'checkConfig' function BEFORE the binary directories can be set up.";

	my $binDir = $GLOBAL_BIN_DIR;
	my $MONKEY_SUPPORT_SUBDIR = catfile($cfg->{monkeyPoo}, "support");
	$cfg->{bin}->{basedir}                   = $binDir;
	$cfg->{bin}->{convertToGenomeBrowserExe} = catfile($binDir,"convert_SAM_or_BAM_for_Genome_Browser.pl");
	$cfg->{bin}->{qsub}                      = catfile($binDir,"qsub"); # queue submission utility
	$cfg->{bin}->{qhold}                     = catfile($binDir,"qhold"); # queue 'hold' utility
	$cfg->{bin}->{qrls}                      = catfile($binDir,"qrls"); # queue 'release' utility (opposite of qhold)
	$cfg->{bin}->{bam2bed}                   = catfile($binDir,"bam2bed");
	$cfg->{bin}->{bigWig}                    = catfile($binDir,"bedGraphToBigWig");
	$cfg->{bin}->{bedmap}                    = catfile($binDir,"bedmap");
	$cfg->{bin}->{bedops}                    = catfile($binDir,"bedops");
	$cfg->{bin}->{bigBed}                    = catfile($binDir,"bedToBigBed");
	$cfg->{bin}->{bedtools}                  = catfile($binDir,"bedtools");
	$cfg->{bin}->{bowtie2}                   = catfile($binDir,"bowtie2");
	$cfg->{bin}->{bwa}                       = catfile($binDir,"bwa");
	$cfg->{bin}->{STAR} = $cfg->{bin}->{star}= catfile($binDir,"STAR"); # <== allow both upper AND lower case for the STAR aligner, which is officially capitalized.
#	$cfg->{bin}->{hisat2_0_1}                = catfile($binDir,"hisat_2_0_1");   # for requesting a specific version
#	$cfg->{bin}->{hisat}                     = catfile($binDir,"hisat_current");
#	$cfg->{bin}->{tophat2_1_1}               = catfile($binDir,"tophat_2_1_1");  # for requesting a specific version
	$cfg->{bin}->{tophat}                    = catfile($binDir,"tophat");
	$cfg->{bin}->{cuffdiff}                  = catfile($binDir,"cuffdiff");
	$cfg->{bin}->{fastqMcf}                  = catfile($binDir,"fastq-mcf");
	$cfg->{bin}->{fastqc}                    = catfile($binDir,"fastqc");
	$cfg->{bin}->{gem}                       = catfile($binDir,"gem.jar");
        $cfg->{bin}->{bcp}                       = catfile($binDir,"BCP_HM");
	$cfg->{bin}->{match}                     = catfile($binDir,"match");
	$cfg->{bin}->{match2moa}                 = catfile($binDir,"match2moa");
	$cfg->{bin}->{moaOverlaps}               = catfile($binDir,"removeMoaOverlaps");
	$cfg->{bin}->{samtools}                  = catfile($binDir,"samtools");
	$cfg->{bin}->{sortBed}                   = catfile($binDir,"sort-bed");
	$cfg->{bin}->{subreadFeatureCounts}      = catfile($binDir,"featureCounts");
	$cfg->{bin}->{featureCountsToPlainMatrix}= catfile($MONKEY_SUPPORT_SUBDIR, "featureCounts_to_plain_matrix.sh");

	$cfg->{bin}->{alt_analyze_py}            = catfile($binDir,"AltAnalyze.py");

	$cfg->{bin}->{peakMotifs}                = "/wynton/group/gladstone/third_party/rsat/rsat/perl-scripts/peak-motifs"; # must be hardcoded to work... dubious. Also requires Perl library "MIME::Lite"

	$cfg->{bin}->{R_EXE}                     = "/wynton/group/gladstone/third_party/monkey_path/R";
	$cfg->{bin}->{RSCRIPT_EXE}               = "/wynton/group/gladstone/third_party/monkey_path/Rscript"; # Some sub-scripts (like 51c.cluster.R) need this, you but should not need to run it manually.
	$cfg->{bin}->{python2_7_exe}             = "/wynton/group/gladstone/third_party/monkey_path/python2.7"; # location of a python executable with version >= 2.7

#	$cfg->{bin}->{atacPython}                = "/home/sthomas/envs/atac/bin/python";
#	$cfg->{bin}->{nucleoatac}                = "/home/sthomas/envs/atac/bin/nucleoatac";

	$cfg->{bin}->{step1Motif}                = catfile($cfg->{monkeyPoo},"s1_getMotifRegex.pl"); # /home/sthomas/tools/utils/s1_getMotifRegex.pl
	$cfg->{bin}->{step2Motif}                = catfile($cfg->{monkeyPoo},"s2_findLocations.pl"); # /home/sthomas/tools/utils/s2_findLocations.pl
}

sub dieIfBothOptionsExist($$$) {
	my ($cfg, $a, $b) = @_;
	(!defined($cfg->{$a}) or !defined($cfg->{$b})) or confess "[ERROR] You have defined BOTH '$a' and '$b' in your config file! These are redundant and/or conflicting. Fix this in your config file!";
}

sub checkConfig {
    my ($cfg) = @_;
    # ================== Fix up the 'studyName' ================
    (defined($cfg->{studyName}))   or confess "[ERROR] Configuration file must contain a 'studyName' entry. This is just a basic name for your study, like 'myRnaSeqProject'. No periods / unusual characters are allowed! ";
    if ($cfg->{studyName} =~ /-/) {
	$cfg->{studyName} =~ s/-/_/g; # Changing all hyphens into UNDERSCORES, because hyphens are confusing to certain elements of the pipeline.
	print STDERR qq{[NOTE]: the hyphens in the 'studyName' will be converted into underscores--study name is now \"$cfg->{studyName}\".\n};
    }
    ($cfg->{studyName} =~ m/${SAFE_JOB_NAME_REGEXP}/) or confess "[ERROR] Study name must contain ONLY alphanumeric characters and the underscore. No periods! You can use underscores. The offending name was: \"$cfg->{studyName}\" ";
    if ($cfg->{studyName} =~ /^[0-9]/) { $cfg->{studyName} = "X".$cfg->{studyName}; ourWarn(qq{WARNING: the 'studyName' cannot START with a number (it will fail in qsub), so we added an 'X' to the beginning---study name is now \"$cfg->{studyName}\".}); }
    # ================== Fix up the 'studyName' ================

    dirExistsInConfigOrDie($cfg, 'monkeyPoo');
    if ($GLOBAL_DEV_MODE) { # MonkeyPoo directory should be changed to "/path/to/monkey/poo/dev_bin" instead of just "/path.../bin/"
	    $cfg->{monkeyPoo} = dirname($cfg->{monkeyPoo}) . "/dev_" . basename($cfg->{monkeyPoo});
	    dirExistsInConfigOrDie($cfg, 'monkeyPoo'); # check it AGAIN if we are in dev mode
    }
    setupBin($cfg); # set up the 'bin' directory. This must happen AFTER the 'monkeyPoo' directory is defined.

    dirExistsInConfigOrDie($cfg, "sampleDir");
    ($cfg->{sampleDir} =~ /^\//)  or confess "[ERROR] 'sampleDir' must be a FULL PATHNAME (so it should start with a '/' -- i.e., /path/to/place) and NOT a relative path! The offending name was: $cfg->{sampleDir} "; # check to make sure it starts with a '/' (and is therefore probably a full path)
    ($cfg->{sampleDir} !~ /[\s]/) or confess "[ERROR] 'sampleDir' must NOT contain whitespaces. The offending name was: $cfg->{sampleDir} ";

    if (defined($OUTPUT_DIR_OVERRIDE)) { $cfg->{resultDir} = $OUTPUT_DIR_OVERRIDE; } # <== the user can specifiy --out=/my/directory/path on the command line to override the path in the study design file.
    (defined($cfg->{resultDir}))  or confess "[ERROR] Configuration file must contain a 'resultDir' entry!";
    chomp($cfg->{resultDir}); # (Sometimes a newline can end up here in odd scenarios)
    length($cfg->{resultDir}) > 0 or confess "[ERROR] Result directory was blank somehow. This may be a result of a 'bad inode' situation---if you DEFINITELY set 'resultDir=something' in your config file, try changing to your home directory and then changing to this directory again (this fixes the 'bad inode' situation caused by directory renaming). "; # check to make sure it starts with a '/' (and is therefore probably a full path)
    ($cfg->{resultDir} =~ /^\//)  or confess "[ERROR] Result directory must be a FULL PATHNAME (i.e., /path/to/place), so that means it has to start with a '/'---it can NOT be a relative path! The offending name is in brackets here: [$cfg->{resultDir}] "; # check to make sure it starts with a '/' (and is therefore probably a full path)
    ($cfg->{resultDir} !~ /[\s]/) or confess "[ERROR] Working directory cannot contain whitespaces. The offending name was: <<<$cfg->{resultDir}>>>. Beware of non-printable characters and newlines! ";
    mkdirOrDie($cfg->{resultDir});

    $cfg->{writeup}                   = catfile($cfg->{resultDir}, "00a.writeup.txt");
    $cfg->{writeup_cite}              = catfile($cfg->{resultDir}, "00b.citations.txt");
    
    # Note: these are all FULL PATHS
    $cfg->{filterResultsDir}          = catfile($cfg->{resultDir},"02a.filtered_fastq");
    $cfg->{fastqcResultsDir}          = catfile($cfg->{resultDir},"02b.filtered_fastq_fastqc");
    $cfg->{mappingResultsDir}         = catfile($cfg->{resultDir},"04a.mapping");
    $cfg->{qcAfterMappingDir}         = catfile($cfg->{resultDir},"04b.mapped_qc_fastqc");
    $cfg->{rseqcQualityDir}           = catfile($cfg->{resultDir},"04c.mapped_qc_rseqc");
    $cfg->{tagResultsDir}             = catfile($cfg->{resultDir},"05.tags");
    $cfg->{densityResultsDir}         = catfile($cfg->{resultDir},"06.tagDensity");
    $cfg->{peakResultsDir}            = catfile($cfg->{resultDir},"07.peaks");
    $cfg->{gemPeakResultsDir}         = catfile($cfg->{resultDir},"07a.gemPeaks");
    $cfg->{bcpPeakResultsDir}         = catfile($cfg->{resultDir},"07b.bcpPeaks");
    $cfg->{atacResultsDir}            = catfile($cfg->{resultDir},"07c.nucleoAtac");
    $cfg->{browserBinResultsDir}      = catfile($cfg->{resultDir},"08.browser_tracks_binned");
    $cfg->{browserWigResultsDir}      = catfile($cfg->{resultDir},"08.browser_wig_and_bam");
    $cfg->{windowResultsDir}          = catfile($cfg->{resultDir},"09.window_figures");
    $cfg->{motifDiscDir}              = catfile($cfg->{resultDir},"10.motif_discovery");
    $cfg->{subreadCountsDir}          = catfile($cfg->{resultDir},"20.subread_counts");
    $cfg->{cuffdiffDir}               = catfile($cfg->{resultDir},"22.cuffdiff");
    $cfg->{cummerbundDir}             = catfile($cfg->{resultDir},"23.cummerbund");
    $cfg->{edgeRDir}                  = catfile($cfg->{resultDir},"50.edgeR_diff_expr");
    $cfg->{monocleDir}                = catfile($cfg->{resultDir},"51a.monocle");
    $cfg->{summaryFigureDir}          = catfile($cfg->{resultDir},"51b.summary_figures");
    $cfg->{clusterFigDir}             = catfile($cfg->{resultDir},"51c.cluster_heatmap");
    $cfg->{altAnalyzeDir}             = catfile($cfg->{resultDir},"52.alt_analyze");

    # Anything in $cfg->{r} gets checked for R syntactic validity, so DO NOT put non-R scripts there!
    $cfg->{r}->{pairPlot}          = catfile($cfg->{monkeyPoo}, "20c.summary.scripts", "R_pair_plot.R");
    $cfg->{r}->{cummerbund}        = catfile($cfg->{monkeyPoo}, "23b.cummerbund.R");
    $cfg->{r}->{edgeR}             = catfile($cfg->{monkeyPoo}, "50b.edgeR.R");
    $cfg->{r}->{monocle}           = catfile($cfg->{monkeyPoo}, "51a.monocle.R");
    $cfg->{r}->{window_figures}    = catfile($cfg->{monkeyPoo}, "09.xAnalyze.RScript");

    $cfg->{subread}->{countsFile}     = "subread.featureCounts.counts.txt"; # <-- literal filename that will be looked for later
    $cfg->{edgeR}->{fpkmFilePath}     = catfile($cfg->{edgeRDir}, "fpkm_naive.txt"); # note: hard coded in the edgeR script for now
    $cfg->{RSEQC}->{rseqc_parent_dir} = "/wynton/group/gladstone/third_party/rseqc/current/build/scripts-2.7"; # "/data/applications/rseqc/rseqc-2.6.1/build/scripts-2.7/"; # Location of the RSEQC python files

    # Set default values for certain unspecified variables.
    setIfUndefined($cfg, 'shouldFilterAdapters', "TRUE");    # Default is to filter adapters. If you don't want this, then specify "shouldFilterAdapters = FALSE".
    setIfUndefined($cfg, 'genomeFasta'         , $NA_VALUE);
    setIfUndefined($cfg, 'tracksShowEveryRead' , "FALSE");   # Default value: do NOT show each read (in other words, generate individual-read BAM files) for the 'browserWigs' part of the Genome Browser (s48_browser_wiggle). For that step, instead, only generate "wiggle" tracks. This is because BAM files are very large. Does not affect the generation of browser 'bin' tracks (step 8, s08a_browser_bins), which is unrelated.

    if (defined($cfg->{skipFastQC})) { fatalConfigErr("(No line number)", $cfg->{filename}, "ERROR: Sorry, the 'skipFastQC' option has been deprecated. Instead of specifying skipFastQC=FALSE, change it to doQC=TRUE . Sorry for the inconvenience!") }

    dieIfBothOptionsExist($cfg, 'skipQC'            ,'doQC');
    dieIfBothOptionsExist($cfg, 'skipWindows'       ,'doWindows');
    dieIfBothOptionsExist($cfg, 'skipTagsAndDensity','doDensity');

    # Handle the older style of "skip a thing that isn't specified"
    if (defined($cfg->{skipQC}))      {                     $cfg->{doQC} = isLiteralStringTrue($cfg->{skipQC})             ? "F" : "T"  } # yes, the literal strings "T" and "F"--not perl values yet!
    if (defined($cfg->{skipWindows})) {                $cfg->{doWindows} = isLiteralStringTrue($cfg->{skipWindows})        ? "F" : "T"  } # also note that it's BACKWARDS: skipQC=T means doQC=F
    if (defined($cfg->{skipTagsAndDensity}))        {  $cfg->{doDensity} = isLiteralStringTrue($cfg->{skipTagsAndDensity}) ? "F" : "T"  }

    if (!defined($cfg->{doQC})) { $cfg->{doQC} = "TRUE"; } # QC things if we didn't explicitly say not to!

    for my $boolOption (@OK_BOOLEAN_KEY_NAMES) { # Only run this AFTER we check for the old-style naming for 'skipQC' and the like, above
	isBooleanStringOrMissing($cfg->{$boolOption}) or confess "[ERROR] Configuration file problem: '$boolOption' was defined as \"$cfg->{$boolOption}\", but it must be the literal text TRUE (T), FALSE (F), or omitted entirely. Blank values are NOT ACCEPTABLE!";
	$cfg->{$boolOption} = isLiteralStringTrue($cfg->{$boolOption});
    }
    
    if (defined($cfg->{species})) {
	    (    $cfg->{species} =~ m/\S/    ) or fatalConfigErr("(No line number)", $cfg->{filename}, "\n\n[ERROR] The 'species' variable is defined, but is BLANK. Unacceptable! Set it to something.\n");
	    ($cfg->{species} =~ m/${SAFE_SPECIES_REGEXP}/) or fatalConfigErr("(No line number)", $cfg->{filename}, "\n\n[ERROR] The 'species' variable contains a SPACE or other illegal character! Unacceptable! Fix this. Right now, it is this: $cfg->{species}\n");
	    (not $cfg->{species} =~ m/(^__|__$)/) or fatalConfigErr("(No line number)", $cfg->{filename}, "\n\n[ERROR] The 'species' variable needs to be set to something VALID, but it looks like it's still our sample one (specifically, it is this: $cfg->{species}), since it has two underscores at the beginning and/or at the end---set it to the real species identifier without underscores at the beginning/end!\n");
    }

    (defined($cfg->{minMapQ})) or confess "[ERROR] MISSING CONFIG FILE OPTION: Configuration file must contain a 'minMapQ' entry! This is the minimum map quality (MAPQ) that a read must posses in order to NOT be disqualified after alignment. It should be between 0 (meaning 'keep everything') and 100. A standard value is 30.";
    ($cfg->{minMapQ} =~ /^[0-9]+$/ and $cfg->{minMapQ} >= 0 and $cfg->{minMapQ} <= 100) or confess "[ERROR] The 'minMapQ' option in the config file must be between 0 and 100, inclusive. The invalid map quality score specified in the file was the following: $cfg->{minMapQ}\n  ";

    if ($cfg->{'shouldFilterAdapters'}) { fileExistsInConfigOrDie($cfg, 'libraryAdapterFile'); } # If the user wants to filter adapters, make sure the adapter file exists!
    else { $cfg->{'libraryAdapterFile'} = "<Since we are not filtering adapters--this file should NOT be specified or used ever. Using this text for anything real is an ERROR.>"; }
    isNA($cfg->{'genomeFasta'})          or fileExistsInConfigOrDie($cfg, 'genomeFasta');
    
    if ($cfg->{'forceRerun'}) { print STDERR "[NOTE] Forcing re-run even of completed steps, due to the 'forceRerun=TRUE' line in the config file."; $FORCE_RERUN = 1; } # The user can either specify '-f' on the command line, or set 'forceRerun=TRUE' in the config file.
    
    if (!defined($cfg->{'bedAnnotForRSEQC'})) { ourWarn("Note: the 'bedAnnotForRSEQC' parameter was NOT DEFINED in the config file. We will skip this QC step until you define that file (or manually set it to 'NA')."); }
    if ($cfg->{doQC} && (!isNA($cfg->{'bedAnnotForRSEQC'})) && defined($cfg->{'bedAnnotForRSEQC'})) {
	($cfg->{'bedAnnotForRSEQC'} =~ /.bed$/i) or confess "[ERROR] The 'bedAnnotForRSEQC' file ('$cfg->{bedAnnotForRSEQC}') did not end with .bed, so it's probably not the right file type! Fix this (or set it to 'NA', or add doQC=FALSE to the config file if you don't care about QC at all)!";
	dieIfFileAccessFails($cfg->{'bedAnnotForRSEQC'}, "Configuration file problem: unless you skip the QC step, you MUST have a bed annotation file defined as 'bedAnnotForRSEQC'--this file is necessary for rseQC to function! You either didn't define a variable, or you passed in the invalid file location '$cfg->{bedAnnotForRSEQC}'! Check your config file and verify that a file is in fact at that location.");
    }

    # ======================== CHECK THE ALIGNER ====================
    (exists($cfg->{'aligner'})) or confess "[ERROR] MISSING CONFIG FILE OPTION: The configuration file must specify an ALIGNER (the 'aligner') option, which was not specified!";
    $cfg->{'aligner'} = lc($cfg->{'aligner'}); # Convert the aligner name to LOWER CASE!
    $cfg->{'aligner'} =~ s/^bowtie$/bowtie2/i; # Plain 'bowtie' always refers to 'bowtie2'. Always refer to it as 'bowtie2' from now on.
    $cfg->{'aligner'} =~ s/^tophat2$/tophat/i; # Tophat, on the other hand, is always version 2, but is called 'tophat'. Always call it plain 'tphat' from now on.
    (grep(/^$cfg->{aligner}$/, @OK_ALIGNERS)) or confess "[ERROR] The configuration file specified the following UNRECOGNIZED aligner: \"$cfg->{aligner}\". We expect the aligner to be something like 'bowtie' or 'tophat' or 'bwa'.\n";

    # ======================== HANDLE THINGS SPECIFIC TO CERTAIN ALIGNER(S) ====================
    if ($cfg->{'aligner'} =~ m/tophat/i) {
	    (defined($cfg->{gtfFile})) or confess "[ERROR] MISSING CONFIG FILE OPTION: If you specify 'tophat' as your aligner, you must **EITHER** (1) specify a valid GTF file with the 'gtfFile' parameter **OR** 2) Specify 'gtfFile=${NA_VALUE}' in the config file to explicitly indicate that there is no GTF file.";
	    verifyThatTophatSegmentJuncsCanRun($cfg) or confess "[ERROR] Failed to successfully launch the 'segment_juncs' executable.";
    }

    if ($cfg->{'aligner'} =~ m/^(tophat|star)/i) {
	    # Aligner can use a GTF file
	    (isNA($cfg->{gtfFile})) or (dieIfFileAccessFails($cfg->{gtfFile}, "The specified GTF annotation file ('gtfFile' setting in the config file) '$cfg->{gtfFile}' was NOT readable on the filesystem. Tophat and STAR requires one of two things: 1) either this file exist, OR 2) 'gtfFile' is set to the special value of 'gtfFile=${NA_VALUE}' to indicate that there is no such file."));
	    # Validate GTF if it's present maybe?? validate_gtf.pl can maybe do this.
    }

    if ($cfg->{'aligner'} =~ m/^(tophat|bowtie)/i) {
	    definedInConfigOrDie($cfg, "bowtie2Index");
	    ($cfg->{bowtie2Index} !~ m/[.]bt2$/) or confess "[ERROR] Configuration file error---the 'bowtie2Index' variable has a '.bt2' suffix on it, which it is not supposed to! Remove any file suffixes from this variable! Example: CORRECT = /path/to/hg19 . INCORRECT = /path/to/hg19.1.bt2 <== do not do that. ";
	    dieIfFileAccessFails("$cfg->{bowtie2Index}.1.bt2", "ALIGNER IS MISSING THE REQUIRED 'bowtie2' FORMAT INDEX FILE (those index files that end in '.bt2'). You probably need to either change the specified bowtie2 index location OR run 'bowtie2-build' on your input genome fasta file to build these indexes in the first place.");
    }
    if ($cfg->{'aligner'} =~ m/^(bwa)/i) { # A bwa index has a TON of files that must exist.
	    definedInConfigOrDie($cfg, "bwaIndex");
	    foreach my $suf (".amb", ".ann", ".bwt", ".pac", ".sa") { dieIfFileAccessFails("$cfg->{bwaIndex}${suf}", "MISSING THE REQUIRED '${suf}' FORMAT BWA INDEX FILE! Remember that the bam index is just the PREFIX for the files. So do not include '.bwt' or anything in the bwaIndex variable."); }
    }
    if ($cfg->{aligner} =~ m/^(star)/i) {
	    definedInConfigOrDie($cfg, "starIndexDir");
	    dirExistsInConfigOrDie($cfg, "starIndexDir"); # <== it needs to also be a DIRECTORY. Example of a valid starIndexDir: /data/info/genome/mm9_ensembl_igenome_with_chr/STAR_index_mm9.chr/
	    foreach my $starfile ("SA", "SAindex", "Genome", "chrLength.txt") {
		    dieIfFileAccessFails("$cfg->{starIndexDir}/${starfile}", "MISSING THE REQUIRED '${starfile}' FILE that the STAR ALIGNER needs. Make sure this file exists!");
	    }
    }
    if ($cfg->{aligner} =~ m/^(hisat)/i) {
	    confess "[ERROR] HISAT is not supported yet";
    }

    setIfUndefined($cfg, 'globalDoPeaks'             , 0);
    setIfUndefined($cfg, 'bowtie2Index'        , $NA_VALUE);
    setIfUndefined($cfg, 'bwaIndex'            , $NA_VALUE);
    setIfUndefined($cfg, 'starIndexDir'        , $NA_VALUE);
    
    # ======================== CHECK THE INPUT files (for ChIP-seq usually)... if any) ====================
    if (defined($cfg->{"input"})) {
	foreach my $k (keys(%{$cfg->{"input"}})) {
	    #next if (isNA($k)); # The "NA" values don't need to be checked--they do not exist!
	    foreach my $k2 (keys(%{$cfg->{"input"}->{$k}})) {
		#next if (isNA($k2)); # <== it's OK if a corresponding input file is NA maybe?
		dieIfFileAccessFails(catfile($cfg->{sampleDir}, $k2), "[ERROR] The ChIP-seq 'input' file for key '$k' and '$k2' and <$cfg->{sampleDir}/${k2}> was missing or otherwise UNREADABLE");
	    }
        }
    }

    (defined($cfg->{'sample'})) or confess "[ERROR] Configuration file does not contain any 'sample' lines. e.g.:\nsample <sampleLabel> <fullPathToSampleFile>. Note that the samples have to start with the LITERAL string 'sample', not the name of the sample!";
    while ( my($k,$th2) = each( %{$cfg->{'sample'}} ) ) {
	    foreach my $sampkey (keys %$th2) {
		    dieIfFileAccessFails(catfile($cfg->{sampleDir}, $sampkey), "[ERROR] The 'sample' file for key '$k' and '$sampkey' <$cfg->{sampleDir}/$sampkey> was unreadable\n ...maybe your sample directory was set incorrectly? It was set to: \"$cfg->{sampleDir}\" \n... you should double check it");
	    }
	    my $th3 = $cfg->{details}->{$k};
	    if (defined($th3->{input})) {
		    my $tch = $cfg->{input};
		    defined($tch->{$th3->{input}}) or confess "[ERROR] The 'Input' line that corresponds to sample '$k' lists a matched input ($th3->{input}) that isn't actually specified in the configuration file.";
	    }
    }

    if (isNA($cfg->{genomeFasta})) {
	    confess "[ERROR] That's weird, we expect you to ALWAYS provide a genome fasta file (like 'hg19.fa' or something), but apparently you did not? Failing in the code here in bananas_agw.pm.";
    }

    $cfg->{genomeName} = basename($cfg->{genomeFasta});
    $cfg->{genomeName} =~ s/[.](fa|fasta)//; # infer the genome name from the (required) input FASTA file. Remove the .fa/.fasta file extensions.
    
    if (not $cfg->{doDensity}) {
	    $cfg->{globalDoPeaks}  and ourWarn("[WARNING]: You are skipping tags and density---but you have data that needs peaks (you have either chip-seq or exome-seq data)! By specifying 'doDensity = FALSE', you skip EVERYTHING that relies on this downstream (browserBins, Sean's motif finding, etc.)");
	    $cfg->{browserBins}    and ourWarn("[WARNING]: You are skipping tags and density---but you have data that needs peaks (you have either chip-seq or exome-seq data)! By specifying 'doDensity = FALSE', you skip EVERYTHING that relies on this downstream (browserBins, Sean's motif finding, etc.)");
	    $cfg->{doWindows}      and ourWarn("[WARNING]: You are skipping tags and density---but you did NOT skip the 'window generation' step. Window generation cannot occur if you skipped tags and density! By specifying 'doDensity = FALSE', you skip EVERYTHING that relies on this downstream (browserBins, Sean's motif finding, etc.)");
	    $cfg->{globalDoPeaks} = 0; # Can't do peak calling without tags/density...
	    $cfg->{browserBins}   = 0; # Can't make browserBins without tags/density...
	    $cfg->{doWindows}     = 0; # If you skip tags and density, you ALSO have to skip window-generation...
    }
    
    if ($cfg->{browserBins} or $cfg->{browserWigs}) {
	    (defined($cfg->{tracksURL}))   or confess "[ERROR] MISSING CONFIG FILE OPTION: browser track URL basename ('tracksURL') needs to be specified in the config file. Or, if you don't want browser tracks, add 'browserBins = FALSE' and 'browserWigs = FALSE' to your config file. You can also specify 'NA' if you don't want a valid URL to be generated.";
	    (defined($cfg->{tracksRsync})) or confess "[ERROR] MISSING CONFIG FILE OPTION: if you want browser tracks to be synced to a remote server, you need to specify a browser track rsync location basename ('tracksRsync') in the config file. OR you can just specify 'NA' to avoid copying the files entirely.";
	    ($cfg->{tracksURL} =~ m/^https:\/\//) or ourWarn("WARNING: in 'tracksURL' specification: the 'tracksURL' parameter (\"$cfg->{tracksURL}\") should probably start with https:// . Note the 's' in httpS---not just http!");
	    if ($cfg->{browserBins}) {
		    dieIfFileAccessFails($cfg->{chromSizesFile}, "For making the 'browserBins' browser tracks, we NEED a valid chromSizesFile file, but the config file does not specify a valid one!");
		    dieIfFileAccessFails($cfg->{genomicBins}, "For making the 'browserBins' browser tracks, we NEED a valid 'genomicBins' file, but the config file does not specify a valid one!"); 	# Used by both Sean's browser and the tag density stuff
	    }
	    if ( defined($cfg->{tracksRsync}) && ($cfg->{tracksRsync} ne '') && !isNA($cfg->{tracksRsync}) ) {
		    # Test RSync to make sure there's a valid destination that we can write to
		    my $tempDir       = "z.test_rsync_permissions.tmpdir";
		    my $syncDest      = $cfg->{tracksRsync} . "/";
		    my $rsyncExitCode = verboseSystem(qq{mkdir -p "$tempDir"; rsync --chmod=u+rx,g+rx,o+rx --perms --dirs $tempDir  $syncDest}); # See if Rsync might run properly.
		    verboseSystem(qq{/bin/rmdir "$tempDir"}); # Remove the temp file once we get the exit code.
		    if (0 != $rsyncExitCode) {
			    fatalConfigErr("(No line number)", $cfg->{filename}, "\n\nERROR DETECTED IN YOUR RSYNC SETTINGS---we appear to possibly not have permissions (or some other error occurred) when we tried to write a test browser track file to this location: <$syncDest>.\n\nHOW TO SOLVE THIS:\nThis probably means that you are running as a non-admin user, so you'll need to probably change two variables: 'tracksURL' and 'tracksRsync' in your config file.\nIf you are on the machine 'Rigel', it should say: tracksRsync = /data/gbwww/browser/client/YOUR_NAME_HERE  and  tracksURL = https://gb.ucsf.edu/bio/browser/client/YOUR_NAME_HERE  .\nDouble check to make sure that the 'client' part is there and that and 'YOUR_NAME_HERE' is set to your username. Make sure to change BOTH 'tracksURL' and 'tracksRsync'!!\n\nHere is what your two variables are CURRENTLY set to:\n'tracksURL' is currently: " . $cfg->{tracksURL} . "\n" . "'tracksRsync' is currently: " . $cfg->{tracksRsync} . "\n\n");
		    }
	    }
    }

    if ($cfg->{doDensity}) {
	    dieIfFileAccessFails($cfg->{genomicBins}, "Configuration file does not contain a valid 'genomicBins' file, but this is required since we aren't skipping the 'tags and density' step!"); 	# Used by both Sean's browser and the tag density stuff
	    dieIfFileAccessFails($cfg->{genomicBinsID}, "Configuration file does not contain a valid genomicBinsID file, but this is required since we arne't skipping the 'tags and density' step!");
    }

    ##(open(SF,$cfg->{exonFile})       && close(SF)) or die "Configuration file does not contain a valid exonFile file!";
    ##(open(SF,$cfg->{transcriptFile}) && close(SF)) or die "Configuration file does not contain a valid transcriptFile file!";
    ##(open(SF,$cfg->{symbolXref})     && close(SF)) or die "Configuration file does not contain a valid symbolXref file, or that file could not be read.!";

    if ($cfg->{doWindows}) {
	    (defined($cfg->{geneWindows}))              or confess "[ERROR] MISSING CONFIG FILE OPTION: Configuration file must contain a 'geneWindows' entry, since we aren't skipping the 'gene windows' step!";
	    dieIfFileAccessFails($cfg->{geneWindows}, "Cannot find the following geneWindows file (which was specified in the configuration file: '$cfg->{geneWindows}'. This is REQUIRED since we aren't skipping the 'gene windows' step.");
    }

    if ($cfg->{globalDoPeaks}) {
	    (!isNA($cfg->{genomeFasta}) && defined($cfg->{genomeFasta}))       or confess "[ERROR] Configuration file does not contain a valid 'genomeFasta' entry (which is MANDATORY if you want the 'peak calling' part of the analysis)! Specify one. ";
	    $cfg->{genomeFasta} =~ /\.(fa|fasta)(\.gz|\.bz2|\.gzip|\.bzip2)?$/ or confess "[ERROR] The genome fasta file (which is MANDATORY if you want the 'peak calling' part of the analysis) DOES NOT end in '.fa' or '.fasta' (or a gzipped/bzipped version of this name. The offending name was: $cfg->{genomeFasta}). ";
	    dieIfFileAccessFails($cfg->{genomeFasta}, "Configuration file does have a valid 'genomeFasta' file (this must be a single fasta file!): it was $cfg->{genomeFasta}!");
	    (defined($cfg->{repeatMask}))              or confess "[ERROR] Configuration file does not contain a valid 'repeatMask' entry (which is MANDATORY if you want the 'peak calling' part of the analysis)! Specify one!";
	    dieIfFileAccessFails($cfg->{repeatMask} , "Configuration file does have a valid 'repeatMask' bed file (which is MANDATORY if you want the 'peak calling' part of the analysis). The missing/bad-permissoins filename was $cfg->{repeatMask}!");
	    dieIfFileAccessFails($cfg->{genomeFastaDir}, "If you are peak-calling with this tool, you absolutely MUST specify a 'genomeFastaDir', and it needs to be the path to a directory with the per-chromosome fasta files in it. For example, you could specify /path/to/hg19_by_chr/ , and inside that would be a bunch of files with names like 'chr1.fa' (or maybe just '1.fa'), 'chr2.fa', 'chr3.fa'...etc.");
    }

}

#########################################################
### 1. Build Job Submission Function ####################
#########################################################

sub sanitizeAndCheckJobName($) {
	my ($jobName) = @_;
	my $MAX_JOB_NAME_LENGTH = 220; # "The job name is up to and including 236 characters in length": PBS Pro documentation at http://www.pbsworks.com/pdfs/PBSProgramGuide13.0.pdf
	$jobName =~ tr/\./\_/; # replace any periods (such as those in the 'jobSample' (like "KO.1" or "Wildtype.2") by SINGLE underscores.
	($jobName =~ m/${SAFE_JOB_NAME_REGEXP}/) or confess "[ERROR] The job name must only contain 'normal' alphanumeric characters (and no spaces, hyphens, or periods)--yours had unusual characters in it (somewhere in this job name -->  $jobName   )";
	(length($jobName) <= $MAX_JOB_NAME_LENGTH) or confess "[ERROR] The job name is TOO LONG, probably because your study name is really long. The offending name was " . length($jobName) . " characters, whereas the maximum is $MAX_JOB_NAME_LENGTH . Here is the offending name: $jobName";
    return $jobName;
}

sub setJobInfo {
    my ($cfg,$jobWithoutStudyName,$jobSample,$dep,$vars,$theBinDir,$script) = @_;
    (length($jobWithoutStudyName) > 0) or confess "[ERROR] setJobInfo: Missing job 'base' name ('jobWithoutStudyName' was blank)";
    (length($jobSample) > 0) or confess "[ERROR] setJobInfo: Missing job 'sample' name ('jobSample' was blank)";
    defined($cfg->{studyName}) or confess "[ERROR] No study name was passed in via the 'cfg' variable! Probably it was called with a wrong argument in the function that called this, which you can see on the backtrace: ";
    my $jobFullName  = sanitizeAndCheckJobName("$cfg->{studyName}__${jobWithoutStudyName}__${jobSample}"); # e.g. "mouseStudy__10_alignment__allsamples" or "mouseStudy__04_quality_control__sample47_b"
    # ============ SET THE DEPENDENCIES ======================
    my $dependencies = "";
    if ($HOLD_UNTIL_ALL_JOBS_SUBMITTED) { # <== global variable
	    # This is an annoying hack---basically, PBS does not allow jobs to depend on jobs that have ALREADY FINISHED,
	    # even if they are still in the queue as 'COMPLETED'. So we make a fake 'start signal' job at the beginning, then hold that job.
	    # EVERY other job now depends on the 'start signal' job, which only gets released (= "qrls") at the end, when all the other jobs have been
	    # added to the queue. So NO jobs can run until every job has been submitted.
	    $dependencies .= '$' . "${SIGNAL_TO_START_VARNAME}";
    }
    foreach my $dkey (sort(keys(%$dep))) {  # <== Now list the dependencies in SORTED ORDER
	    my $dval = $dep->{$dkey};
	    next if (isNA($dval)); # Skip any "NA" prerequisites.
	    (defined($dval) && ("" ne $dval)) or confess "[ERROR making qsub jobs]: The dependency for job '$jobFullName',\n...specifically in step '$dkey', was undefined!\n...This probably indicates an error in the dependency names (something like 'step2' instead of 'step02' for example)!";
	    ($jobFullName ne $dval)            or confess "[ERROR making qsub jobs]: You set the job name ($jobFullName) to THE SAME NAME as one of its prerequisites! This is not valid--a job cannot depend on itself. Fix it! ";
	    (exists($remember{$REM_ALL_JOBS}->{ $dval }) && defined($remember{$REM_ALL_JOBS}{$dval})) or confess "[ERROR] A dependency was added to '$dval', but that does NOT appear to have been a valid job name that had been previously added!";
	    $dependencies .= ',' . '$' . $dval; # note the LITERAL DOLLAR SIGN, which makes this dependency into a shell variable name (yes, really, that is how they are stored).
    }
    $dependencies =~ s/^[:]//; # Remove any extraneous leading colon, if there is one.
    # ============ SET QSUB VARIABLES with -v ==================
    my $variables = "-v bananas='${BANANAS_VERSION}'"; # "-v" option supplies the variable list. Note that there is ALWAYS at least one variable (the bananas version)
    my $expvars = "";
    foreach my $key (sort(keys(%$vars))) {
	my $value = $vars->{$key};
	# "$NA_VALUE" is a valid value to pass on to the output, so do not omit it!
	(defined($value))    or confess "[ERROR making qsub jobs]: The variable for key '$key' was 'undef' (in other words, you never defined it, yet it WAS submitted to 'setJobInfo somehow'! This is not allowed! Figure out where this variable was passed in and why it wasn't set.";
	("" ne $value)       or confess "[ERROR making qsub jobs]: the value '$value' was set to the empty string (''). Currently we do NOT expect this to ever happen, so it probably indicates an error. If you need to indicate a null value, try using \$bananas_agw::NA_VALUE instead."; # We expect certain values (such as 'inputFile2' to INTENTIONALLY be blank, but most files should not be blank.
	($value !~ "[\'\`]") or confess "[ERROR making qsub jobs]: SINGLE QUOTES / BACKTICKS NOT ALLOWED]: There was a single quotation mark or backtick in the key-value pair $key='$value' (in other words, <$value> had an invalid character in it!)! This will break the qsub command--quotes and backticks are NOT allowed in inputs!";
	($value !~ "[,]")    or confess "[ERROR making qsub jobs]: COMMA NOT ALLOWED]: Torque/qsub is VERY PARTICULAR and cannot accept a comma in an argument (even inside double quotes). You have to REMOVE the comma, use some kind of placeholder text, and then manually re-add it in the script that you called. The specific value that seemed dubious (not including the quotes) was: \"$value\". I recommend using <COMMA> and then replacing that with a ',' in the being-called script.";
	#$variables .= qq{,}.qq{$key='$value'}; # <== Variables are comma-delimited and are surrounded by single quotation marks. This is buggy; we use '-V' and environment variables instead!
	$expvars   .= qq{export $key='$value'\n}; # Literally just using BASH exporting. More robust than "-v='VAR'". Only one problem---'export' commands from earlier commands can theoretically be seen by later ones. So beware!
    }

    my $scriptFullPath = catfile($theBinDir,$script);
    validateScriptOrDie($scriptFullPath, $cfg);
    my $depStr = ($dependencies ne "") ? "-hold_jid $dependencies" : ""; # Note that if $dependencies is non-blank, it will always begin with a ":"
#    my $depStr = ($dependencies ne "") ? "-hold_jid \$(echo $dependencies | awk '{print \$3}')" : ""; # Note that if $dependencies is non-blank, it will always begin with a ":"
    my $qsubCmd;
    if ($RUN_DIRECTLY_WITHOUT_TORQUE) { $qsubCmd = "${expvars}${scriptFullPath}"; } # <== note the literal backtick! There will be a matching one later. Also: the saving-to-shell-variables part is SPACE-SENSITIVE, so DO NOT add spaces here whatever you do.
    else {                              $qsubCmd = "${expvars}${jobFullName}=`" . get_qsub_cmd($cfg, $vars) . "-terse -V -N ${jobFullName} ${depStr} ${variables} ${scriptFullPath}`"; } # <== note the literal backtick! There will be a matching one later. Also: the saving-to-shell-variables part is SPACE-SENSITIVE, so DO NOT add spaces here whatever you do.
    $globalNumJobsSubmitted++;
    debugPrint("[DEBUG] Adding job ${globalNumJobsSubmitted} via qsub: $qsubCmd\n", "cyan");
    (!exists($cfg->{jobs}->{$jobWithoutStudyName}->{$jobSample})) or confess "bananas.pm [ERROR]: Trying to add a job that ALREADY EXISTS! Can't double-add a job with the same name AND sub-sample name! The illegal combination of job name / sample name was: \"$jobWithoutStudyName\" (job) and \"$jobSample\" (sample).";
    $cfg->{jobs}->{$jobWithoutStudyName}->{$jobSample} = { "jobName"=>$jobFullName, "qsub"=>$qsubCmd }; # <== '49' so it gets alphabetized after all ALIGNMENT is done, guaranteed. Jobs are done in alphabetical order.
    (!exists($remember{$REM_ALL_JOBS}{$jobFullName})) or confess "bananas.pm [ERROR]: more than one job had the name '$jobFullName'---this is a PROGRAMMING ERROR that needs to be fixed...";
    $remember{$REM_ALL_JOBS}{$jobFullName} = 1; # Remember that we have seen this job! Used for error checking the dependencies.
    $remember{$REM_DEPENDENCIES_STR_COLON_DELIM}{$jobFullName} = $dependencies; # Remember the dependency string. Looks like "dep1,dep2". Can be blank.
    return($jobFullName, $qsubCmd); # <== this info is used by a few things later
}

sub getUserPriv { # Gets user privilege. See "QSETTINGS"
	# No arguments!
	# Note that this is hard-coded for the particular system right now---if we are NOT on Rigel, then just return 100 for everyone.
	if (not(Sys::Hostname::hostname() =~ m/^rigel/i)) {
		# If we aren't on rigel, then EVERYONE should act like a 'privileged' user
		return $ADMIN_GROUP_NAME;
	}
	# Ok, I guess we are on the RIGEL machine.
	#my $username = $ENV{LOGNAME} || $ENV{USER} || getpwuid($<); chomp($username);
	my @gids = POSIX::getgroups(); # getgroups() is from the POSIX module
	for my $priority (@UNIX_BIO_GRP_ID_ARRAY) {
		if (grep(/^$priority$/, @gids)) { return $ADMIN_GROUP_NAME; } # apparently the user belongs to a privileged group--let them use more CPUs, etc.
	}
	return $NORMAL_GROUP_NAME;
}

sub get_qsub_cmd($;$) { # returns something like: qsub -q Bio -A "bioqueue" (with the appropriate groups set)
	my ($cfg, $hr) = @_; # hr is an (optional) hash ref of arguments to qsub! Example: "{h_rt=>"3600", mem_free=4g}"
	# -q General -A "genqueue"
	#    qsub -I -q Bio -A "bioqueue" -l scratch=5g
	#    qsub -q Bio     -A "bioqueue" [/path/to/script]
	# or qsub -q General -A "genqueue" [/path/to/script]
	#if (!defined($hr)) { my %h = (); $hr = \%h; }; # just define the hash ref if it wasn't actually passed in...
	my $qsub         = getBinPath($cfg,"qsub");
	my $userPriv     = getUserPriv();
#	my $qdest_param  = $QSETTINGS{dest}{$userPriv};
#	my $qgroup_param = $QSETTINGS{grouplist}{$userPriv};
	my $mem      = "-l mem_free="      . (defined($hr->{pbs_mem})      ? $hr->{pbs_mem}      : $QSETTINGS{default_mem}{$userPriv});
#	my $walltime = "-l h_rt=" . (defined($hr->{sge_h_rt}) ? $hr->{sge_h_rt} : $QSETTINGS{default_h_rt}{$userPriv});
	my $walltime = "-l h_rt=259200";
	my $ncpus    = "-pe smp "    . (defined($hr->{pbs_ncpus})    ? $hr->{pbs_ncpus}    : $QSETTINGS{default_ncpus}{$userPriv});
	#my $nodes    = "-pe smp 1"; # <== ????????? unclear if this is useful
	my $stderr   = ""; #"-e /dev/null"
	my $stdout   = ""; #"-o /dev/null"
	$mem      =~ m/^-l mem_free=\d+[gm]$/                or confess "[ERROR]: 'mem_free' parameter needs to be a number followed by 'g' or 'm'! Yours is: $mem";
	$walltime =~ m/^-l h_rt=\d+$/ or confess "[ERROR]: 'h_rt' parameter must be in units of seconds and look like this: 3600. Yours is: $walltime";
	$ncpus    =~ m/^-pe smp [1-9][0-9]*$/                   or confess "[ERROR]: 'ncpus' parameter must be numeric integer greater than 0! Yours is: $ncpus";
	return qq{$qsub $mem $walltime $ncpus $stderr $stdout};
}

sub verifyThatTophatSegmentJuncsCanRun($) {
	# Note: this only checks LOCALLY, not on the remote server!
	my ($cfg) = @_;
	my $tophat_binary = getBinPath($cfg, "tophat");
	my $segment_juncs_binary = `readlink -f $tophat_binary`; chomp($segment_juncs_binary); # find the TRUE location of tophat, not the symlink
	$segment_juncs_binary =~ s/tophat([\d.])*$/segment_juncs/; # should be wherever tophat is, but named something different
	#my $sj_exit_code = system($segment_juncs_binary);
	(-e $segment_juncs_binary) or confess "[ERROR] Could not find the 'segment_juncs' binary in expected location (alongside Tophat) in the following seemingly-nonexistent location: $segment_juncs_binary";
	my $sj_exit_text = `$segment_juncs_binary 2> /dev/stdout`; chomp($sj_exit_text);
	if ($sj_exit_text =~ m/Usage/) {
		# ok, looks good
	} else {
		confess "[ERROR]: 'segment_juncs' (required for 'tophat' to run) exited with a failure message!!. Tophat (which we found in <$tophat_binary) REQUIRES this component in order to run, so make sure you have it as a usable program.\n     Here is the segment_juncs message: $sj_exit_text\n     (If that says something about 'libboost', then it means your server needs to have this command run as root: 'yum install boost boost-devel'";
	}
	return 1;
}

sub buildJobSubmissionList {
	# This function organizes all jobs that will need to be run and creates a script
	# that will be run in order to process all of the jobs in the most efficient way possible.
	# To add a new job, specify the base job name, dependencies, variables, script directory, and called script.

	# Note: jobs need to be added IN ORDER--pre-requisites must come BEFORE their dependencies here!
	# The qsub scheduler is NOT clever enough to understand arbitrary ordering of prerequisites.
	# Maybe these need to explicitly be sorted alphabetically?
  	my ($cfg) = @_;
	if (-e $cfg->{writeup}  &&  (-s $cfg->{writeup} > 0)) { # check that the file exists AND is size > 0
		my $time = int(time());
		verboseSystem("mv " . $cfg->{writeup}      . "  " . $cfg->{writeup}      . ".old.${time}.backup.tmp");
		verboseSystem("mv " . $cfg->{writeup_cite} . "  " . $cfg->{writeup_cite} . ".old.${time}.backup.tmp");
	} # If there is already a writeup text file of non-zero size, move the old writeup to a backup file.
	my $sampleHash = $cfg->{sample};
	my $inputHash  = (exists($cfg->{input})) ? $cfg->{input} : undef; # Set to undef if there's no input
	# build filtering jobs
	foreach my $h ($sampleHash,$inputHash) {
		next if (!defined($h)); # Skip any undefined 'inputHash' cases.
		foreach my $k (keys %$h) {  # 'k' is the sample name--for example "DRUG.3" or "Wildtype.1"
			my @sKeys   = sort(keys(%{$h->{$k}})); # sorted (so the first pair for paired end should always come first!) list of input filenames. Length 1 (single end) or 2 (paired end) only.
			(scalar(@sKeys) == 1) or (scalar(@sKeys) == 2) or confess("[ERROR] in programming: somehow sKeys was not length 1 or 2, which are the only valid options! It cannot be 0 or 3, for example.");
			my $isPaired = (2==scalar(@sKeys)) ? 1 : 0;
			my $dep = {}; # <== Note: no dependencies! This is the first real step!
			my $vars = { 'force'=>$FORCE_RERUN, 'verbose'=>$GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG
				     , 'fastqMcf'             => getBinPath($cfg, "fastqMcf")
				     , 'fastqDir'             => $cfg->{'sampleDir'}
				     , 'filterDir'            => $cfg->{'filterResultsDir'}
				     , 'adapterFile'          => $cfg->{'libraryAdapterFile'}
				     , 'shouldFilterAdapters' => $cfg->{'shouldFilterAdapters'} # If this is FALSE, then we don't filter--we just symlink the input file to the output file
				     , 'inputFile1'           =>               $sKeys[0]
				     , 'inputFile2'           => ($isPaired) ? $sKeys[1] : $NA_VALUE # <== for paired end, get the second file this way. This MUST be specified as $NA_VALUE if there is nothing there! Do not omit it!
				     , 'downstreamBamFile'    => getBamFullPath($cfg, $k) # <-- the bam file that will supposedly EVENTUALLY be made from these filtered files
				   };
			(-e catfile($vars->{'fastqDir'}, $vars->{'inputFile1'})) or confess("PROGRAMMING [ERROR] [buildJobSubmissionList in bananas_agw.pm]: The file <$vars->{fastqDir}/$vars->{'inputFile1'}> was expected to exist (since it was in the configuration file), but we did NOT find it on the filesystem!");
			(!$isPaired) or (-e catfile($vars->{'fastqDir'}, $vars->{'inputFile2'})) or confess("PROGRAMMING [ERROR] [buildJobSubmissionList in bananas_agw.pm]: The file <$vars->{fastqDir}/$vars->{'inputFile2'}> was expected to exist (since it was in the configuration file), but we did NOT find it on the filesystem!");
			setJobInfo($cfg, $STEP_FILTER, "${k}",$dep,$vars,$cfg->{monkeyPoo},"02.filter.pl");
		}
	}
	writeup($cfg, qq{Trimming of known adapters and low-quality regions of reads was performed using Fastq-mcf [CITE_MCF].}, {"CITE_MCF"=>qq{Fastq-mcf (in EA-utils): Erik Aronesty (2011). "EA-Utils: Command-line tools for processing biological sequencing data"; http://code.google.com/p/ea-utils}});
	my $ncpus_for_alignment = $QSETTINGS{ncpus_for_alignment}{getUserPriv()};
	# build mapping jobs (map to the genome with tophat or bowtie)
	foreach my $h ($sampleHash,$inputHash) {
		next if (!defined($h)); # Skip any undefined 'inputHash' cases.
		foreach my $k (keys %$h) { # 'k' is the sample name--for example "DRUG.3" or "Wildtype.1"
			my @sKeys   = sort(keys(%{$h->{$k}}));
			my $dep     = { 'filterJob' => $cfg->{jobs}->{$STEP_FILTER}->{$k}->{jobName} }; # Depends on filtering being done already.
			my $isPaired       = (2==scalar(@sKeys)) ? 1 : 0;
			my $sampleType     = $cfg->{'details'}->{$k}->{'type'};
			my $sampleCategory = $cfg->{'details'}->{$k}->{'name'};	# exprimental category, like "WILDTYPE" or "CONTROL"
			my $sampleRep      = $cfg->{'details'}->{$k}->{'replicate'}; # Replicate number, like the "2" from "WILDTYPE.2" (which is "$k" here)
			#my $dataType       = $cfg->{'details'}->{$k}->{'type'}; # like 'exo' or 'rna' or 'chip'
			# ================= Set the aligner executable path, and warn the user about weird aligner choices ========
			my $aligner_requirement_string = undef;
			if ($cfg->{aligner} =~ m/^star/) { $aligner_requirement_string = "star_aligner_needs"; }
			else { $aligner_requirement_string = "tophat_needs"; }
			my $pbs_mem = $QSETTINGS{$aligner_requirement_string}{pbs_mem}; # default amount of ram needed is whatever tophat wants
			# =========================================================================================================
			my $vars = { 'force'=>$FORCE_RERUN, 'verbose'=>$GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG
				     , 'pbs_ncpus'   => $ncpus_for_alignment, 'pbs_mem' => $pbs_mem, 'pbs_walltime'=> "1076404"
				     , 'filterDir'   => $cfg->{'filterResultsDir'}
				     , 'mappingDir'  => $cfg->{'mappingResultsDir'}
				     , 'alignerPath' => getBinPath($cfg, $cfg->{'aligner'})
				     , 'alignerExeWithoutPath' => $cfg->{'aligner'}
				     , 'samtools'    => getBinPath($cfg, "samtools")
				     , 'gtfFile'     => $cfg->{'gtfFile'}                    # Optional, only used by Tophat (probably). Can be "NA" if there is not one
				     , 'minMapQ'     => $cfg->{'minMapQ'}
				     , 'bowtie2Index'=> $cfg->{'bowtie2Index'}
				     , 'bwaIndex'    => $cfg->{'bwaIndex'}
				     , 'starIndexDir'=> $cfg->{'starIndexDir'}
				     , 'sampleName'  => $k
				     , 'inputFile1'  =>               $sKeys[0]             
				     , 'inputFile2'  => ($isPaired) ? $sKeys[1] : $NA_VALUE  # <== important to still pass $NA_VALUE if it's unpaired
				     , 'sortMethod'  => "BY_COORD"                           # check 04.mapping.pl for valid values; should be "BY_COORD" or "BY_NAME" only!
				     , 'finalBam'    => getBamFullPath($cfg, $k) #catfile($cfg->{'mappingResultsDir'}, "${k}_${genomeShortName}_q${mapQ}.bam")
				   };
			my ($jn,$qs) = setJobInfo($cfg, "$STEP_MAPPING", "${k}", $dep, $vars, $cfg->{monkeyPoo}, "04.mapping.pl");
			$remember{$REM_ALIGN_BAM}{$k} = {"job"=>$jn, "id"=>$k, "bam"=>$vars->{'finalBam'}, "isPaired"=>$isPaired, "category"=>$sampleCategory, "replicate"=>$sampleRep }; # remember that we generate this BAM file, and remember some of the metadata that we will need later!
			$remember{$REM_EXPERIMENT_CATEGORY}{$sampleCategory}{$sampleRep} = { "bam"=>$vars->{'finalBam'} }; # Remember this aligned bam file path, and that it was in this experimental category! We'll use this information for Cuffdiff.
		}
	}

	if (isSplicedAligner($cfg->{aligner}) && $cfg->{browserBins}) {
		# If you're running a spliced aligner and ALSO browser bins, we output a warning that the bins do not always
		# reflect actual splicing behavior. Note that this is NOT AN ERROR, but is actually intended behavior.
		($GLOBAL_VERBOSE) && configWarning("(No line number)", $cfg->{filename}, "Note that the 'browserBins' tracks are not 100% accurate visual representations of read behavior in a spliced alignment. For example, in the event that a read spans three exons, it will only be displayed at its location in the first exon. This does not generate any errors when tallying up gene-level counts, but may be slightly misleading if you expect the browser bins to always correspond exactly to the read coverage! If you want a base-level view of all your reads, check the 'browserWigs' wiggle (actually '.bigWig' / '.bw') tracks.");
	}
	
	writeup($cfg, qq{Alignment of the provided samples to the reference genome was performed using:\n}
		. qq{     a) For SPLICED reads: Tophat 2.0.13 [CITE_TOP]\n}
		. qq{     b) For UNSPLICED reads: Bowtie 2.2.4 [CITE_BOW]\n}
		. qq{     Exact additional alignment parameters:\n}
		. qq{          --min-anchor=5 --segment-length=25 --no-coverage-search --segment-mismatches=2 --splice-mismatches=2 \n}
		,{"CITE_TOP"=>qq{TopHat: Kim D, Pertea G, Trapnell C, Pimentel H, Kelley R, Salzberg SL. "TopHat 2: accurate alignment of transcriptomes in the presence of insertions, deletions and gene fusions." Genome Biology 2011, 14:R36.}
		  , "CITE_BOW" =>qq{Bowtie: Langmead B, Salzberg S. "Fast gapped-read alignment with Bowtie 2." Nature Methods. 2012, 9:357-359.}});
	
	# ========= [START] SANITY CHECK THE INPUT CATEGORY (SAMPLE GROUP) NAMES: MAKE SURE THERE AREN'T ANY ACCIDENTAL CHANGES IN CAPITALIZATION =========
	my %seenCategories = (); # hash: key = lower-case version of the name, value = the original-case version of the name
	foreach my $ccc (keys(%{$remember{$REM_EXPERIMENT_CATEGORY}})) {
		if (exists($seenCategories{lc($ccc)}) && $seenCategories{lc($ccc)} ne $ccc) {
			confess "[FATAL CAPITALIZATION ERROR in config file]: You had two experimental categories with the same name but DIFFERENT capitalization! This is NOT ALLOWED, because it's probably an error. Probably the capitalization of <$ccc> and <" . $seenCategories{lc($ccc)} . "> should be the same. Fix this and re-submit your job.";
		}
		$seenCategories{lc($ccc)} = $ccc; # Save the category name, but note that the KEY is the LOWER-CASED version of the name. This lets us detect wrong-capitalization versions of the otherwise-same experimental groups.
	}
	# ========= [DONE] SANITY CHECK THE INPUT CATEGORY (SAMPLE GROUP) NAMES =========
	# Tags
	if ($cfg->{doDensity}) {
		foreach my $h ($sampleHash,$inputHash) {
			next if (!defined($h)); # Skip any undefined 'inputHash' cases.
			foreach my $k (keys %$h) {
			        my @sKeys    = sort(keys(%{$h->{$k}}));
			        my $isPaired = (2==scalar(@sKeys)) ? 1 : 0;
				my $inFile   = getBamFullPath($cfg, $k); #catfile($cfg->{'mappingResultsDir'}, "${k}_${genomeShortName}_q${mapQ}.bam")
				  #my $inFile   = catfile($cfg->{'mappingResultsDir'},"${k}_$cfg->{genomeName}_q$cfg->{minMapQ}.bam");
				my $dep  = { 'mappingJob' => $cfg->{jobs}->{$STEP_MAPPING}->{$k}->{jobName} };
				my $vars = { 'force' => $FORCE_RERUN, 'verbose' => $GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG };
				$vars->{genome} = $cfg->{genomeName};
				$vars->{seqType} = $cfg->{details}->{$k}->{type};
				$vars->{bam2bed}  = getBinPath($cfg, "bam2bed");
				$vars->{samtools} = getBinPath($cfg, "samtools");
				$vars->{sortBed}  = getBinPath($cfg, "sortBed");
				$vars->{tagsDir} = $cfg->{tagResultsDir};
				$vars->{inputFile} = $inFile;
				$vars->{sampleName} = $k;
				$vars->{isPaired} = $isPaired;
				setJobInfo($cfg, "$STEP_TAGS", "${k}",$dep,$vars,$cfg->{monkeyPoo},"05.tags.pl");
			}
		}
		writeup($cfg, "Tags & density calculation", {});
	}

	# Tag mapping stats
	if ($cfg->{doDensity}) {
		my $dep;
		foreach my $k (sort keys(%{$cfg->{jobs}->{$STEP_TAGS}})) {
			$dep->{$k} = $cfg->{jobs}->{$STEP_TAGS}->{$k}->{jobName};
		}
		my $vars = { 'force' => $FORCE_RERUN, 'verbose' => $GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG };
		$vars->{genome}     = $cfg->{genomeName};
		$vars->{minMapQ}    = $cfg->{minMapQ};
		$vars->{tagsDir}    = $cfg->{tagResultsDir};
		$vars->{mappingDir} = $cfg->{mappingResultsDir};
		setJobInfo($cfg, $STEP_TAG_STATS, "all",$dep,$vars,$cfg->{monkeyPoo},"05.xSummary.pl");
		writeup($cfg, "Tags mapping stats (summary)", {});
	}
	
	## ChIP/exo specific processing: Tag density
	if ($cfg->{doDensity}) {
		foreach my $k (keys %$sampleHash) {
			my $hasMatchingInput  = sampleHasMatchingInput($cfg->{details}->{$k});
			my $matchingInputName = getMatchingInput($cfg->{details}->{$k}); # can be NA if there is no matching input
			#print "The matching input ($hasMatchingInput), since the config was <" . $cfg->{details}->{$k}->{input} . "> for sample $k\nwas this input: $matchingInputName\n";
			my $dep = { 'sampleTagsJob' =>                       $cfg->{jobs}->{$STEP_TAGS}->{$k}->{jobName} ,
				    'inputTagsJob'  => ($hasMatchingInput) ? $cfg->{jobs}->{$STEP_TAGS}->{$matchingInputName}->{jobName} : $NA_VALUE }; # Depends on the SAMPLE tag job AND (if it exists) the INPUT tag job
			my $vars = { 'force' => $FORCE_RERUN, 'verbose' => $GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG };
			$vars->{genome}     = $cfg->{genomeName};
			$vars->{seqType}    = $cfg->{details}->{$k}->{type};
			$vars->{bedmap}     = getBinPath($cfg, "bedmap");
			$vars->{binsFile}   = $cfg->{genomicBins};
			$vars->{binsFileID} = $cfg->{genomicBinsID};
			$vars->{tagsDir}    = $cfg->{tagResultsDir};
			$vars->{densityDir} = $cfg->{densityResultsDir};
			$vars->{sampleName} = $k;
			$vars->{inputName}  = $matchingInputName; # <== allowed to be 'NA'
			setJobInfo($cfg, $STEP_DENSITY, "${k}", $dep, $vars, $cfg->{monkeyPoo}, "06.density.pl");
		}
		writeup($cfg, "Tag density", {});
	}

	#print "do epaks is: $cfg->{globalDoPeaks}\n";
	#print "do dens is: $cfg->{doDensity}\n";
	#print "d qc is: $cfg->{doQC}\n";
	if ($cfg->{globalDoPeaks}) {
		# Generate the peaks for each sample. Depends on 05 and 06
	        foreach my $k (keys %$sampleHash) {
		    if ($cfg->{details}->{$k}->{needsPeakCall}) {
			my $dep = { "densityJob" => $cfg->{jobs}->{$STEP_DENSITY}->{$k}->{jobName} }; # Each of these only depends on the ONE density job.
			my $vars = { 'force' => $FORCE_RERUN, 'verbose' => $GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG };
			$vars->{genome}     = $cfg->{genomeName};
			$vars->{seqType}    = $cfg->{details}->{$k}->{type};
			$vars->{bedmap}     = getBinPath($cfg, "bedmap");
			$vars->{bedops}     = getBinPath($cfg, "bedops");
			$vars->{sortBed}    = getBinPath($cfg, "sortBed");
			$vars->{tagsDir}    = $cfg->{tagResultsDir};
			$vars->{peaksDir}   = $cfg->{peakResultsDir};
			$vars->{densityDir} = $cfg->{densityResultsDir};
			$vars->{sampleName} = $k;
			setJobInfo($cfg, $STEP_PEAKS, "${k}", $dep, $vars, $cfg->{monkeyPoo}, "07.peaks.pl");
		    }
		}

		my $pooData = catfile($cfg->{monkeyPoo}, $poo_data_subdir);
		# Generate the gem peaks for each sample. Depends on 04
		foreach my $k (keys %$sampleHash) {
		    if ($cfg->{details}->{$k}->{needsGemPeaks}) {
			my $readDist          = undef;
			if ($cfg->{details}->{$k}->{needsExo}) { $readDist = catfile($pooData, $GEM_READ_DIST_CHIP_EXO); }
			else {	                                 $readDist = catfile($pooData, $GEM_READ_DIST_DEFAULT); }
			my $genomeShortName   = $cfg->{genomeName};
			my $hasMatchingInput  = sampleHasMatchingInput($cfg->{details}->{$k});
			my $matchingInputName = getMatchingInput($cfg->{details}->{$k}); # can be NA if there is no matching input
			my $dep = { 'sampleTagsJob' =>                       $cfg->{jobs}->{$STEP_MAPPING}->{$k}->{jobName} ,
				    'inputTagsJob'  => ($hasMatchingInput) ? $cfg->{jobs}->{$STEP_MAPPING}->{$matchingInputName}->{jobName} : $NA_VALUE }; # Depends on the SAMPLE tag job AND (if it exists) the INPUT tag job
			my $vars = { 'force' => $FORCE_RERUN, 'verbose' => $GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG };
			#$vars->{pbs_mem}         = $QSETTINGS{gem}{pbs_mem}
			#$vars->{pbs_ncpus}       = $QSETTINGS{gem}{pbs_ncpus}
			$vars->{genome}          = $genomeShortName;
			$vars->{seqType}         = $cfg->{details}->{$k}->{type};
			$vars->{sampleName}      = $k;
			$vars->{gemJar}          = getBinPath($cfg, "gem"); # this is really gem.jar, not a "script" per se
			$vars->{chrSizes}        = $cfg->{chromSizesFile};
			$vars->{genomeFastaPath} = $cfg->{genomeFastaDir}; # <== genomeFastaDir is the directory with PER-CHROMOSOME fasta files (e.g., it isn't "hg19.fa", it's "chr1.fa" or something)
			$vars->{peaksDir}        = $cfg->{gemPeakResultsDir};
			$vars->{expBam}          = getBamFullPath($cfg, $k); #catfile($cfg->{'mappingResultsDir'}, ("${k}_${genomeShortName}_q" . ${cfg}->{minMapQ} . ".bam"));
			$vars->{ctrlBam}         = ($hasMatchingInput) ? getBamFullPath($cfg, ${matchingInputName}) : $NA_VALUE; #($hasMatchingInput) ? catfile($cfg->{'mappingResultsDir'}, ("${matchingInputName}_${genomeShortName}_q" . ${cfg}->{minMapQ} . ".bam")) : $NA_VALUE;
			$vars->{readDist}        = $readDist;
			setJobInfo($cfg, $STEP_GEM_PEAKS, "${k}", $dep, $vars, $cfg->{monkeyPoo}, "07a.peaksGem.pl");
		    }
		}

		# Generate bcp peaks for each sample.
		foreach my $k (keys %$sampleHash) {
			my $hasMatchingInput  = sampleHasMatchingInput($cfg->{details}->{$k});
			my $matchingInputName = getMatchingInput($cfg->{details}->{$k}); # can be NA if there is no matching input
			if ($cfg->{details}->{$k}->{needsBcpPeaks} and $hasMatchingInput) {
				my $genomeAbbrev = $cfg->{genomeName};
				my $dep = { 'sampleTagsJob' => $cfg->{jobs}->{$STEP_TAGS}->{$k}->{jobName} ,
					    'inputTagsJob'  => $cfg->{jobs}->{$STEP_TAGS}->{$matchingInputName}->{jobName} };
				my $vars = { 'force' => $FORCE_RERUN, 'verbose' => $GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG };
				$vars->{bcp}         = getBinPath($cfg, "bcp");
				$vars->{genome}      = $genomeAbbrev;
				$vars->{seqType}     = $cfg->{details}->{$k}->{type};
				$vars->{sampleName}  = $k;
				$vars->{peaksDir}    = $cfg->{bcpPeakResultsDir};
                                $vars->{chromSizes}  = $cfg->{chromSizesFile};
				$vars->{expTagsBed} = catfile($cfg->{'tagResultsDir'}, ("${k}_${genomeAbbrev}_tags.bed"));
				$vars->{ctlTagsBed} = catfile($cfg->{'tagResultsDir'}, ("${matchingInputName}_${genomeAbbrev}_tags.bed"));
				setJobInfo($cfg, $STEP_BCP_PEAKS, "${k}", $dep,$vars,$cfg->{monkeyPoo},"07b.peaksBCP.pl");
			}
		}

		# Generate atacSeq nucleosome positioning data if necessary.
		# This uses the BAM files to figure out where, in ATAC-seq, the nucleosomes are probably located.
		# (Note that this is a different question from "which chromatin is *accessible*".)
		foreach my $k (keys %$sampleHash) {
			if ($cfg->{details}->{$k}->{needsATAC}) {
				# ATAC-seq data: for each ATAC-seq sample, find out where the accessible chromatin is and where the nucleosomes are.
				my $dep = { 'peakCalling' => $cfg->{jobs}->{$STEP_PEAKS}->{$k}->{jobName} }; # Depends on peak-calling
				my $vars = { 'force'=>$FORCE_RERUN, 'verbose'=>$GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG
					     , 'genome'      => $cfg->{genomeName}
					     , 'seqType'     => $cfg->{details}->{$k}->{type}
					     , 'sampleName'  => $k
					     , 'atacDir'     => $cfg->{atacResultsDir}
					     , 'peaksDir'    => $cfg->{peakResultsDir}
					     , 'mappingDir'  => $cfg->{mappingResultsDir}
					     , 'genomeFasta' => $cfg->{genomeFasta}
					     , 'atacPython'  => getBinPath($cfg, "atacPython")
					     , 'nucleoatac'  => getBinPath($cfg, "nucleoatac")
					     , 'bedops'      => getBinPath($cfg, "bedops")
					   };
				warn "NOTE: THIS NEVER ACTUALLY EXECUTES ANYTHING! THe nucleo-atac job NEVER RUNS.";
				warn "Message about this from (presumably) Sean: # this is suspended pending evaluation.  --AlexW";
				# setJobInfo($cfg, $STEP_NUCLEOATAC, "${k}", $dep,$vars,$cfg->{monkeyPoo},"07c.NucleoAtac.pl");
			}
		}
		
		# Now generate the peak SUMMARY stats...
		my $dep;
		foreach my $k (sort keys %{$cfg->{jobs}->{$STEP_PEAKS}}) {
			$dep->{$k} = $cfg->{jobs}->{$STEP_PEAKS}->{$k}->{jobName}; # depends on all the peaks being done
		}
		my $vars = { 'force' => $FORCE_RERUN, 'verbose' => $GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG
			     ,'genome'  => $cfg->{genomeName}
			     ,'bedmap'  => getBinPath($cfg, "bedmap")
			     ,'tagsDir' => $cfg->{tagResultsDir}
			     ,'peaksDir'=> $cfg->{peakResultsDir}  };
		setJobInfo($cfg, $STEP_PEAK_SUMMARY, "summary", $dep,$vars,$cfg->{monkeyPoo},"07.xSummary.pl");
		writeup($cfg, "Peak calling?", {});
	}

	## Motif discovery module (only occurs if there are chip/exo data, and only on chip/exo data)
	if ($cfg->{globalDoPeaks}) {
		my $numSamplesForMotifFinding = 0;
		foreach my $k (keys %$sampleHash) {
			if ($cfg->{details}->{$k}->{type} ne "rna") { # Alex's note: note that we do NOT do motif discovery for type 'rna'! Just 'chip' and 'exo'.
				$numSamplesForMotifFinding++;
				my $infile = "";
				if ($cfg->{'details'}->{$k}->{'type'} eq "exo") {
					$infile = catfile($cfg->{'peakResultsDir'},"${k}_$cfg->{genomeName}_footprints.bed");
				} else {
					$infile = catfile($cfg->{'peakResultsDir'},"${k}_$cfg->{genomeName}_peaks.bed");
				}
				my $dep  = { "peaksJob"=>$cfg->{jobs}->{$STEP_PEAKS}->{$k}->{jobName} }; # Depends on "07a_peaks"
				my $vars = { 'force'=>$FORCE_RERUN, 'verbose'=>$GLOBAL_VERBOSE, 'debug'=>${GLOBAL_DEBUG}
					     , 'genome' => $cfg->{'genomeName'}
					     , 'bedops'     => getBinPath($cfg, "bedops")
					     , 'bedtools'   => getBinPath($cfg, "bedtools")
					     , 'peakMotifs' => getBinPath($cfg, "peakMotifs")
					     , 'fasta'      => $cfg->{'genomeFasta'}
					     , 'repMask'    => $cfg->{'repeatMask'}
					     , 'motifsDir'  => $cfg->{'motifDiscDir'}
					     , 'infile'     => ${infile}
					     , 'ZLM_LICENSE' => $ZLM_LICENSE_HARD_CODED
					     , 'sampleName' => $k };
				setJobInfo($cfg, "${STEP_MOTIFS}", "${k}", $dep,$vars,$cfg->{monkeyPoo},"10.motifDiscovery.pl");
			}
		}
		if ($numSamplesForMotifFinding > 0) { # summarize motif discovery, put all discovered motifs into a single transfac format file
			my $dep;
			foreach my $k (sort keys %{$cfg->{jobs}->{$STEP_MOTIFS}}) {
				$dep->{$k} = $cfg->{jobs}->{$STEP_MOTIFS}->{$k}->{jobName}; # Depend on ALL motif jobs being done!
			}
			my $vars;
			$vars->{force}       = $FORCE_RERUN;
			$vars->{studyName}   = $cfg->{studyName};
			$vars->{motifDir}    = $cfg->{motifDiscDir};
			$vars->{peakDir}     = $cfg->{peakResultsDir};
			$vars->{bedmap}      = getBinPath($cfg, "bedmap");
			$vars->{bedops}      = getBinPath($cfg, "bedops");
			$vars->{bedtools}    = getBinPath($cfg, "bedtools");
			$vars->{genomeFasta} = $cfg->{genomeFasta};
			$vars->{step1}       = getBinPath($cfg, "step1Motif");
			$vars->{step2}       = getBinPath($cfg, "step2Motif");
			$vars->{sortBed}     = getBinPath($cfg, "sortBed");
			setJobInfo($cfg, "${STEP_MOTIF_SUMMARY}", "all",$dep,$vars,$cfg->{monkeyPoo},"10.xSummary.pl");
			writeup($cfg, "Motifs?", {});
		}
	}
	## end motif discovery module
    
	if ($cfg->{browserWigs}) { # Alex's UCSC files -- generated if "browserWigs" is specified
		my $gExe = catfile($cfg->{bin}->{basedir}, "genomeCoverageBed"); (-x $gExe) or confess "[ERROR]: FAILED TO FIND the required 'genomeCoverageBed' executable in path '$gExe'. Check to make sure this file exists **and** is executable by the current user!";
		my $wExe = catfile($cfg->{bin}->{basedir}, "wigToBigWig");       (-x $wExe) or confess "[ERROR]: FAILED TO FIND the required 'wigToBigWig' executable in path '$wExe'. Check to make sure this file exists **and** is exectuable by the current user!";
		my $dep                = getHashPointerOfAllBamJobs(); # hash pointer -- require ALL bam alignment jobs are DONE!
		my $bamMultiFileString = join(" ", getArrayOfAlignedBamFilenames("all")); # From the "remember" hash
		my $vars = {   'force'=>$FORCE_RERUN, 'verbose'=>$GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG
			       , 'pbs_ncpus' => 1, 'pbs_mem' => "8g", 'pbs_walltime'=> "172808"
			       , 'outputDir'                 =>  $cfg->{browserWigResultsDir} # the monkey output directory. NOT the http-accessible one! That is set by "scpDestination"
			       , 'convertToGenomeBrowserExe' =>  getBinPath($cfg, "convertToGenomeBrowserExe") # executable location
			       , 'inputMultiFileString'      =>  ${bamMultiFileString}
			       , 'binDir'                    =>  $cfg->{bin}->{basedir}
			       , 'serverBaseURL'             =>  $cfg->{tracksURL}
			       , 'tracksShowEveryRead'       =>  $cfg->{tracksShowEveryRead}
			       , 'rsyncDestinationBase'      =>  $cfg->{tracksRsync} # <== copy files to THIS destination, the filesystem path that corresponds to the 'serverBaseURL' above.
			       , 'studyName'                 =>  $cfg->{studyName}     };
		setJobInfo($cfg, $STEP_BROWSER_WIGGLE, "all", $dep,$vars,$cfg->{monkeyPoo},"48.browser_agw.pl");
		writeup($cfg, "Generating UCSC-compatible browser (wiggle) tracks", {});
	}
	
	if ($cfg->{browserBins}) { # Sean's binned UCSC files -- generated if "browserBins" is specified
		($cfg->{doDensity}) or confess "Bug--this requires density to be calculated beforehand.";

		foreach my $k (keys %$sampleHash) {
			my $hasMatchingInput  = sampleHasMatchingInput($cfg->{details}->{$k});
			my $matchingInputName = getMatchingInput($cfg->{details}->{$k}); # can be NA if there is no matching input
			my $dep = { "densityJob" => $cfg->{jobs}->{$STEP_DENSITY}->{$k}->{jobName} }; # Each of these only depends on the ONE density job.

			($cfg->{globalDoPeaks} || !$cfg->{details}->{$k}->{needsPeakCall}) or confess "Bug--this requires peaks to be calculated beforehand.";

			if ($cfg->{details}->{$k}->{needsPeakCall})                       {    $dep->{peaksJob} = $cfg->{jobs}->{$STEP_PEAKS}->{$k}->{jobName};      } # basic peak calling dependencies
			if ($cfg->{details}->{$k}->{needsGemPeaks})                       { $dep->{gemPeaksJob} = $cfg->{jobs}->{$STEP_GEM_PEAKS}->{$k}->{jobName};  } # gem peak calling dependencies
			if ($cfg->{details}->{$k}->{needsBcpPeaks} and $hasMatchingInput) { $dep->{bcpPeaksJob} = $cfg->{jobs}->{$STEP_BCP_PEAKS}->{$k}->{jobName};  } # bcp peak calling dependencies
			if ($cfg->{details}->{$k}->{needsATAC}) { # atac (chromatin accessibility) dependencies
			    # this is commented out for evaluation
			    # $dep->{atacSeqJob} = $cfg->{jobs}->{$STEP_NUCLEOATAC}->{$k}->{jobName};
			}
			my $vars = { 'force' => $FORCE_RERUN, 'verbose' => $GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG
				     , 'pbs_ncpus' => 1, 'pbs_mem' => "8g", 'pbs_walltime'=> "172808"
				     , 'trackFileSuffix' => $TRACKFILE_SUFFIX
				     , 'genome'      => $cfg->{genomeName}
				     , 'seqType'     => $cfg->{details}->{$k}->{type} # 'chip', 'rna', 'exo'... etc...
				     , 'bigWig'      => getBinPath($cfg, "bigWig")
				     , 'bigBed'      => getBinPath($cfg, "bigBed")
				     , 'peaksDir'    => $cfg->{peakResultsDir}
				     , 'gemPeaksDir' => $cfg->{gemPeakResultsDir}
				     , 'bcpPeaksDir' => $cfg->{bcpPeakResultsDir}
				     , 'atacSeqDir'  => $cfg->{atacResultsDir}
				     , 'densityDir'  => $cfg->{densityResultsDir}
				     , 'browserDir'  => $cfg->{browserBinResultsDir}
				     , 'studyName'   => $cfg->{studyName}
				     , 'sampleName'  => ${k}
				     , 'inputName'   => $matchingInputName # just to check if BCP applicable
				     , 'chromSizes'  => $cfg->{chromSizesFile}
				     , 'genomicBins' => $cfg->{genomicBins}
				     , 'url'         => $cfg->{tracksURL}                     }; # END definition of '$vars'
			setJobInfo($cfg, $STEP_BROWSER_BINS, "${k}", $dep, $vars, $cfg->{monkeyPoo}, "08.browser.pl");
			print "Set up this dependency:";
			print $dep->{densityJob} . "\n";
		}
		
		# After the browser tracks have all been added to the queue, generate the SUMMARY...
		my $dep;	# hash reference
		foreach my $k (sort keys %{$cfg->{jobs}->{$STEP_BROWSER_BINS}}) {
			$dep->{$k} = $cfg->{jobs}->{$STEP_BROWSER_BINS}->{$k}->{jobName}; # Add dependencies
		}
		my $vars = {'force'=>$FORCE_RERUN, 'verbose'=>$GLOBAL_VERBOSE, 'debug'=>${GLOBAL_DEBUG}
			    , 'pbs_ncpus' => 1, 'pbs_mem' => "8g", 'pbs_walltime'=> "86408"
			    , 'trackFileSuffix' => $TRACKFILE_SUFFIX
			    , 'rsyncDir'   => $cfg->{'tracksRsync'} . "/" . $cfg->{'studyName'}
			    , 'browserDir' => $cfg->{'browserBinResultsDir'} };
		setJobInfo($cfg, $STEP_BROWSER_BIN_SUMMARY, "all",$dep,$vars,$cfg->{monkeyPoo},"08.xSummary.pl");

		writeup($cfg, "Browser bins", {});
	}

	# Generate windows
	if ($cfg->{doWindows}) {
	    foreach my $k (keys %$sampleHash) { # Note: does NOT go through the 'inputhash'! Inputs don't need windows I guess.
		    # (Sean's note: it doesn't go through input because the samples have already been input normalized...)
		    my $dep  = { 'densityJob' => $cfg->{jobs}->{$STEP_DENSITY}->{$k}->{jobName} };
		    my $vars = { 'force' => $FORCE_RERUN, 'verbose' => $GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG
				 ,'pbs_ncpus' => 1, 'pbs_mem' => "8g", 'pbs_walltime'=> "172821"
			       };
		    $vars->{genome}      = $cfg->{genomeName};
		    $vars->{bedmap}      = getBinPath($cfg, "bedmap");
		    $vars->{sortBed}     = getBinPath($cfg, "sortBed");
		    $vars->{densityDir}  = $cfg->{densityResultsDir};
		    $vars->{windowDir}   = $cfg->{windowResultsDir};
		    $vars->{sampleName}  = $k;
		    $vars->{geneWindows} = $cfg->{geneWindows};
		    $vars->{analyze}     = $cfg->{r}->{window_figures};
		    setJobInfo($cfg, $STEP_WINDOWS, "${k}",$dep,$vars,$cfg->{monkeyPoo},"09.windows.pl");
	    }
	    writeup($cfg, "Windows", {});
        }

	my $isSubreadBeingDone = 0; # <== many downstream scripts shouldn't run if subread_featureCounts isn't being run.
	# ================================= [START] RUN SUBREAD FEATURE COUNTS ON ALL MAPPED BAM FILES **TOGETHER** ========================
	if (isNA($cfg->{'gtfFile'})) {
		verboseSkipPrint("[OK] [NOTE] We are skipping the 'subread featureCount'-tallying step since no GTF file was specified (the 'gtfFile' was set to 'NA'). This means that many downstream programs (edgeR, etc...) also won't be run, since they depend on featureCounts having run properly.\n");
	} else {
		# NB: "subread" is the name of a program (i.e. a "sub-read" is not actually a type of read)
		# Prerequisites: *all* aligned BAM files (in step 04a)
		# Note: runs single-end and paired-end files SEPARATELY in the unusual (but possible) case that both exist!
		# We then wrote a script to generates a single merged output file from the (possibly two) subread runs.
		my $dep = getHashPointerOfAllBamJobs();
		my $seBams = " " . join(" ", getArrayOfAlignedBamFilenames("single_end")); # All SINGLE END bam files with spaces between each one. Or " " if there aren't any! Makes sure there is at least SOME character here; qsub cannot submit the TOTALLY empty (zero length) string for some reason!
		my $peBams = " " . join(" ", getArrayOfAlignedBamFilenames("paired_end")); # All PAIRED END bam files with spaces between each one. Or " " if there aren't any! Makes sure there is at least SOME character here; qsub cannot submit the TOTALLY empty (zero length) string for some reason!
		my $vars = { 'force'=>$FORCE_RERUN, 'verbose'=>$GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG,
			     'pbs_ncpus' => 6, 'pbs_mem' => "12g", 'pbs_walltime'=> "57620",
			     'subreadFeatureCountsExe'     => getBinPath($cfg, "subreadFeatureCounts") ,
			     'subreadCountsFilename'       => $cfg->{'subread'}->{'countsFile'}   ,
			     'outputDir'                   => $cfg->{'subreadCountsDir'}          ,
			     'gtf'                         => $cfg->{gtfFile}                     , # Mandatory
			     'singleEndBamFilesSpaceDelim' => $seBams                             , # <== Ok if this is " " (one space)
			     'pairedEndBamFilesSpaceDelim' => $peBams }; # <== Ok if this is " " (one space)
		setJobInfo($cfg, "$STEP_SUBREAD", "all", $dep,$vars,$cfg->{monkeyPoo},"20a.count.subread.pl");
		$isSubreadBeingDone = 1; # Ok, apparently we ARE running subread for this job! Requires a GTF file.
		
		writeup($cfg,
			qq{Reads were assigned to genes using "featureCounts" [CITE_FEATURE], part of the Subread suite (http://subread.sourceforge.net/)\n}
			. qq{     * Gene-level counts were arrived at using Ensembl's gene annotation, in GTF format.}
			, {"CITE_FEATURE"=>qq{FeatureCounts: Liao Y, Smyth GK and Shi W. "FeatureCounts: an efficient general-purpose program for assigning sequence reads to genomic features." Bioinformatics, 30(7):923-30, 2014.}});
	}
	# ================================= [DONE] RUNNING SUBREAD ========================
    
	if ($cfg->{doExpression} && !isNA($cfg->{'gtfFile'})) {	
		# ================================= [START] RUN CUFFDIFF ON ALL MAPPED BAM FILES **TOGETHER** ========================
		# Depends on: *all* BAM files being aligned (in step 04a), and nothing else.
		my ($labelsCommaDelim, $bamsByGroup) = getGroupLabelsAndBamStringByGroup("<COMMA>", " ", "<COMMA>"); # <== Will eventually look something like "exp1.bam,exp2.bam  ctrl1.bam,ctrl2.bam,ctrl3.bam"
		my $dep = getHashPointerOfAllBamJobs();
		my $vars = { 'force'=>$FORCE_RERUN, 'verbose'=>$GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG,
			     'pbs_ncpus' => 8, 'pbs_mem' => "12g", 'pbs_walltime'=> "172821", # cuffdiff takes forever to run
			     'cuffdiffExe'    => getBinPath($cfg, "cuffdiff") ,
			     'outputDir'      => $cfg->{cuffdiffDir}     , # <== the output directory, which will be created
			     'gtfFile'        => $cfg->{gtfFile}         ,
			     'genomeFasta'    => $cfg->{genomeFasta}     ,
			     'labels'         => $labelsCommaDelim       , # <== labels is delimited by "<COMMA>" and spaces, NOT the actual comma ',', because of qsub weirdness
			     'bamFilesByGroup'=> $bamsByGroup # <== bamFilesByGroup is delimited by "<COMMA>" and spaces, NOT the actual comma ',', because of qsub weirdness
			   };
		setJobInfo($cfg, "$STEP_CUFFDIFF", "all",$dep,$vars,$cfg->{monkeyPoo},"22.cuffdiff.pl");

		writeup($cfg, "In order to generate certain cummeRbund figures, we run CuffDiff---however, note that we are NOT using CuffDiff for calculating differential expression!");
		# ================================= [DONE] RUN CUFFDIFF ON ALL MAPPED BAM FILES **TOGETHER** ========================

		# ================================= [START] NOW RUN CUMMERBUND ON THE CUFFDIFF OUTPUT FILES **TOGETHER** ========================
		# Step added Feb 10, 2015. Requires R and the 'cummeRbund' package
		my $depCummerbund  = { 'cuffdiff'=> $cfg->{jobs}->{$STEP_CUFFDIFF}->{"all"}->{jobName} };
		my $varsCummerbund = { 'force'=>$FORCE_RERUN, 'verbose'=>$GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG,
				       'pbs_ncpus' => 1, 'pbs_mem' => "12g", 'pbs_walltime'=> "64821",
				       'cuffdiffDir' => $cfg->{cuffdiffDir}     ,
				       'script'      => $cfg->{'r'}->{cummerbund} ,
				       'outputDir'   => $cfg->{cummerbundDir}   ,
				       'R_EXE'       => getBinPath($cfg, "R_EXE")                 
				     };
		setJobInfo($cfg, $STEP_CUMMERBUND, "all",$depCummerbund,$varsCummerbund,$cfg->{monkeyPoo},"23a.cummerbund.pl");
		writeup($cfg, qq{Certain additional figures (indicated with "Cuffdiff" in the title) were generated using Cuffdiff + cummeRbund [CITE_CUMMERBUND].}, {"CITE_CUMMERBUND"=>qq{CummeRbund: L. Goff, C. Trapnell and D. Kelley (2012). "cummeRbund: Analysis, exploration, manipulation, and visualization of Cufflinks high-throughput sequencing data." R package version 2.6.1.}});
		# ================================= [DONE] RUNNING CUMMERBUND ===================================================================
	} else {
		(isNA($cfg->{'gtfFile'})) and verboseSkipPrint("We are skipping the 'cuffdiff' step no GTF file was specified (the 'gtfFile' was set to 'NA').\n");
	}
    
	# ================================= [START] MAKE A PAIRS PLOT OF THE 'SUBREAD FEATURECOUNTS' 'CPM' VALUES ===========
	if (not $isSubreadBeingDone) {
		verboseSkipPrint("We are skipping the 'pairs plot' step since it depends on the 'subread' step being completed, which did not occur.\n");
	} else {
		# Only depends on "subread featureCounts" being completed.
		my $dep  = { 'subread' => $cfg->{jobs}->{$STEP_SUBREAD}->{"all"}->{jobName} };
		my $vars = { 'force'=>$FORCE_RERUN, 'verbose'=>$GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG
			     ,'pbs_ncpus' => 1, 'pbs_mem' => "4g", 'pbs_walltime'=> "86421",
			     ,'R_EXE'                  => getBinPath($cfg, "R_EXE")
			     ,'subreadCountsFilename'  => $cfg->{'subread'}->{'countsFile'}
			     ,'subreadDir'             => $cfg->{'subreadCountsDir'}  # <== the input directory, where we'll look for a file with a filename HARD-CODED in the R script
			     ,'outputDir'              => $cfg->{'subreadCountsDir'}  # <== (unusually, the output directory is the *same* as the input subreadDir)
			     ,'pairScript'             => $cfg->{'r'}->{'pairPlot'}
			   };
		setJobInfo($cfg, "$STEP_SUBREAD_SUMMARY_FIGS", "all", $dep,$vars,$cfg->{monkeyPoo},"20b.summary.figures.pl");
	}
	# ================================= [DONE] MADE A PAIRS PLOT OF THE 'SUBREAD FEATURECOUNTS' 'CPM' VALUES ===========

	# ================================= Differential expression calculation using edgeR ===============================
	my $isEdgeRBeingDone = 0;
	if (not $isSubreadBeingDone) {
		verboseSkipPrint("We are skipping the 'edgeR differential expression' step since it depends on the 'subread' step being completed, which did not occur.\n");
	} else {
		$isEdgeRBeingDone = 1;
		my $dep  = { 'subread' => $cfg->{jobs}->{$STEP_SUBREAD}->{"all"}->{jobName} }; # edgeR only depends on "subread featureCounts" being completed.
		my $MAJOR_DELIM = ":"; # colon is the "major" delimiter. Do not use a comma, as 'torque' (the scheduler) doesn't like commas for command line args!
		my $MINOR_DELIM = "~"; # tilde is the "minor" delimiter. Do not use a comma, as 'torque' (the scheduler) doesn't like commas for command line args!
		my ($groupsColonDelim, $bamsPerGroup) = getGroupLabelsAndBamStringByGroup($MAJOR_DELIM, $MAJOR_DELIM, $MINOR_DELIM); # <== Will eventually look something like "exp:ctrl" "exp1.bam:exp2.bam  ctrl1.bam:ctrl2.bam:ctrl3.bam"
		if (1 == 1) {
			# Error checking. Limiting the scope to prevent us from using these variables, which we should not use again!
			my @groupArr = split($MAJOR_DELIM, $groupsColonDelim); # Colon-delimited
			my @bamsGrouped = split($MAJOR_DELIM, $bamsPerGroup);
			(scalar(@groupArr) == scalar(@bamsGrouped)) or confess("[ERROR] The number of groups indicated (" . scalar(@groupArr) . ": from the string '$groupsColonDelim') was NOT the same as the number of bams-clustered-by-group (" . scalar(@bamsGrouped) . ": from the string '$bamsPerGroup'). This is probably a programming error.");
		}
		my $vars = { 'force'=>$FORCE_RERUN, 'verbose'=>$GLOBAL_VERBOSE, 'debug'=>${GLOBAL_DEBUG}
			     ,'pbs_ncpus' => 1, 'pbs_mem' => "12g", 'pbs_walltime'=> "172821",
			     ,'RSCRIPT_EXE'                 => getBinPath($cfg, "RSCRIPT_EXE")
			     ,'edgeRScript'           => $cfg->{'r'}->{'edgeR'}
			     ,'subreadCountsFilePath' => catfile($cfg->{subreadCountsDir}, $cfg->{subread}->{countsFile})
			     ,'outputDir'             => $cfg->{'edgeRDir'}
			     ,'groupList'             => $groupsColonDelim
			     ,'bamsPerGroup'          => $bamsPerGroup
			     ,'majorDelim'            => $MAJOR_DELIM
			     ,'minorDelim'            => $MINOR_DELIM };
		setJobInfo($cfg, "$STEP_EDGER_DIFF_EXPR", "all", $dep,$vars,$cfg->{monkeyPoo},"50a.edgeR.pl");

		writeup($cfg, qq{We calculated differential expression P-values using edgeR [CITE_EDGER], an R package available through Bioconductor.\n}
			. qq{    * We first filter out any genes where there were not at least two samples with a CPM (counts per million) between 0.5 and 5000.}
			. qq{* CPM below 0.5 indicates non-detectable gene expression}
			. qq{* CPM above 5000 is typically only seen in mitochondrial genes (if these high-expression genes were not excluded, their counts would disproportionately affect the normalization).}
			. qq{* We then exclude these genes and re-normalize the remaining ones using "calcNormFactors" in edgeR.}
			. qq{* As a result, the sum of CPM in each sample is always exactly 1,000,000.}
			. qq{* Genes that would have had too-high or too-low CPM are reported as "NA".}
			. qq{* edgeR performs the calculation of P-values for the differential expression between samples, using a negative binomial distribution for gene expression.}
			. qq{* Finally we use the built-in R function "p.adjust" to calculating the FDR (false discovery rate) for each P-value, using the Benjamini-Hochberg method [CITE_FDR].}
			, {"CITE_EDGER"=>qq{EdgeR: Robinson MD, McCarthy DJ and Smyth GK. "edgeR: a Bioconductor package for differential expression analysis of digital gene expression data." Bioinformatics (2010), 26, 139-140.}, "CITE_FDR"=>qq{FDR: Benjamini Y, Hochberg Y. "Controlling the false discovery rate: a practical and powerful approach to multiple testing." Journal of the Royal Statistical Society Series B (1995), 57, 289-300.}});
	}
	# ================================= [DONE] Using edgeR ===========    



	# ================================= [START] MAKE A PAIRS PLOT OF THE 'SUBREAD FEATURECOUNTS' 'CPM' VALUES ===========
	if (not $isSubreadBeingDone) {
		verboseSkipPrint("We are skipping the 'AltAnalyze' step since it depends on the 'subread' step, which did not occur.\n");
	} else { # AltAnalyze currently only depends on "subread featureCounts" being completed. Not EdgeR.

			my $dep  = { 'subread'=> $cfg->{jobs}->{$STEP_SUBREAD}->{"all"}->{jobName} }; # <== depends on SUBREAD
			my $vars = { 'force'=>$FORCE_RERUN, 'verbose'=>$GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG
				     ,'python2_exe'                => getBinPath($cfg, "python2_7_exe")
				     ,'alt_analyze_py'             => getBinPath($cfg, "alt_analyze_py")
				     ,'featureCountsToPlainMatrix' => getBinPath($cfg, "featureCountsToPlainMatrix")
				     ,'subreadCountsFullPath'      => catfile($cfg->{subreadCountsDir}, $cfg->{subread}->{countsFile})
				     ,'species'                    => (defined($cfg->{species}) ? $cfg->{species} : $NA_VALUE)
				     ,'outputDir'                  => $cfg->{altAnalyzeDir}
				     ,'PYTHONPATH'                 => getMostExtensivePythonPath()
				   };
		setJobInfo($cfg, "$STEP_ALT_ANALYZE", "all", $dep,$vars,$cfg->{monkeyPoo},"52.alt_analyze.pl");
	}
	# ================================= [DONE] MADE A PAIRS PLOT OF THE 'SUBREAD FEATURECOUNTS' 'CPM' VALUES ===========

	# ================================= Post-edgeR processing with various scripts ===============================
	if (not $isEdgeRBeingDone) {
		verboseOkPrint("We are skipping the after-edgeR post-processing scripts, which depend on edgeR.\n");
	} else {
		my $dep  = { 'edgeR'=>$cfg->{jobs}->{$STEP_EDGER_DIFF_EXPR}->{"all"}->{jobName} };
		my $vars = { 'force'=>$FORCE_RERUN, 'verbose'=>$GLOBAL_VERBOSE, 'debug'=>${GLOBAL_DEBUG}
			     ,'pbs_ncpus' => 1, 'pbs_mem' => "12g", 'pbs_walltime'=> "28821",
			     ,'subreadRawCountsFile' => catfile($cfg->{'subreadCountsDir'}, $cfg->{'subread'}->{'countsFile'})
			     ,'edgeRNaiveFpkmFile'   => $cfg->{'edgeR'}->{'fpkmFilePath'} # A specifically formatted FPKM matrix file generated by edgeR
			     ,'monocleDir'           => $cfg->{'monocleDir'} # (Full path)
			     ,'summaryFigureDir'     => $cfg->{'summaryFigureDir'} # (Full path)
			     ,'nGroups'              => getNumExperimentalGroups()
			   };

		if ($SHOULD_ATTEMPT_TO_RUN_MONOCLE) {
			setJobInfo($cfg, $STEP_MONOCLE            , "all", $dep,$vars,$cfg->{monkeyPoo}, "51a.monocle.R");
			writeup($cfg, qq{We applied Monocle [CITE_MONO] to the FPKM values calculated by edgeR.\n}
				, {"CITE_MONO"=>qq{Monocle: Trapnell C, Cacchiarelli D, Grimsby J, Pokharel P, Li S, Morse M, Lennon NJ, Livak KJ, Mikkelsen TS, and Rinn JL. "The dynamics and regulators of cell fate decisions are revealed by pseudo-temporal ordering of single cells." Nature Biotechnology (2014).} } );
		}
		setJobInfo($cfg, $STEP_SUMMARY_AFTER_EDGER, "all", $dep,$vars,$cfg->{monkeyPoo}, "51b.summary.R");
		setJobInfo($cfg, $STEP_CLUSTER_AFTER_EDGER, "all", $dep,$vars,$cfg->{bin}->{RSCRIPT_EXE}," ",$cfg->{monkeyPoo}, "51c.cluster.R"); # output goes into the summaryFigureDir
		
		writeup($cfg, qq{We also generated PCA plots.\n}
			, {"CITE_PCA"=>qq{PCA plots: No citation yet.} } );
	}
	# ================================= [DONE] Using edgeR ============================================

	# Build quality control jobs (FastQC and RSEQC). Nothing depends on FastQC or RSEQC output, so this can be added last. Remember that normally all jobs must be added IN DEPENDENCY ORDER, or else Torque errors out.
	if ($cfg->{doQC}) {
		foreach my $h ($sampleHash,$inputHash) {
			next if (!defined($h)); # Skip any undefined 'inputHash' cases.
			while ( my ($k, $tHash) = each(%$h) ) { # $tHash (hash reference) is the value
				my @sKeys    = sort(keys(%$tHash)); # these are the plain filenames (without the path)
				my $jobPairIndex = 0;
				foreach my $filePrefix (@sKeys) {
					# ================================= RUN FASTQC ON UNMAPPED (but filtered!) FASTQ FILES =============
					$jobPairIndex++; # 1 or 2 depending on the pair in question
					my $dep  = { 'filterJob' => $cfg->{jobs}->{$STEP_FILTER}->{$k}->{jobName} };
					# somehow this dependency is not always found??
					my $vars = {'force'=>$FORCE_RERUN, 'verbose'=>$GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG,
						    'pbs_ncpus' => 1, 'pbs_mem' => "2g", 'pbs_walltime'=> "8989",
						    'outputDir'   =>  $cfg->{fastqcResultsDir} , # the PARENT output directory
						    'fastqc'      =>  getBinPath($cfg, "fastqc")    , # binary location
						    'inputFile'   =>  catfile($cfg->{filterResultsDir}, (noExtensions($filePrefix).".gz")) # <== Must be the file's FULL PATH. The input file will have been FILTERED at this point, and is always re-compressed with ".gz", no matter what the initial input compression was. This should ALWAYS be .gz, assuming that the "02.filtering" step generates a .gz file.
						   };
					setJobInfo($cfg, "${STEP_FASTQC_UNMAPPED}_${jobPairIndex}", "${k}",$dep,$vars,$cfg->{monkeyPoo},"49.qc.pl");
					# ================================= DONE RUNNING FASTQC ON UNMAPPED (but filtered!) FASTQ FILES ====
				}
			}
		}
		writeup($cfg, qq{Sequence quality control was assessed using the program FastQC [CITE_FASTQC].}, {"CITE_FASTQC"=>qq{FastQC: Andrews S. FastQC. http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ .}});
		my $shouldRunRSEQC = (defined($cfg->{bedAnnotForRSEQC}) && !isNA($cfg->{bedAnnotForRSEQC})); # We will run RSEQC *unless* the bed annotation file is 'NA' (or undefined), in which case we skip it
		foreach my $sample (keys(%{$remember{$REM_ALIGN_BAM}})) {
			my $prereqJob = $remember{$REM_ALIGN_BAM}{$sample}{"job"}; defined($prereqJob) or confess "[ERROR] failed to locate the required pre-requisite job variable!";
			my $bamFile   = $remember{$REM_ALIGN_BAM}{$sample}{"bam"}; defined($bamFile)   or confess "[ERROR] failed to locate the required bam variable!";
			my $dep       = { 'prereqAlignJob' => $prereqJob };
			# ================================= RUN FASTQC ON MAPPED BAM FILES INDIVIDUALLY ========================
			my $fqcVars = {'force'=>$FORCE_RERUN, 'verbose'=>$GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG,
				       'pbs_ncpus' => 1, 'pbs_mem' => "2g", 'pbs_walltime'=> "8989",
				       'outputDir'   =>  $cfg->{qcAfterMappingDir}, # the PARENT output directory
				       'fastqc'      =>  getBinPath($cfg, "fastqc")    , # binary location
				       'inputFile'   =>  $bamFile # <== Must be the file's FULL PATH
				      };
			setJobInfo($cfg, $STEP_FASTQC_ALIGNED, $sample,$dep,$fqcVars,$cfg->{monkeyPoo},"49.qc.pl");
			# ================================= DONE RUNNING FASTQC ON MAPPED BAM FILES INDIVIDUALLY ==============
			# ================================= RUN RSEQC ON MAPPED BAM FILES INDIVIDUALLY ========================
			if ($shouldRunRSEQC) {
				($GLOBAL_DRY_RUN or (-e "$cfg->{RSEQC}->{rseqc_parent_dir}/geneBody_coverage.py")) or confess "[ERROR] RSEQC wants to run, but we could not find any indication that its scripts actually exists in the alleged location, specifically: $cfg->{RSEQC}->{rseqc_parent_dir}";
				my $rseqVars = { 'force'=>$FORCE_RERUN, 'verbose'=>$GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG
						 ,'pbs_ncpus' => 1, 'pbs_mem' => "8g", 'pbs_walltime'=> "172800" # super long time is allocated here...
						 ,'outputDir'         => $cfg->{rseqcQualityDir}        # the output directory, which will be created
						 ,'inputBAM'          => $bamFile # <== Must be the input bam file's FULL PATH
						 ,'bedAnnot'          => $cfg->{bedAnnotForRSEQC} # location of bed-format annotation files for RSEQC
						 ,'python2_7_exe'     => getBinPath($cfg, "python2_7_exe")
						 ,'samtools'          => getBinPath($cfg, "samtools")
						 ,'scriptParentDir'   => $cfg->{RSEQC}->{rseqc_parent_dir} # location where all the RSEQC executables are
						 ,'nBamLines'         => 150000 # RSEQC is super slow, so only examine this subset of lines
						 ,'shouldUseAllReads' => 0 # Do not examine ALL reads in the bam file--that's super slow! If this is 1, then the 'nBamLines' is ignored.
						 ,'PYTHONPATH'        => getMostExtensivePythonPath()
					       };
				setJobInfo($cfg, $STEP_RSEQC_SINGLE, $sample, $dep,$rseqVars,$cfg->{monkeyPoo},"47.rse.qc.pl");
			}
			# ================================= DONE RUNNING RSEQC ON MAPPED BAM FILES INDIVIDUALLY ===============
		}
	    
		# ================================= RUN RSEQC ON ALL MAPPED BAM FILES **TOGETHER** ========================
		if ($shouldRunRSEQC) {
			# This seems to fail frequently for some reason... but it doesn't generate any useful error messages!
			my $dependOnAllBams; # hash pointer
			foreach my $sample (keys(%{$remember{$REM_ALIGN_BAM}})) {
				my $prereqJob = $remember{$REM_ALIGN_BAM}{$sample}{"job"}; defined($prereqJob) or confess "[ERROR] failed to locate the required pre-requisite job variable!";
				$dependOnAllBams->{$prereqJob} = $prereqJob; # save this as a prerequisite!
			}
			my $multiVars = {  'force'=>$FORCE_RERUN, 'verbose'=>$GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG
					   ,'pbs_ncpus' => 1, 'pbs_mem' => "16g", 'pbs_walltime'=> "345600"
					   ,'outputDir'       => $cfg->{rseqcQualityDir}         # the output directory, which will be created
					   ,'inputBamDir'     => $cfg->{mappingResultsDir}       # we specify a DIRECTORY of bam files
					   ,'bedAnnot'        => $cfg->{bedAnnotForRSEQC}        # location of bed-format annotation files for RSEQC
					   ,'python2_7_exe'   => getBinPath($cfg, "python2_7_exe") # binary location for python 2.7 or later
					   ,'scriptParentDir' => $cfg->{RSEQC}->{rseqc_parent_dir} # location where all the RSEQC executables are
					   ,'PYTHONPATH'        => getMostExtensivePythonPath()
					};
			setJobInfo($cfg, $STEP_RSEQC_TOGETHER, "all",$dependOnAllBams,$multiVars,$cfg->{monkeyPoo},"46.qc.multiple.bam.pl");
			# ================================= DONE RUNNING RSEQC ON ALL MAPPED BAM FILES TOGETHER ========================
			writeup($cfg, qq{Additional quality control figures were generated using RSeQC [CITE_RS].}, {"CITE_RS"=>qq{Wang L, Wang S, Li W. "RSeQC: quality control of RNA-seq experiments." Bioinformatics (2012) 28 (16): 2184-2185. doi: 10.1093/bioinformatics/bts356 . PubMed ID 22743226.}});
		}
	} # Done with QC (FastQC an RSEQC) quality control
	# put further analyses here...
}

#########################################################
### 2. generate writeup #################################
#########################################################
sub writeup($$\%) {
	# Generates two files (or appends to them, if they exist): $cfg->{writeup} and $cfg->{writeup_cite}.
	# Fill out a legible writeup. Note that the writeup should only be modified
	# from *THIS* script! It should not be updated by any of the 'qsub' jobs, because then
	# it could be edited in a non-deterministic order. (It would be OK if they wrote to a bunch of temp files
	# that were concatenated in the end, but that's not how it works right now.)
	my ($cfg, $msg, $citationHashPtr) = @_;
	my $writeup = $cfg->{writeup};
	my $cite    = $cfg->{writeup_cite};
	my $openedFile = 0;
	if ($GLOBAL_DRY_RUN) {
		*OF    = *STDOUT;
		*CITEF = *STDOUT;
		$msg   = "DRY RUN WRITEUP >>>> $msg";
	} else {
		open (OF,">> $writeup") or confess "[ERROR] Cannot APPEND to the output 'writeup' file, named '$writeup'. Check permissions and/or filenames.";
		open (CITEF,">> $cite") or confess "[ERROR] Cannot APPEND to the output 'citation' file, named '$cite'. Check permissions and/or filenames.";
		$openedFile = 1;
	}
	print OF "$msg" . "\n";
	while (my($k,$v) = each(%$citationHashPtr)) {
		print CITEF "Citation: $k = $v" . "\n";
	}
	if ($openedFile) { close(OF); close(CITEF); } # If we opened a file, now CLOSE it. Don't try to close STDOUT!
}

#########################################################
### 3. run jobs #########################################
#########################################################

sub appendLinesInPlace($$$) {
	# Appends lines to '$strRef' (reference, not directly a string!!) and counts how many newlines were appended
	# You pass in 3 things:
	# 1. a *REFERENCE* (pointer) to the line-count variable (first of 2 arguments). Should be numeric and already initialized to 0 (or whatever).
	# 2. ...and a *REFERENCE* (pointer) to the string variable to update (second of 2 arguments). Should also already be a blank string (or whatever).
	# 3. ...and finally, the actual new text to append (with newlines, which are counted and used to update the line-count variable!)
	my ($newlineCountRef, $strRef, $newText) = @_;
	$$strRef .= $newText;
	my $numNewlinesInNewText = () = $newText =~ /\n/g; # <== count number of matches of "\n" in the output string
	$$newlineCountRef += $numNewlinesInNewText; 	# no return value
}
sub appendLinesInPlaceNL($$$) { return appendLinesInPlace($_[0], $_[1], $_[2]."\n"); } # NL = newline. This ALSO appends a newline, but is otherewise the same as 'appendLinesInPlace'

sub runJobs {	      # Generate AND RUN the qsub "list of jobs to submit" file.
	my ($cfg)   = @_;
	my $outfile = "$cfg->{studyName}_jobs.sh"; # Name of the shell script that will we run below    
	open (OF,">$outfile") or confess "[ERROR] Cannot write to the following 'qsub job list' file: \"$outfile\"";
	my $lnum = 0; # For debugging problematic qsub commands, we can output an 'echo' command that runs after every qsub job. These commands start on line 4 of the script, which is why this is initialized to 4.
	my $OUT_HEADER_PRINT = '';
	appendLinesInPlaceNL(\$lnum, \$OUT_HEADER_PRINT, "#!/bin/bash");
	appendLinesInPlaceNL(\$lnum, \$OUT_HEADER_PRINT, "set -e           # <= abort if ANY command returns non-zero value");
	appendLinesInPlaceNL(\$lnum, \$OUT_HEADER_PRINT, "set -u           # <= abort on undefined variables");
	appendLinesInPlaceNL(\$lnum, \$OUT_HEADER_PRINT, "set -o pipefail  # <= show failed exit codes properly");
	appendLinesInPlaceNL(\$lnum, \$OUT_HEADER_PRINT, "export SHELL=/bin/bash"); # Make sure to use the BASH shell, not just plain '/bin/sh'!


	if ($HOLD_UNTIL_ALL_JOBS_SUBMITTED && !$RUN_DIRECTLY_WITHOUT_TORQUE) {
		my $WALLTIME_FOR_WAIT_JOBS = "300"; # these really only take like 1 second to run
		my $IMMEDIATELY_HOLD_JOB   = "-terse -h "; # Starts the job in the HELD state, until it is 'qrls'-ed.
		# Note that this job gets HELD immediately upon submission! It can't run until it gets 'qrls'-ed later.
		appendLinesInPlaceNL(\$lnum, \$OUT_HEADER_PRINT, qq{${SIGNAL_TO_START_VARNAME}=`echo "sleep 1" | } . get_qsub_cmd($cfg, {'pbs_walltime'=>$WALLTIME_FOR_WAIT_JOBS, 'pbs_stderr'=>"/dev/null", 'pbs_stdout'=>"/dev/null" }) . qq{ ${IMMEDIATELY_HOLD_JOB} -V -N $cfg->{studyName}_Start_signal`});
	}
	print OF $OUT_HEADER_PRINT;
	#printf OF ($GLOBAL_VERBOSE) ? "set -o verbose  # <== verbose \n" : "# (Not setting verbose mode) \n";
	foreach my $stepName (sort keys(%{$cfg->{jobs}})) { # <== it is CRITICALLY IMPORTANT that the jobs are added to the script in sorted order!
		foreach my $sampleName (sort keys %{$cfg->{jobs}->{$stepName}} ) {
			my $qcmd  = $cfg->{jobs}->{$stepName}->{$sampleName}->{qsub};
			my $jname = $cfg->{jobs}->{$stepName}->{$sampleName}->{jobName};
			my $OUTPRINT = '';
			($GLOBAL_VERBOSE) and appendLinesInPlaceNL(\$lnum, \$OUTPRINT, qq{echo ''});
			if ($GLOBAL_VERBOSE && !$RUN_DIRECTLY_WITHOUT_TORQUE) {
				appendLinesInPlaceNL(\$lnum, \$OUTPRINT, qq{echo '${jname} dependencies:' });
				my @dependenciesArr = split(",", $remember{$REM_DEPENDENCIES_STR_COLON_DELIM}{$jname});
				foreach my $d (@dependenciesArr) {
					$d =~ s/^[\$]//; # Remove the leading '$' from each variable so it doesn't get auto-evaluated when we $echo it
					appendLinesInPlaceNL(\$lnum, \$OUTPRINT, qq{echo "    * The dependency variable \"$d\" is => \$$d) : Confirming that this is a real job with 'qstat -f -j \$$d' "}); # note that the 'd-with-dollar-sign' gets EVALUATED since it has a dollar sign. So this will print something like "Dependency result was: 5928.machine" and not the actual dependency name.
					appendLinesInPlaceNL(\$lnum, \$OUTPRINT, qq{qstat -j \$$d > /dev/null });

					appendLinesInPlaceNL(\$lnum, \$OUTPRINT, qq{echo "         * 'qstat' result: \$? (should be 0)"});
				}
			}
			($GLOBAL_VERBOSE) and appendLinesInPlaceNL(\$lnum, \$OUTPRINT, qq{echo '[OK] About to add the following job:  ${jname}'...});
			my $qdebug = $qcmd;
			$qdebug =~ s/\'/\'\\\'\'/g; # ehcos the commands with **quoted single-quotes** (so they show up properly) for diagnostic purposes
			#appendLinesInPlaceNL(\$lnum, \$OUTPRINT, qq{echo 'Submitting a command that looks mostly like this (multi-line) command: $qdebug'} . qq{\n} . qq{echo ''}); # Each job has a 'qsub' command associated with it, which we finally print here.
			appendLinesInPlaceNL(\$lnum, \$OUTPRINT, qq{$qcmd}); # Each job has a 'qsub' command associated with it, which we finally print here.
			($GLOBAL_VERBOSE) and appendLinesInPlaceNL(\$lnum, \$OUTPRINT, qq{echo }.'$?'.qq{ "was the exit code for the submission of job $jname ..."}); # Causes the shell to print the above command's exit code (shell variable is $?) to STDOUT
			if ($GLOBAL_VERBOSE && !$RUN_DIRECTLY_WITHOUT_TORQUE) {
				my $submitText = ($RUN_DIRECTLY_WITHOUT_TORQUE) ? "directly ran (without TORQUE)" : "submitted to SGE";
				appendLinesInPlaceNL(\$lnum, \$OUTPRINT, qq{echo "[OK] we ${submitText} the job '${jname}' (Step: '$stepName' for sample '$sampleName') (after line $lnum of <$outfile>) -> result = \"\$${jname}\""\n\n}); # extra newlines to mark the end of this job
			}
			print OF $OUTPRINT;
		}
	}

	if ($HOLD_UNTIL_ALL_JOBS_SUBMITTED) {
		my $FOOTER = '';  # It is crucial to RELEASE the blocking-everything job (with 'qrls')!
		appendLinesInPlaceNL(\$lnum, \$FOOTER, getBinPath($cfg, "qrls") . " \$${SIGNAL_TO_START_VARNAME}"); # <== run everything, by running the 'start signal' variable that all other jobs depend on (i.e. stop holding the job)
		print OF $FOOTER;
	}
	close(OF);
	
	print STDERR "\n";
	if ($GLOBAL_DRY_RUN) {                 print STDERR "Dry run: we are not ACTUALLY submitting any jobs. Check the file \"$outfile\" to see what WOULD have been submitted.\n"; }
	elsif ($RUN_DIRECTLY_WITHOUT_TORQUE) { print STDERR "Because '--no-torque' was specified on the command line, we are NOT submitting these jobs to the queue management software, but instead running a single job at a time, VERY SLOWLY in sequential order.\n"; }
	else {                                 print STDERR "Now submitting jobs in the file \"$outfile\" to the Torque job-management queue...\n"; }
	print STDERR "***********************************************************************\n";
	print STDERR "* Note: 1. You can type 'qstats' or 'qstat' to see your job status.\n";
	print STDERR "*       2a. To cancel ONE of your jobs, you can type 'qdel JOB_ID_GOES_HERE'!\n";
	print STDERR "*       2b. To cancel EVERY SINGLE job of yours, use 'qselect -u \$USER | xargs qdel'. DANGEROUS!\n";
	print STDERR "*       3. To see historical / completed jobs, try: qstat -x | less -S\n";
	print STDERR "*       4. To see extremely detailed job info, try 'qstat -f'\n";
	print STDERR "***********************************************************************\n";
	mkdirOrDie($LOGDIR); # The directory that all the torque output goes to. Must come BEFORE the system call.
	verboseSystem("chmod 711 $outfile"); # Allow the job-submission script to be executed by this user. We do this even in a dry run.
	my $result = verboseSystem("./$outfile");
	(0 == $result) or confess "[ERROR] in attempting to run the ./$outfile script! Exit code was '$result'. "; # <== This actually RUNS the jobs
	(!$GLOBAL_DRY_RUN && !$RUN_DIRECTLY_WITHOUT_TORQUE) && verboseOkPrint("Here are the first few jobs of yours (command is: qstat -u $ENV{USER}):\n" . `sleep 2; qstat -u \"$ENV{USER}\" | head -n 15` . "\n (Full job names available by running 'qstat -f -u $ENV{USER}' instead of 'qstat')...\n");
}

1;	 # <== Perl requires any module to end with a true value. Hence, this 1.
