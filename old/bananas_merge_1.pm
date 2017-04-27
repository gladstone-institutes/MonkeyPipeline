#!/usr/bin/perl
# bananas.pm version 2 - parallelized for rigel
# Authors: Sean Thomas and Alex Williams, 2013
# This program is designed to automate many of the tasks that accompany next generation sequencing tasks. 

# Alex's changed package name so as to not inadvertently interfere with Sean's.
# Thing to note:
#        * when jobs are queued up, they are added IN ALPHABETICAL ORDER
#        * a job that begins alphabetically earlier therefore CANNOT depend on a later job ('a' cannot depend on 'z')
#        * thus, jobs are numbered in such a way that alphabetization is related to prerequisites

package bananas_agw;

use strict; use warnings; use Carp; # Carp has "confess"
use Data::Dumper;
use File::Basename;
use Scalar::Util qw(looks_like_number);
use File::Spec::Functions qw(catfile); # File::Spec has "catfile"
use IO::Compress::Gzip qw(gzip $GzipError) ;
use List::Util qw(max);
use Exporter;

### 0.i Export functions
our (@ISA, @EXPORT, @EXPORT_OK);
@ISA       = qw(Exporter);
@EXPORT_OK = qw(GLOBAL_VERBOSE GLOBAL_DRY_RUN GLOBAL_DEBUG);
@EXPORT    = qw(loadConfig checkConfig buildJobSubmissionScript generateWriteup runJobSubmissionScript);

### 0.ii Option choices
our ($GLOBAL_VERBOSE)      = 0; # 0 = quiet, 1 = verbose
our ($GLOBAL_DEBUG)        = 0; # 0 = normal, 1 = ultra-verbose diagnostic messages
our ($GLOBAL_DRY_RUN)      = 0; # 0 = normal, 1 = dry run (dry = "don't actually do anything")
our ($GLOBAL_COLOR)        = 0; # 0 = only black and white messages. 1 = COLOR text to console
our ($FORCE_RERUN)         = 0; # 0 = normal, 1 = re-run all jobs, even completed ones
our ($OUTPUT_DIR_OVERRIDE) = undef; # If defined, this means the user has specified their own resultDir on the command line, rather than in the config file.

# List recognized "key" names from the config file. Anything in the config file but NOT in here will raise a warning.
my @OK_KEY_NAMES = qw(studyName sampleDir resultDir forceRerun minMapQ libraryAdapterFile bowtie2Index genomicBins genomicBinsID chromSizesFile geneWindows exonFile transcriptFile symbolXref genomeFasta repeatMask monkeyPoo aligner gtfFile tracksURL tracksSCP skipFastQC skipQC skipBrowser browserBins browserWigs skipWindows skipTagsAndDensity bedAnnotForRSEQC);
my @OK_ALIGNERS  = qw(bowtie tophat);

my $BANANAS_VERSION = "1.00a";

my $globalNumJobsSubmitted = 0; # counter

# =================== "REMEMBER" stores various file relationships (e.g., which BAM files were generated) which are not stored anywhere else ===================
my $REM_ALIGN_BAM = "alignbam"; # where we remember which bam files we aligned
my $REM_EXPERIMENT_CATEGORY = "category"; # remember which experimental categories exist
my $REM_SCRIPTS_CHECKED = "scripts-checked-for-ok-syntax";
my %remember = (  $REM_ALIGN_BAM          =>{} # new empty hash
		, $REM_EXPERIMENT_CATEGORY=>{} # new empty hash
		, $REM_SCRIPTS_CHECKED    =>{} # new empty hash
    ); # A hash that stores important file relationships that we need to remember. Somewhat ad hoc at the moment.
# ==============================================================================



#########################################################
# Code for the sub-scripts to (for example) evaluate arguments

# envLooksTrue is also used by the sub-scripts in the bin directory--do not delete it!
sub envLooksTrue($) {
    # Returns whether $ENV{'the_input_thing'} is a TRUE perl value, and NOT "false" or "undef" (case-insensitive).
    # Example correct usage:   my $boolean = envLooksTrue("varname");      # Example wrong usage:     my $wrongAnswer = envLooksTrue($ENV{"varname"} <- WRONG ) <-- WRONG! Do not pass "$ENV" in here!
    return (defined($ENV{$_[0]}) && $ENV{$_[0]} && (uc($ENV{$_[0]}) ne "FALSE") && (uc($ENV{$_[0]}) ne "UNDEF"));
}

# mkdirOrDie is also used by the sub-scripts in the bin directory--do not delete it!
sub mkdirOrDie($) { my ($dir) = @_; ((-d $dir) or mkdir($dir)) or confess "Failure in our attempt to make a directory! Unable to find or create the directory '$dir'!"; }

# systemAndLog is also used by the sub-scripts in the bin directory--do not delete it!i
sub systemAndLog($;$$) {
    # Takes a command and (optionally) a boolean (should log to STDERR?) and a filename (if defined, append to this log file?)
    # Returns the exit code of the command, just like plain 'system'
    my ($cmd, $shouldLogToStderr, $logToThisFile) = @_;
    if (!defined($shouldLogToStderr)) { $shouldLogToStderr = 0; }
    my $date = `date`; chomp($date);
    my $msg = "[${date}] Running this command:\n     ${cmd}\n";
    my $exitCode = system($cmd);
    $msg .= "     (Returned exit code $exitCode)\n";
    if ($exitCode != 0) { $msg .= "[WARNING]: Exit code was nonzero!\n"; }
    if ($shouldLogToStderr) { print STDERR $msg; }
    if (defined($logToThisFile)) {
	open(LOG,">>$logToThisFile") or die "ERROR: Failed to open log file $logToThisFile. Probably this is a programming error. Fix it!";
	if ($logToThisFile) { print LOG $msg; }
	close(LOG);
    }
    return($exitCode);
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
    defined(%ENV) or confess "The global 'ENV' hash variable (which should have been passed in by the 'qsub' queue submission program) was *not* defined in script submission! This indicates a fatal error somewhere in the queue submission process. Possibly this job was run by some means OTHER than qsub?";
    for my $requiredVarName (@varList) {
	defined($ENV{$requiredVarName}) or confess "The required variable '$requiredVarName' was NOT passed into this script by bananas.pm! This is a bug in bananas--fix it! Also check your capitalization--variable names are CASE SENSITIVE.";
    }
    return 1; # <-- If we get here, then all the required variables were defined.
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

sub dieIfFileAccessFails($;$) {
    my ($filename, $optionalMessage) = @_;
    if (!$GLOBAL_DRY_RUN) {
	if (!defined($optionalMessage)) { $optionalMessage = "FILE HANDLING FAILED"; }
	open(THEFILE,$filename) or confess "$optionalMessage: Cannot open file with name '$filename'! This is a serious problem, and we must quit now.";
	close(THEFILE);
    }
}

# aid with checking configuration file key values
sub addKeyValuePairToConfig($$$$$) {
    my ($cfgPointer, $key, $value, $lineNum, $configFileName) = @_;
    (!exists(${$cfgPointer}->{$key})) or confess "\nDUPLICATE ASSIGNMENT OF '$key' DETECTED (filename: \"$configFileName\")\n>>>> The '$key' was defined a SECOND time on line $lineNum... ";
    ("$key" ~~ @OK_KEY_NAMES) or warn "\n>>>> KEY WARNING: the key '$key' (on line $lineNum in the config file $configFileName) was NOT FOUND in the list of recognized keys.";
    ${$cfgPointer}->{$key} = $value;
}

# Same as regular "print" in most regards, except it ONLY actually prints if the $GLOBAL_VERBOSE flag is set. Only accepts one string argument. No arrays!
sub verbosePrint($;$) { # If it is a dry run, then it prepends [DRY RUN] to the statement.
    my ($text, $color) = @_;
    ($GLOBAL_VERBOSE || $GLOBAL_DEBUG) && printColorErr($text, $color);  # Prints to STDERR always!
}

sub debugPrint($;$) {
    my ($text, $color) = @_;
    ($GLOBAL_DEBUG) && printColorErr($text, (defined($color)?$color:"yellow on_red"));  # Prints yellow on_red if no color was specified
}

sub isLiteralStringTrue($) {    # Returns "TRUE" if input is both defined and the string "TRUE", in any capitalization, or 0 otherwise.
    my ($inputString) = @_;
    return ((defined($inputString) && ("TRUE" eq uc($inputString))) ? "TRUE" : 0);
}

sub isBooleanStringOrMissing { 
    my ($value) = @_;     # Returns "TRUE" if and only if the input is TRUE/true, FALSE/false, or undefined. Used for checking validity of inputs.
    return ((!defined($value) || (uc($value) eq "FALSE") || (uc($value) eq "TRUE")) ? "TRUE" : 0);
}

sub sanitizeJobName($) {
    my ($name) = @_;
    $name =~ tr/\./\_/; # Replaces any periods by underscores to make a valid job name. If the job starts with a number, warn the user and quit
    if ($name =~ /^\d/) { confess "Job names CANNOT start with a number--they have to start with a letter. Your invalid job name was: $name. "; }
    return($name);
}

sub getHashPointerOfAllBamJobs() { # (uses a global variable: the %remember hash)
    # Gets a valid "dependencies" hash that requires that all the BAM jobs have completed
    # This is used for any job that requires the BAM files be aligned, but doesn't have any other requirements.
    my $dependOnAllBams; # hash pointer
    my %allBamHash = %{$remember{$REM_ALIGN_BAM}}; # <-- global variable!
    foreach my $sample (keys(%allBamHash)) {
	my $prereqJob = $allBamHash{$sample}{"job"}; defined($prereqJob) or confess "failed to locate the required pre-requisite job variable!";
	$dependOnAllBams->{$prereqJob} = $prereqJob; # save this as a prerequisite!
    }
    return $dependOnAllBams; # hash pointer
}

sub noExtensions($) {
    my ($filename) = @_;
    if (!defined($filename)) { return undef; }
    $filename =~ s/[.](bz2|bzip2|gzip|gz|zip)$//i;
    return($filename);
}

# Runs a system command, and prints it out if we are using "verbose" mode (i.e., monkey.pl --verbose)
# DRY RUN behavior: Does NOT actually run the system call if this is a dry run!
sub verboseSystem($) { 
    my ($cmd) = @_;
    verbosePrint(">>>>> System call >>>>> $cmd\n", "green");
    my $exitCode = ($GLOBAL_DRY_RUN) ? 0 : system($cmd); # Only run the actual system command if it is NOT a dry run!
    verbosePrint("      (Returned exit code $exitCode)\n", "green");
    return($exitCode);    # <-- Mandatory that we RETURN the system call result!
}

sub loadConfig {
    # loads configuration file information into the $cfg hash
    my ($file) = @_;
    my $cfg; # <-- hash
    my %dupeCheckHash = (); # used to check whether an ID is a dupe  
    my $lineNum = 0;
    open(CF,$file) or die "ERROR: INVALID CONFIG FILE SPECIFIED: Cannot open the configuration file \"$file\" !";
    while (my $line = <CF>) { 
	chomp($line);
	$lineNum++;
	$line =~ s/\s+[#].*//;         # <-- Remove non-full-line comments! Note that comments need a space before them! Example: turns "mydata something    # comment here" into "mydata something"
	$line =~ s/[ \t\s]+/\t/g;      # <-- Replace RUNS of spaces/tabs with a single tab. Note that ALL SPACES ARE TURNED INTO TABS here! That means you can't have errant spaces in (say) filenames and whatnot.
	$line = trimWhitespace($line); # <-- Remove any other leading / trailing whitespace.
	next if ($line =~ /^[#]/);     # <-- Skip any lines that START with a '#' (comment character
	next if ($line =~ /^(\s)*$/);  # <-- Skip any ENTIRELY whitespace lines, in other words, there's nothing left after our modifications above.
	my @t = split(/\t/,$line);
	if ($line =~ /^(sample|input)\s/i) { #(scalar @t == 5) or (scalar(@t) == 6)) { # should probably just check to see if the line starts with "^sample\t" or "^input\t"
	    (scalar(@t) == 5 or scalar(@t) == 6) or confess "[ERROR ON LINE $lineNum OF CONFIG FILE <$file>]: The line began with 'sample' or 'input', but it did NOT have 5 or 6 elements! All lines beginning with 'sample' or 'input' must have exactly 5 (for single end) or 6 (for paired end) values. The text of the offending line was:\n$line";
	    
	    my $sampleOrInput  = lc($t[0]); # either "sample" or "input" -- the literal text!
	    my $idName         = $t[1]; # like "myDrugTest.1"
	    my $dataType       = lc($t[2]); # rna, chip, or exo. Lower-cased!
	    my $inputField     = $t[3]; # if it's CHIP-SEQ, then the "input" is a corresponding file that matches this sample. Or it can be "NA"
	    my $firstPairFile  = $t[4];
	    my $secondPairFile = (scalar(@t) >= 6) ? $t[5] : undef;

	    ($idName !~ m/[.].*[.]/  )         or confess "[ERROR ON LINE $lineNum OF CONFIG FILE <$file>]: The sample ID CANNOT have more than one decimal point in it: <$idName>! ";
	    ($idName =~ m/[.][0-9]+$/)         or confess "[ERROR ON LINE $lineNum OF CONFIG FILE <$file>]: The name \"$idName\" DID NOT end in a period and then a replicate number. Examples of valid names: WT.1 WT.2 DRUG.1 DRUG.2 DRUG.3. ";
	    (!exists($dupeCheckHash{$idName})) or confess "[ERROR ON LINE $lineNum OF CONFIG FILE <$file>]: The ID name for each sample must be UNIQUE, but the ID \"$idName\" on line $lineNum had already been used on line $dupeCheckHash{$idName}. Fix this! ";
            $dupeCheckHash{$idName} = $lineNum; # <-- Save the line number that this ID occurred on. Used to report any erroneous duplicates.

	    my @id_replicate_split = split(/\./, $idName);
	    (scalar(@id_replicate_split) == 2) or confess "[ERROR ON LINE $lineNum OF CONFIG FILE <$file>]: The sample/input name must be of the format 'name.replicateNumber'. The offending name was \"$idName\" . An example of a good name is \"WILDTYPE.2\" or \"DRUG.1\". Please start counting replicates from ONE!";
	    $cfg->{details}->{$idName}->{name}      = $id_replicate_split[0]; # example: "WILDTYPE" or "CONTROL"
	    $cfg->{details}->{$idName}->{replicate} = $id_replicate_split[1]; # example: "1" (for first replicate of WILDTYPE)
	    ($cfg->{details}->{$idName}->{name} =~ m/^[-_A-Za-z0-9]+$/) or confess "[ERROR ON LINE $lineNum OF CONFIG FILE <$file>]: The sample name <$cfg->{details}->{$idName}->{name}> is NOT acceptable, because it must only consist of standard alphanumeric characters plus the '-' and '_'. No special characters are allowed, as they may result in problematic filenames downstream. Please change this name to something without special characters, spaces, or punctuation";

	    # get data type <chip|rna|exo>
	    $cfg->{details}->{$idName}->{type} = $dataType;
	    if    ($dataType eq "chip") {   # chip-seq data
		$cfg->{doPeaks}      = 1;
		$cfg->{doExpression} = 0; # <-- don't do differential expression
		$cfg->{doExo}        = 0; # this variable remains unused as well
	    } elsif ($dataType eq "rna") {  # rna-seq data
		$cfg->{doPeaks}      = 0; # <-- don't peak-call RNA-seq!
		$cfg->{doExpression} = 1; # <-- DO differential expression with cuffdiff / whatever
		$cfg->{doExo}        = 0; # <-- don't exo-call RNA-seq!
	    } elsif ($dataType eq "exo") {  # ChIP-exo data
		$cfg->{doPeaks}      = 1;
		$cfg->{doExpression} = 0;
		$cfg->{doExo}        = 1;
	    } else { die "INPUT SPECIFICATION ERROR in file <$file> on line $lineNum: the sample type (which must be either 'chip', 'exo', or 'rna') is unrecognized. The invalid type was: <$dataType>"; }
	    
	    if ((uc($inputField) ne "NA") and ($inputField ne '') and ($inputField =~ /\S/)) { # If this field isn't blank or NA or has at least one non-space character (\S)
		$cfg->{details}->{$idName}->{input} = $t[3];
	    }

	    $cfg->{$sampleOrInput}->{$idName}->{$firstPairFile} = 1; # indicate that we have obtained a sample with this specific name
	    if (defined($secondPairFile)) { $cfg->{$sampleOrInput}->{$idName}->{$secondPairFile} = 1; } # indicate that we have obtained a paired end file!!!
	} elsif ($line =~ /^([^=]+)=(.*)$/) { # look for something with an EQUALS SIGN between it, e.g. "genome = myGenome.fa"
	    my $key   = trimWhitespace($1); (length($key) > 0)   or confess "[ERROR ON LINE $lineNum OF CONFIG FILE <$file>]: The key was BLANK or otherwise invalid. The offending line was:\n   $line\nPlease double-check the config file";
	    my $value = trimWhitespace($2); (length($value) > 0) or confess "[ERROR ON LINE $lineNum OF CONFIG FILE <$file>]: The value for key '$key' was BLANK or otherwise invalid. Note that every assignment MUST have a value, even if it is a placeholder like 'NONE'. The offending line was:\n   $line\nPlease double-check the config file";
	    addKeyValuePairToConfig(\$cfg, $key, $value, $lineNum, $file);
	} elsif (scalar(@t) == 2) { # look for any line with TWO ELEMENTS on it (i.e. tab-delimited key and value)
	    my $key = $t[0]; my $value = $t[1];
	    addKeyValuePairToConfig(\$cfg, $key, $value, $lineNum, $file);
	} else {
	    confess "[ERROR ON LINE $lineNum OF CONFIG FILE <$file>]: A line in the config file could not be interpreted. Perhaps it is incorrectly formatted, or has a 'weird' character in it somewhere. The offending line was:\n   $line\nPlease double-check the config file";
	}
    }
    close(CF);
    return($cfg);
}

sub trimWhitespace { my $str=shift; $str =~ s/^\s+|\s+$//g; return($str) }; # trims leading and trailing whitespace

sub checkConfig {
    my ($cfg) = @_;

    (defined($cfg->{monkeyPoo})) or die "Configuration file does not contain a 'monkeyPoo' entry--please specify this directory! ";
    (-d $cfg->{monkeyPoo})       or die "The specified 'monkeyPoo' directory ($cfg->{monkeyPoo}) does not exist. Double check that directory! ";

    (defined($cfg->{studyName}))              or die "Configuration file must contain a 'studyName' entry. This is just a basic name for your study, like 'myRnaSeqProject'. No periods / hyphens / unusual characters are allowed! ";
    ($cfg->{studyName} !~ /[\/\\\s\t\-\,\.]/) or die "Study name cannot contain whitespaces/slashes/backslashes/hyphens/periods. No periods! You can use underscores. The offending name was: \"$cfg->{studyName}\" ";

    if ($cfg->{studyName} =~ /^[0-9]/) {
	$cfg->{studyName} = "X".$cfg->{studyName}; # Add an "X" if the studyName *starts* with a number. "qsub" can't use job names that start with numbers.
	warn(qq{WARNING: Since the 'studyName' cannot START with a number, we added an 'X' to the beginning---study name is now \"$cfg->{studyName}\".});
    }

    (defined($cfg->{sampleDir}))  or die "Configuration file must contain a 'sampleDir' entry!";
    ($cfg->{sampleDir} =~ /^\//)  or die "Sample directory must be a FULL PATHNAME (i.e., /path/to/place) and NOT a relative path! The offending name was: $cfg->{sampleDir}!"; # check to make sure it starts with a '/' (and is therefore probably a full path)
    ($cfg->{sampleDir} !~ /[\s]/) or die "Sample directory must NOT contain whitespaces. The offending name was: $cfg->{sampleDir}!";
    (-d $cfg->{sampleDir})        or die "Sample directory must ALREADY exist, but this directory appears not to: $cfg->{sampleDir}!";

    if (defined($OUTPUT_DIR_OVERRIDE)) { $cfg->{resultDir} = $OUTPUT_DIR_OVERRIDE; } # <-- the user can specifiy --out=/my/directory/path on the command line to override the path in the study design file.
    (defined($cfg->{resultDir}))  or die "Configuration file must contain a 'resultDir' entry!";
    ($cfg->{resultDir} =~ /^\//)  or die "Result directory must be a FULL PATHNAME (i.e., /path/to/place) and NOT a relative path! The offending name was: $cfg->{resultDir}!"; # check to make sure it starts with a '/' (and is therefore probably a full path)
    ($cfg->{resultDir} !~ /[\s]/) or die "Working directory cannot contain whitespaces. The offending name was: $cfg->{resultDir}";
    mkdirOrDie($cfg->{resultDir});

    #$cfg->{writeup} = catfile($cfg->{resultDir}, "00.writeup.txt");

    $cfg->{filterResultsDir}  = catfile($cfg->{resultDir},"02a.filtered_fasta");
    $cfg->{fastqcResultsDir}  = catfile($cfg->{resultDir},"02b.filtered_fasta_fastqc");
    $cfg->{mappingResultsDir} = catfile($cfg->{resultDir},"04a.mapping");
    $cfg->{qcAfterMappingDir} = catfile($cfg->{resultDir},"04b.mapped_qc_fastqc");
    $cfg->{rseqcQualityDir}   = catfile($cfg->{resultDir},"04c.mapped_qc_rseqc");
    $cfg->{tagResultsDir}     = catfile($cfg->{resultDir},"05.tags");
    $cfg->{densityResultsDir} = catfile($cfg->{resultDir},"06.tagDensity");
    $cfg->{peakResultsDir}    = catfile($cfg->{resultDir},"07.peaks");
    $cfg->{browserResultsDir} = catfile($cfg->{resultDir},"08.browser_tracks");
    $cfg->{windowResultsDir}  = catfile($cfg->{resultDir},"09.window_figures");
    $cfg->{motifDiscDir}      = catfile($cfg->{resultDir},"10.motif_discovery");
    $cfg->{cuffdiffDir}       = catfile($cfg->{resultDir},"22.cuffdiff");

    # check tools
    my $binDir = "/data/tools/bin"; # <-- note: hard-coded.
    $cfg->{bin}->{basedir}    = $binDir;
    $cfg->{bin}->{convertToGenomeBrowserExe} = catfile($binDir,"convert_SAM_or_BAM_for_Genome_Browser.pl");
    $cfg->{bin}->{bam2bed}    = catfile($binDir,"bam2bed");
    $cfg->{bin}->{bigWig}     = catfile($binDir,"bedGraphToBigWig");
    $cfg->{bin}->{bedmap}     = catfile($binDir,"bedmap");
    $cfg->{bin}->{bedops}     = catfile($binDir,"bedops");
    $cfg->{bin}->{bigBed}     = catfile($binDir,"bedToBigBed");
    $cfg->{bin}->{bedtools}   = catfile($binDir,"bedtools");
    $cfg->{bin}->{bowtie2}    = catfile($binDir,"bowtie2");
    $cfg->{bin}->{peakMotifs} = "/data/applications/rsat/rsat/perl-scripts/peak-motifs"; # must be hardcoded to work
    $cfg->{bin}->{fastqMcf}   = catfile($binDir,"fastq-mcf");
    $cfg->{bin}->{fastqc}     = catfile($binDir,"fastqc");
    $cfg->{bin}->{match}      = catfile($binDir,"match");
    $cfg->{bin}->{match2moa}  = catfile($binDir,"match2moa");
    $cfg->{bin}->{moaOverlaps}= catfile($binDir,"removeMoaOverlaps");
    $cfg->{bin}->{samtools}   = catfile($binDir,"samtools");
    $cfg->{bin}->{sortBed}    = catfile($binDir,"sort-bed");
    $cfg->{bin}->{step1Motif} = catfile($binDir,"s1_getMotifRegex");
    $cfg->{bin}->{step2Motif} = catfile($binDir,"s2_findLocations");
    $cfg->{bin}->{tophat}     = catfile($binDir,"tophat");
    $cfg->{bin}->{cuffdiff}   = catfile($binDir,"cuffdiff");
    $cfg->{bin}->{rseqc_parent_dir} = "/usr/local/bin"; # Location of the RSEQC python files
    $cfg->{bin}->{python2_7_exe}    = "/bin/python2.7"; # location of a python executable with version >= 2.7

    my $th = $cfg->{bin};
    while ( my ($exeName, $exePath) = each(%$th)) {
	(-e $exePath) or confess "[ERROR in bananas_agw.pm]: The '$exeName' executable was expected to be in '$exePath', but it cannot be located there! Check to make sure this file actually exists. If that is the wrong location for this executable, you may have to change the expected location in the code in 'bananas.pm'. ";
    }
    
    (defined($cfg->{libraryAdapterFile})) or die "MISSING CONFIG FILE OPTION: Configuration file must contain a 'libraryAdapterFile' entry!";
    dieIfFileAccessFails($cfg->{libraryAdapterFile}, "Cannot find the following libraryAdapterFile: $cfg->{libraryAdapterFile}!");

    (defined($cfg->{minMapQ})) or die "MISSING CONFIG FILE OPTION: Configuration file must contain a 'minMapQ' entry! This is the minimum map quality (MAPQ) that a read must posses in order to NOT be disqualified after alignment. It should be between 0 (meaning 'keep everything') and 100. A standard value is 30.";
    ($cfg->{minMapQ} >= 0 && $cfg->{minMapQ} <= 100) or die "The 'minMapQ' option in the config file must be between 0 and 100, inclusive. The invalid map quality score specified in the file was: $cfg->{minMapQ}";

    (!defined($cfg->{genomeFasta}) || (-e $cfg->{genomeFasta})) or die "INCORRECT FILE PATH IN CONFIG FILE: you specified 'genomeFasta' as the file '$cfg->{genomeFasta}', but that file appears not to exist";
    (!defined($cfg->{skipBrowser}) || (!defined($cfg->{browserBins}) && !defined($cfg->{browserWigs}))) or die "You have defined BOTH a 'skipBrowser' entry in the config file AND you have also defined either a 'browserBins' entry and/or a 'browserWigs' entry. These are mututally exclusive. Recommendation: REMOVE the 'skipBrowser' line from your config file!";

    (!defined($cfg->{skipQC}) || !defined($cfg->{skipFastQC})) or die "You have defined BOTH 'skipQC' and 'skipFastQC'! These are redundant--remove 'skipFastQC' from your config file!";
    if (defined($cfg->{skipFastQC})) { $cfg->{skipQC} = $cfg->{skipFastQC}; } # This is just another name for skipQC

    for my $boolOption ('browserBins', 'browserWigs', 'skipBrowser', 'skipWindows', 'skipTagsAndDensity', 'skipQC','forceRerun') {
	isBooleanStringOrMissing($cfg->{$boolOption}) or die "Configuration file problem: '$boolOption' was defined as \"$cfg->{$boolOption}\", but it must be the literal text TRUE, FALSE, or omitted entirely. Blank values are NOT ACCEPTABLE!";
	$cfg->{$boolOption} = isLiteralStringTrue($cfg->{$boolOption});
    }

    if ($cfg->{'skipBrowser'}) { $cfg->{'browserWigs'} = 0; $cfg->{'browserBins'} = 0; } # 'skipBrowser' should not be used anymore, but in case it is, this is a (deprecated) shortcut to disabling both methods of browser track generation

    if (!defined($cfg->{bedAnnotForRSEQC})) { $cfg->{bedAnnotForRSEQC} = "(bed annotation for RSEQC was NOT DEFINED in the config file!"; }
    ((-e $cfg->{bedAnnotForRSEQC}) or $cfg->{skipQC}) or confess "Configuration file problem: unless you skip the QC step, you MUST have a bed annotation file defined as 'bedAnnotForRSEQC'--this file is necessary for rseQC to function! You either didn't define a variable, or you passed in the invalid file location '$cfg->{bedAnnotForRSEQC}'! Check your config file and verify that a file is in fact at that location.";
    ($cfg->{bedAnnotForRSEQC} =~ /.bed$/i) or confess "The 'bedAnnotForRSEQC' file ('$cfg->{bedAnnotForRSEQC}') did not end with .bed, so it's probably not the right file type! Fix this (or add skipQC=TRUE to the config file if you don't care about QC at all)!";

    # ======================== CHECK THE ALIGNER ====================
    (exists($cfg->{aligner})) or die "MISSING CONFIG FILE OPTION: The configuration file must specify an ALIGNER (the 'aligner') option, which was not specified!";
    $cfg->{aligner} = lc($cfg->{aligner}); # Convert the aligner name to LOWER CASE!
    ($cfg->{aligner} ~~ @OK_ALIGNERS) or die "The configuration file specified the following UNRECOGNIZED aligner: \"$cfg->{aligner}\". We expect the aligner to be something like 'bowtie' or 'tophat'.\n";
    if ("tophat" eq $cfg->{aligner}) {
	(defined($cfg->{gtfFile})) or die "MISSING CONFIG FILE OPTION: If you specify 'tophat' as your aligner, you must ALSO specify a valid gtf file with the 'gtfFile' parameter! You can also specify gtfFile=NONE in order to specifically indicate that there is no GTF file.";
	((uc($cfg->{gtfFile}) eq "NONE") or (-e $cfg->{gtfFile})) or die "The specified GTF-format annotation file (the 'gtfFile' setting in the config file) '$cfg->{gtfFile}' was NOT found on the filesystem. Tophat requires this file either exist OR be set to the special value of 'NONE' to indicate that there is no file..";
    }

    # ======================== CHECK THE INPUT files (for ChIP-seq usually)... if any) ====================
    if (defined($cfg->{input})) {
        my $th = $cfg->{input};
	while ( my($k,$th2) = each(%$th) ) {
	    foreach my $k2 (keys %$th2) { dieIfFileAccessFails(catfile($cfg->{sampleDir}, $k2), "The input file cfg->{sampleDir}/${k2} was missing or otherwise UNREADABLE"); }
        }
    }

    (defined($cfg->{sample})) or confess "Configuration file does not contain any 'sample' lines. e.g.:\nsample <sampleLabel> <fullPathToSampleFile>!";
    my $thsamp = $cfg->{sample};
    while ( my($k,$th2) = each(%$thsamp) ) {
	foreach my $sampkey (keys %$th2) { dieIfFileAccessFails(catfile($cfg->{sampleDir}, $sampkey), "SAMPLE FILE UNREADABLE -- maybe your sample directory was set incorrectly? It was set to: \"$cfg->{sampleDir}\"... double check it"); }
	my $th3 = $cfg->{details}->{$k};
	if (defined($th3->{input})) {
	    my $tch = $cfg->{input};
	    defined($tch->{$th3->{input}}) or confess "The 'Input' line that corresponds to sample '$k' lists a matched input ($th3->{input}) that isn't actually specified in the configuration file.";
	}
    }

    if (defined($cfg->{bowtie2Index})) {
	$cfg->{genomeName} = basename($cfg->{bowtie2Index});
    } else {
        confess "Configuration file does not contain a bowtie2Index entry, but that is MANDATORY! Even for tophat alignments! "
    }

    if ($cfg->{skipTagsAndDensity}) {
	($cfg->{doPeaks}) and print STDERR "[WARNING]: You are skipping tags and density---but you have data that needs peaks! By specifying 'skipTagsAndDensity', you skip EVERYTHING that relies on this downstream (browserBins, Sean's motif finding, etc.)\n";
	($cfg->{browserBins}) and print STDERR "[WARNING]: You are skipping tags and density---but you have data that needs peaks! By specifying 'skipTagsAndDensity', you skip EVERYTHING that relies on this downstream (browserBins, Sean's motif finding, etc.)\n";
	(!$cfg->{skipWindows}) and print STDERR "[WARNING]: You are skipping tags and density---but you have data that needs peaks! By specifying 'skipTagsAndDensity', you skip EVERYTHING that relies on this downstream (browserBins, Sean's motif finding, etc.)\n";

	$cfg->{doPeaks} = 0;     # Can't do peak calling without tags/density...
	$cfg->{browserBins} = 0; # Can't make browserBins without tags/density...
	$cfg->{skipWindows} = 1; # If you skip tags and density, you ALSO have to skip window-generation...
    }
    
    if ($cfg->{browserBins}) {
	dieIfFileAccessFails($cfg->{chromSizesFile}, "Configuration file does not contain a valid chromSizesFile file!");
	(defined($cfg->{tracksURL})) or confess "MISSING CONFIG FILE OPTION: browser track URL basename ('tracksURL') needs to be specified in the config file. Or, if you don't want browser tracks, add 'browserBins = FALSE' to your config file.";
	(defined($cfg->{tracksSCP})) or confess "MISSING CONFIG FILE OPTION: browser track scp location basename ('tracksSCP') needs to be specified in the config file.";
	dieIfFileAccessFails($cfg->{genomicBins}, "Configuration file does not contain a valid 'genomicBins' file!"); 	# Used by both Sean's browser and the tag density stuff
    }

    unless ($cfg->{skipTagsAndDensity}) {
	dieIfFileAccessFails($cfg->{genomicBins}, "Configuration file does not contain a valid 'genomicBins' file!"); 	# Used by both Sean's browser and the tag density stuff
	dieIfFileAccessFails($cfg->{genomicBinsID}, "Configuration file does not contain a valid genomicBinsID file!");
    }
    ##(open(SF,$cfg->{exonFile})       && close(SF)) or die "Configuration file does not contain a valid exonFile file!";
    ##(open(SF,$cfg->{transcriptFile}) && close(SF)) or die "Configuration file does not contain a valid transcriptFile file!";
    ##(open(SF,$cfg->{symbolXref})     && close(SF)) or die "Configuration file does not contain a valid symbolXref file, or that file could not be read.!";

    unless ($cfg->{skipWindows}) {
	(defined($cfg->{geneWindows}))              or die "MISSING CONFIG FILE OPTION: Configuration file must contain a 'geneWindows' entry!";
	dieIfFileAccessFails($cfg->{geneWindows}, "Cannot find the following geneWindows file (which was specified in the configuration file: '$cfg->{geneWindows}'");
    }

    if ($cfg->{doPeaks}) {
	(defined($cfg->{genomeFasta}))                                     or die "Configuration file does not contain a valid 'genomeFasta' entry! Specify one!";
	$cfg->{genomeFasta} =~ /\.(fa|fasta)(\.gz|\.bz2|\.gzip|\.bzip2)?$/ or die "The genome fasta file DOES NOT end in '.fa' or '.fasta' (or a gzipped/bzipped version of this name. The offending name was: $cfg->{genomeFasta}). ";
	dieIfFileAccessFails($cfg->{genomeFasta}, "Configuration file does have a valid 'genomeFasta' file (this must be a single fasta file!): it was $cfg->{genomeFasta}!");
	(defined($cfg->{repeatMask}))              or die "Configuration file does not contain a valid 'repeatMask' entry! Specify one!";
	dieIfFileAccessFails($cfg->{repeatMask} , "Configuration file does have a valid 'repeatMask' bed file: it was $cfg->{repeatMask}!");
    }

}

#########################################################
### 1. Build Job Submission Function ####################
#########################################################

sub setJobInfo {
    # job name, dependencies, variables, script name
    my ($jobName,$dep,$vars,$theBinDir,$script) = @_;
    ($jobName !~ /\./) or confess "The job name CANNOT CONTAIN A PERIOD / decimal point, but job '$jobName' did! Quitting. ";

    # ============ SET THE DEPENDENCIES ======================
    my $dependencies = (scalar(keys(%$dep)) > 0) ? "-W depend=afterok" : "";
    foreach my $k (sort(keys(%$dep))) {  # <-- Now list the dependencies
	($jobName ne $dep->{$k}) or confess "[ERROR in setting up qsub jobs]: You set the job name ($jobName) to THE SAME NAME as one of its prerequisites! This is not valid--a job cannot depend on itself. Fix it! ";
	(defined($dep->{$k}) and (length($dep->{$k}) > 0)) or confess "[ERROR in setting up qsub jobs]: You had a prerequisite for dependency <${k}> that was either BLANK or UNDEFINED. Fix this!";
	$dependencies .= ':' . '$' . $dep->{$k}; # note the LITERAL DOLLAR SIGN, which makes this dependency into a shell variable name (yes, really, that is how they are stored)
    }
    # ============ SET QSUB VARIABLES with -v ==================
    my $variables = "-v bananas='${BANANAS_VERSION}'"; # "-v" option supplies the variable list. Note that there is ALWAYS at least one variable (the bananas version)
    foreach my $key (sort(keys(%$vars))) {
	my $value = $vars->{$key};
	next if (!defined($value)); # Undefined values MUST be skipped---this is **critical**, otherwise they get included as the (literal) text 'undef'--which is NOT the same as Perl undefined!
	($value !~ "[\'\`]") or confess "[ERROR in setting up qsub jobs: SINGLE QUOTES / BACKTICKS NOT ALLOWED]: There was a single quotation mark or backtick in the key-value pair $key='$value' (in other words, <$value> had an invalid character in it!)! This will break the qsub command--quotes and backticks are NOT allowed in inputs!";
	($value !~ "[,]") or confess "[ERROR in setting up qsub jobs: COMMA NOT ALLOWED]: Torque/qsub is VERY PARTICULAR and cannot accept a comma in an argument (even inside double quotes). You have to REMOVE the comma, use some kind of placeholder text, and then manually re-add it in the script that you called. I recommend using <COMMA> and then replacing that with a ',' in the being-called script.";
	$variables .= qq{,}.qq{$key='$value'}; # <-- Variables are comma-delimited and are surrounded by single quotation marks.
    }
    # ============ SET THE SCRIPT PATH TO EXECUTE (usually a perl or shell script) ==================
    my $scriptFullPath = catfile($theBinDir,$script);
    (-e $scriptFullPath) or confess("The 'setJobInfo' function called for the script <$scriptFullPath>, but there did not appear to be any file actually there! Double check this! ");
    if ($scriptFullPath =~ /.pl$/i and (!$remember{$REM_SCRIPTS_CHECKED}{$scriptFullPath})) { # Is this a perl script that we have NOT checked to make sure its syntax is OK? If so, check it.
	(0 == system(("perl", "-c", $scriptFullPath))) or confess("[ERROR in syntax checking while setting up qsub jobs]: Looks like perl was not able to parse the syntax of the script <$scriptFullPath> (in other words, perl -c $scriptFullPath FAILED to return exit code zero. That script probably has an error that needs to be fixed");
	$remember{$REM_SCRIPTS_CHECKED}{$scriptFullPath} = 1; # Remember that this script was checked---don't bother checking it again
    }
    # ============ DONE! ==================
    my $qsubCmd = "${jobName}=`qsub -N $jobName $dependencies $variables $scriptFullPath `"; # <-- note the literal backtick! There will be a matching one later. Also: the saving-to-shell-variables part is SPACE-SENSITIVE, so DO NOT add spaces here whatever you do.
    $globalNumJobsSubmitted++;
    debugPrint("[DEBUG] Added job ${globalNumJobsSubmitted} via qsub: $qsubCmd\n", "cyan");
    return($jobName, $qsubCmd);
}

sub buildJobSubmissionList {
    # This function organizes all jobs that will need to be run and creates a script
    # that will be run in order to process all of the jobs in the most efficient way possible.
    # To add a new job, specify the base job name, dependencies, variables, script directory, and called script.

    # Note: jobs need to be added IN ORDER--pre-requisites must come BEFORE their dependencies here!
    # The scheduler is NOT clever enough to understand arbitrary ordering of prerequisites.

    my ($cfg) = @_;
    my $sampleHash = $cfg->{sample};
    my $inputHash  = (exists($cfg->{input})) ?  $cfg->{input}  : undef; # Set to undef if there's no input

    # build filtering jobs
    foreach my $h ($sampleHash,$inputHash) {
	next if (!defined($h)); # Skip any undefined 'inputHash' cases.
	foreach my $k (keys %$h) {
	    my $tHash   = $h->{$k};
	    my @sKeys   = sort(keys(%$tHash)); # sorted (so the first pair for paired end should always come first!) list of input filenames. Length 1 (single end) or 2 (paired end) only.
	    (scalar(@sKeys) == 1) or (scalar(@sKeys) == 2) or confess("Programming error: somehow sKeys was not length 1 or 2.");
	    my $jobName = sanitizeJobName("$cfg->{studyName}_s02_${k}");
	    my $dep; # <-- Note: no dependencies! This is the first real step!
	    my $vars = { 'force'=>$FORCE_RERUN, 'verbose'=>$GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG
			     , 'fastqMcf'    => ($cfg->{bin}->{fastqMcf})
			     , 'fastqDir'    => ($cfg->{sampleDir})
			     , 'filterDir'   => ($cfg->{filterResultsDir})
			     , 'adapterFile' => ($cfg->{libraryAdapterFile})
			     , 'inputFile1'  => $sKeys[0]
			     , 'inputFile2'  => ((2==scalar(@sKeys)) ? $sKeys[1] : undef) };
	    (-e catfile($vars->{'fastqDir'}, $vars->{'inputFile1'})) or confess("The file <$vars->{fastqDir}/$vars->{'inputFile1'}> was expected to exist, but we did NOT find it! This is a programming error in bananas");
	    (!defined($vars->{'inputFile2'})) or (-e catfile($vars->{'fastqDir'}, $vars->{'inputFile2'})) or confess("The file <$vars->{fastqDir}/$vars->{'inputFile2'}> was expected to exist, but we did NOT find it! This is a programming error in bananas");
	    my ($jn,$qs) = setJobInfo($jobName,$dep,$vars,$cfg->{monkeyPoo},"02.filter.pl");
	    $cfg->{jobs}->{'02_filter'}->{$k} = { "jobName"=>$jn, "qsub"=>$qs, "filtered_sample_file_path" => catfile($cfg->{filterResultsDir}, (noExtensions($sKeys[0]) . ".gz")) };
	    #verbosePrint("Filter-er: anticipating the creation of the file " . $cfg->{jobs}->{'02_filter'}->{$k}->{filtered_sample_file_path} . " ...\n", "yellow on_red");
	}
    }

    # build mapping jobs (map to the genome with tophat or bowtie)

    foreach my $h ($sampleHash,$inputHash) {
	next if (!defined($h)); # Skip any undefined 'inputHash' cases.
	foreach my $k (keys %$h) { # 'k' is the sample name--for example "DRUG.3" or "Wildtype.1"
	    my $tHash   = $h->{$k};
	    my @sKeys   = sort(keys(%$tHash));
	    my $jobName = sanitizeJobName($cfg->{studyName} . "_s04a_" . $k);
	    my $dep     = { 'filterJob' => $cfg->{jobs}->{'02_filter'}->{$k}->{jobName} }; # Depends on filtering being done already.

	    my $sampleType     = $cfg->{details}->{$k}->{type};
	    my $sampleCategory = $cfg->{details}->{$k}->{name};      # exprimental category, like "WILDTYPE" or "CONTROL"
	    my $sampleRep      = $cfg->{details}->{$k}->{replicate}; # Replicate number, like the "2" from "WILDTYPE.2" (which is "$k" here)

	    my $genomeShortName = basename($cfg->{bowtie2Index});
	    my $alignerExePath = (($cfg->{aligner} eq "tophat") && ($sampleType eq "rna")) ? "$cfg->{bin}->{tophat}" : "$cfg->{bin}->{bowtie2}";
	    my $vars = { 'filterDir'  => $cfg->{filterResultsDir}  ,
			 'mappingDir' => $cfg->{mappingResultsDir} ,
			 'aligner'    => $alignerExePath           ,
			 'samtools'   => $cfg->{bin}->{samtools}   ,
			 'gtfFile'    => $cfg->{gtfFile}           ,  # Optional, only used by Tophat (probably). Can be "undef"
			 'minMapQ'    => $cfg->{minMapQ}           ,
			 'bowtie2Index'=>$cfg->{bowtie2Index}      ,
			 'force'      => $FORCE_RERUN    ,
			 'verbose'    => $GLOBAL_VERBOSE ,
			 'debug'=>$GLOBAL_DEBUG          ,
			 'sampleName' => $k              ,
			 'inputFile1' => $sKeys[0]       ,
			 'inputFile2' => ((2==scalar(@sKeys)) ? $sKeys[1] : undef) ,   # <-- for paired end, get the second file this way
			 'sortMethod' => "BY_COORD"                                ,        # check 04.mapping.pl for valid values; should be "BY_COORD" or "BY_NAME" only!
			 'finalBam'   => catfile($cfg->{mappingResultsDir}, ("${k}_${genomeShortName}_q" . ${cfg}->{minMapQ} . ".bam"))
	    };

	    $remember{$REM_ALIGN_BAM}{$k} = {"job"=>$jobName, "bam"=>$vars->{finalBam}, "category"=>$sampleCategory, "replicate"=>$sampleRep }; # remember that we generate this BAM file, and remember some of the metadata that we will need later!
	    $remember{$REM_EXPERIMENT_CATEGORY}{$sampleCategory}{$sampleRep} = { "bam"=>$vars->{finalBam} }; # Remember this aligned bam file path, and that it was in this experimental category! We'll use this information for Cuffdiff.

	    my ($jn,$qs) = setJobInfo($jobName,$dep,$vars,$cfg->{monkeyPoo},"04.mapping.pl");
	    $cfg->{jobs}->{'04a_mapping'}->{$k} = { "jobName"=>$jn, "qsub"=>$qs };
	}
    }

    # ========= [START] SANITY CHECK THE INPUT CATEGORY (SAMPLE GROUP) NAMES: MAKE SURE THERE AREN'T ANY ACCIDENTAL CHANGES IN CAPITALIZATION =========
    my %seenCategories = (); # hash: key = lower-case version of the name, value = the original-case version of the name
    foreach my $ccc (keys(%{$remember{$REM_EXPERIMENT_CATEGORY}})) {
	if (exists($seenCategories{lc($ccc)}) && $seenCategories{lc($ccc)} ne $ccc) { confess "[ERROR in config file]: You had two experimental categories with the same name but DIFFERENT capitalization! This is NOT ALLOWED, because it's probably an error. Probably the capitalization of <$ccc> and <" . $seenCategories{lc($ccc)} . "> should be the same"; }
	$seenCategories{lc($ccc)} = $ccc; # Save the category name, but note that the KEY is the LOWER-CASED version of the name. This lets us detect wrong-capitalization versions of the otherwise-same experimental groups.
    }
    # ========= [DONE] SANITY CHECK THE INPUT CATEGORY (SAMPLE GROUP) NAMES =========

    # tags
    if (not $cfg->{skipTagsAndDensity}) {
	foreach my $h ($sampleHash,$inputHash) {
	    next if (!defined($h)); # Skip any undefined 'inputHash' cases.
	    foreach my $k (keys %$h) {
		my $inFile = catfile($cfg->{mappingResultsDir},"${k}_$cfg->{genomeName}_q$cfg->{minMapQ}.bam");
		my $dep  = { 'mappingJob' => $cfg->{jobs}->{'04a_mapping'}->{$k}->{jobName} };
		my $vars = { 'force' => $FORCE_RERUN, 'verbose' => $GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG };
		$vars->{genome} = $cfg->{genomeName};
		$vars->{seqType} = $cfg->{details}->{$k}->{type};
		$vars->{bam2bed} = $cfg->{bin}->{bam2bed};
		$vars->{sortBed} = $cfg->{bin}->{sortBed};
		$vars->{tagsDir} = $cfg->{tagResultsDir};
		$vars->{inputFile} = $inFile;
		$vars->{sampleName} = $k;

		my $jobName = sanitizeJobName($cfg->{studyName} . "_s05a_tags" . $k);
		my ($jn,$qs) = setJobInfo($jobName,$dep,$vars,$cfg->{monkeyPoo},"05.tags.pl");
		$cfg->{jobs}->{'05a_tags'}->{$k} = { "jobName"=>$jn, "qsub"=>$qs };
	    }
	}
    }

    # Tag mapping stats
    if (not $cfg->{skipTagsAndDensity}) {
	my $dep;
	my $th = $cfg->{jobs}->{'05a_tags'};
	foreach my $k (sort keys %$th) {
	    $dep->{$k} = $cfg->{jobs}->{'05a_tags'}->{$k}->{jobName};
	}
	my $vars = { 'force' => $FORCE_RERUN, 'verbose' => $GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG };
	$vars->{genome} = $cfg->{genomeName};
	$vars->{minMapQ} = $cfg->{minMapQ};
	$vars->{tagsDir} = $cfg->{tagResultsDir};
	$vars->{mappingDir} = $cfg->{mappingResultsDir};

	my $jobName = sanitizeJobName($cfg->{studyName} . "_s05b_summary"); # job names must be strictly in dependency alphabetical order!
	my ($jn,$qs) = setJobInfo($jobName,$dep,$vars,$cfg->{monkeyPoo},"05.xSummary.pl");
	$cfg->{jobs}->{'05b_tags_summary'}->{'xSummary'} = { "jobName"=>$jn, "qsub"=>$qs };
    }
    
    ## ChIP/exo specific processing
    # Tag density
    if (not $cfg->{skipTagsAndDensity}) {
	foreach my $k (keys %$sampleHash) {
	    my $jobName              = sanitizeJobName($cfg->{studyName} . "_s06_" . $k);
	    my $hasMatchingInput = (defined($cfg->{details}->{$k}->{input}));
	    my $dep = { 'sampleTagsJob' => $cfg->{jobs}->{'05a_tags'}->{$k}->{jobName} , 
			'inputTagsJob' => ($hasMatchingInput) ? $cfg->{jobs}->{'05_tags'}->{$vars->{inputName}}->{jobName}  : undef
	    }; # Depends on the SAMPLE tags being done AND maybe the INPUT tags being done
	    my $vars = { 'force' => $FORCE_RERUN, 'verbose' => $GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG };
	    $vars->{genome} = $cfg->{genomeName};
	    $vars->{seqType} = $cfg->{details}->{$k}->{type};
	    $vars->{bedmap} = $cfg->{bin}->{bedmap};
	    $vars->{binsFile} = $cfg->{genomicBins};
	    $vars->{binsFileID} = $cfg->{genomicBinsID};
	    $vars->{tagsDir} = $cfg->{tagResultsDir};
	    $vars->{densityDir} = $cfg->{densityResultsDir};
	    $vars->{sampleName} = $k;
	    $vars->{inputName}  = ($hasMatchingInput) ? $cfg->{details}->{$k}->{input}  :  undef;
	    my ($jn,$qs) = setJobInfo($jobName,$dep,$vars,$cfg->{monkeyPoo},"06.density.pl");
	    $cfg->{jobs}->{'06_density'}->{$k}->{jobName} = $jn;
	    $cfg->{jobs}->{'06_density'}->{$k}->{qsub} = $qs;
	}
    }

    if ($cfg->{doPeaks}) {
	# Generate the peaks for each sample. Depends on 05 and 06
	foreach my $k (keys %$sampleHash) {
	    # Note that we DO NOT generate peaks for type 'rna'!
	    if ($cfg->{details}->{$k}->{type} eq "chip" || $cfg->{details}->{$k}->{type} eq "exo") {
		my $jobName = sanitizeJobName($cfg->{studyName} . "_s07a_" . $k);
		my $dep = { "densityJob" => $cfg->{jobs}->{'06_density'}->{$k}->{jobName} }; # Each of these only depends on the ONE density job.
		my $vars = { 'force' => $FORCE_RERUN, 'verbose' => $GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG };
		$vars->{genome}     = $cfg->{genomeName};
		$vars->{seqType}    = $cfg->{details}->{$k}->{type};
		$vars->{bedmap}     = $cfg->{bin}->{bedmap};
		$vars->{bedops}     = $cfg->{bin}->{bedops};
		$vars->{sortBed}    = $cfg->{bin}->{sortBed};
		$vars->{tagsDir}    = $cfg->{tagResultsDir};
		$vars->{peaksDir}   = $cfg->{peakResultsDir};
		$vars->{densityDir} = $cfg->{densityResultsDir};
		$vars->{sampleName} = $k;
		my ($jn,$qs) = setJobInfo($jobName,$dep,$vars,$cfg->{monkeyPoo},"07.peaks.pl");
		$cfg->{jobs}->{'07a_peaks'}->{$k} = { "jobName"=>$jn, "qsub"=>$qs };
	    }
	}

	# Now generate the peak SUMMARY stats...
	# Note from Alex: I notice there is no explicit check for 'rna' here --- this code runs regardless, but I"m not sure it actually does anything in 'rna'-type samples.
	my $jobName = $cfg->{studyName} . "_s07b_xSummary";
	my $th = $cfg->{jobs}->{'07a_peaks'};
	my $dep;
	foreach my $k (sort keys %$th) {
	    $dep->{$k} = $cfg->{jobs}->{'07a_peaks'}->{$k}->{jobName}; # depends on all the peaks being done
	}
	my $vars = { 'force' => $FORCE_RERUN, 'verbose' => $GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG   ,
		     'genome'  => $cfg->{genomeName}       ,
		     'bedmap'  => $cfg->{bin}->{bedmap}    ,
		     'tagsDir' => $cfg->{tagResultsDir}    ,
		     'peaksDir'=> $cfg->{peakResultsDir}  };
	my ($jn,$qs) = setJobInfo($jobName,$dep,$vars,$cfg->{monkeyPoo},"07.xSummary.pl");
	$cfg->{jobs}->{"07b_xSummary"}->{'xSummary'} = { "jobName"=>$jn, "qsub"=>$qs };
    }

    ## Motif discovery module (only occurs if there are chip/exo data, and only on chip/exo data)
    if ($cfg->{doPeaks}) {
	my $numSamplesForMotifFinding = 0;
        foreach my $k (keys %$sampleHash) {
            if ($cfg->{details}->{$k}->{type} eq "chip" || $cfg->{details}->{$k}->{type} eq "exo") { 	# Alex's note: note that we do NOT do motif discovery for type 'rna'! Just 'chip' and 'exo'.
		$numSamplesForMotifFinding++;
		my $jobName = sanitizeJobName($cfg->{studyName} . "_s10_" . $k);
		my $dep  = { "peaksJob"=>$cfg->{jobs}->{'07a_peaks'}->{$k}->{jobName} }; # Depends on "07a_peaks"
		my $vars = { 'force' => $FORCE_RERUN, 'verbose' => $GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG };
                $vars->{genome} = $cfg->{genomeName};
                $vars->{bedops} = $cfg->{bin}->{bedops};
                $vars->{bedtools} = $cfg->{bin}->{bedtools};
                $vars->{peakMotifs} = $cfg->{bin}->{peakMotifs};
                $vars->{fasta} = $cfg->{genomeFasta};
                $vars->{repMask} = $cfg->{repeatMask};
                if ($cfg->{details}->{$k}->{type} eq "exo") {
                    $vars->{infile} = catfile($cfg->{peakResultsDir},"${k}_$cfg->{genomeName}_footprints.bed");
                } else {
		    $vars->{infile} = catfile($cfg->{peakResultsDir},"${k}_$cfg->{genomeName}_peaks.bed");
		}
                $vars->{motifsDir} = $cfg->{motifDiscDir};
                $vars->{sampleName} = $k;

                my ($jn,$qs) = setJobInfo($jobName,$dep,$vars,$cfg->{monkeyPoo},"10.motifDiscovery.pl");
		$cfg->{jobs}->{"10a_motifs"}->{$k} = { "jobName"=>$jn, "qsub"=>$qs };
            }
        }

	if ($numSamplesForMotifFinding > 0) { # summarize motif discovery, put all discovered motifs into a single transfac format file
	    my $jobName = $cfg->{studyName} . "_s10b_motif_summary";
	    my $th = $cfg->{jobs}->{'10a_motifs'};
	    my $dep;
	    foreach my $k (sort keys %$th) {
		$dep->{$k} = $cfg->{jobs}->{'10a_motifs'}->{$k}->{jobName};
	    }
	    my $vars;
	    $vars->{force} = $FORCE_RERUN;
	    $vars->{studyName} = $cfg->{studyName};
	    $vars->{motifDir} = $cfg->{motifDiscDir};
	    $vars->{peakDir} = $cfg->{peakResultsDir};
	    $vars->{bedmap} = $cfg->{bin}->{bedmap};
	    $vars->{bedops} = $cfg->{bin}->{bedops};
	    $vars->{bedtools} = $cfg->{bin}->{bedtools};
	    $vars->{genomeFasta} = $cfg->{genomeFasta};
	    $vars->{step1} = $cfg->{bin}->{step1Motif};
	    $vars->{step2} = $cfg->{bin}->{step2Motif};
	    $vars->{sortBed} = $cfg->{bin}->{sortBed};
	    my ($jn,$qs) = setJobInfo($jobName,$dep,$vars,$cfg->{monkeyPoo},"10.xSummary.pl");
	    $cfg->{jobs}->{"10b_motif_summary"}->{'xSummary'} = { "jobName"=>$jn, "qsub"=>$qs }; # <-- 10.2 comes AFTER 10.1
	}
    }
    ## end motif discovery module


    
    if ($cfg->{browserWigs}) {     # Alex's UCSC files -- generated if "browserWigs" is specified
	my $gExe = catfile($cfg->{bin}->{basedir}, "genomeCoverageBed");
	my $wExe = catfile($cfg->{bin}->{basedir}, "wigToBigWig");
	(-x $gExe) or confess "[ERROR]: FAILED TO FIND A RUNNABLE EXECUTABLE. Cannot find the required 'genomeCoverageBed' executable in path '$gExe'. Check to make sure this file exists **and** is executable by the current user!";
	(-x $wExe) or confess "[ERROR]: FAILED TO FIND A RUNNABLE EXECUTABLE. Cannot find the required 'wigToBigWig' executable in path '$wExe'. Check to make sure this file exists **and** is exectuable by the current user!";
	my $bamMultiFileString = "";
	my $dep = getHashPointerOfAllBamJobs(); # hash pointer -- require ALL bam alignment jobs are DONE!
	foreach my $sampleName (keys(%{$remember{$REM_ALIGN_BAM}})) {
	    my $singleBamFile   = $remember{$REM_ALIGN_BAM}{$sampleName}{"bam"}; defined($singleBamFile) or confess "failed to locate the required bam variable!";
	    $bamMultiFileString .= ${singleBamFile} . " "; # space-delimited!
	}
	my $jobName = sanitizeJobName($cfg->{studyName} . "_s48_agwbrowser_" . "everything");
	my $vars = {   'force'=>$FORCE_RERUN, 'verbose'=>$GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG,
		       'outputDir'                  =>  $cfg->{browserResultsDir}, # the PARENT output directory
		       'convertToGenomeBrowserExe'  =>  $cfg->{bin}->{convertToGenomeBrowserExe}    , # binary location
		       'inputMultiFileString'       =>  $bamMultiFileString  ,
		       'binDir'                     =>  $cfg->{bin}->{basedir} };
	#printf Dumper \$bamMultiFileString;
	my ($jn,$qs) = setJobInfo($jobName,$dep,$vars,$cfg->{monkeyPoo},"48.browser_agw.pl");
	$cfg->{jobs}->{'48_browser_wiggle'}->{"everything"} = { "jobName"=>$jn, "qsub"=>$qs }; # '49' so it gets alphabetized after all alignment is done, guaranteed.
    }

    if ($cfg->{browserBins}) {      # Sean's binned UCSC files -- generated if "browserBins" is specified
	foreach my $k (keys %$sampleHash) {
            my $jobName = sanitizeJobName($cfg->{studyName} . "_s08a_" . $k);
	    my $dep;
	    $dep->{densityJob} = $cfg->{jobs}->{'06_density'}->{$k}->{jobName};
	    if ($cfg->{details}->{$k}->{type} eq "chip" || $cfg->{details}->{$k}->{type} eq "exo") {
		$dep->{peaksJob} = $cfg->{jobs}->{'07a_peaks'}->{$k}->{jobName};
	    }
	    my $vars = { 'force' => $FORCE_RERUN, 'verbose' => $GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG };
	    $vars->{genome} = $cfg->{genomeName};
	    $vars->{seqType} = $cfg->{details}->{$k}->{type};
	    $vars->{bigWig} = $cfg->{bin}->{bigWig};
	    $vars->{bigBed} = $cfg->{bin}->{bigBed};
	    $vars->{peaksDir} = $cfg->{peakResultsDir};
	    $vars->{densityDir} = $cfg->{densityResultsDir};
	    $vars->{browserDir} = $cfg->{browserResultsDir};
	    $vars->{studyName} = $cfg->{studyName};
	    $vars->{sampleName} = $k;
	    $vars->{chromSizes} = $cfg->{chromSizesFile};
	    $vars->{genomicBins} = $cfg->{genomicBins};
	    $vars->{url} = $cfg->{tracksURL};
	    my ($jn,$qs) = setJobInfo($jobName,$dep,$vars,$cfg->{monkeyPoo},"08.browser.pl");
	    $cfg->{jobs}->{'08a_browser'}->{$k} = { "jobName"=>$jn, "qsub"=>$qs };
	}
	# After the browser tracks have all been added to the queue, generate the SUMMARY...
        my $jobName = sanitizeJobName($cfg->{studyName} . "_s08b_browser_summary");
        my $dep; # hash reference
        my $th = $cfg->{jobs}->{'08a_browser'};
        foreach my $k (sort keys %$th) {
	    $dep->{$k} = $cfg->{jobs}->{'08a_browser'}->{$k}->{jobName}; # Add dependencies
        }
	my $vars = { 'force'            =>  $FORCE_RERUN
			 , 'verbose'    =>  $GLOBAL_VERBOSE 
			 , 'scp'        => catfile($cfg->{tracksSCP},$cfg->{studyName})
			 , 'browserDir' => $cfg->{browserResultsDir} };
	my ($jn,$qs) = setJobInfo($jobName,$dep,$vars,$cfg->{monkeyPoo},"08.xSummary.pl");
	$cfg->{jobs}->{'08b_browser_summary'}->{$k} = { "jobName"=>$jn, "qsub"=>$qs };
    }

    # Generate windows
    unless ($cfg->{skipWindows}) {
	foreach my $k (keys %$sampleHash) { # Note: does NOT go through the 'inputhash'! Inputs don't need windows I guess.
            my $jobName = sanitizeJobName("$cfg->{studyName}_s09_${k}");
	    my $analyze = catfile($cfg->{monkeyPoo},"09.xAnalyze.RScript");
            my $dep     = { 'densityJob' => $cfg->{jobs}->{'06_density'}->{$k}->{jobName} };
	    my $vars = { 'force' => $FORCE_RERUN, 'verbose' => $GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG };
	    $vars->{genome} = $cfg->{genomeName};
	    $vars->{bedmap} = $cfg->{bin}->{bedmap};
	    $vars->{sortBed} = $cfg->{bin}->{sortBed};
	    $vars->{densityDir} = $cfg->{densityResultsDir};
	    $vars->{windowDir} = $cfg->{windowResultsDir};
	    $vars->{sampleName} = $k;
	    $vars->{geneWindows} = $cfg->{geneWindows};
	    $vars->{analyze} = $analyze;
	    my ($jn,$qs) = setJobInfo($jobName,$dep,$vars,$cfg->{monkeyPoo},"09.windows.pl");
	    $cfg->{jobs}->{'09_window'}->{$k} = { "jobName"=>$jn, "qsub"=>$qs };
	}
    }

    # Build quality control jobs (FastQC and RSEQC). Nothing depends on FastQC output, by the way, so this can be added last.
    unless ($cfg->{skipQC}) {
	foreach my $h ($sampleHash,$inputHash) {
	    next if (!defined($h)); # Skip any undefined 'inputHash' cases.
	    while ( my ($k, $tHash) = each(%$h) ) { # $tHash (hash reference) is the value
		my @sKeys    = sort(keys(%$tHash)); # these are the plain filenames (without the path)
		#my $isPaired = (scalar(@sKeys)==2); # paired end files, or not
		my $jobPairIndex = 0;
		foreach my $filePrefix (@sKeys) {
		    # ================================= RUN FASTQC ON UNMAPPED (but filtered!) FASTQ FILES =============
		    $jobPairIndex++; # 1 or 2 depending on the pair in question
		    my $jobName = sanitizeJobName($cfg->{studyName} . "_s49.2_unaligned_pair${jobPairIndex}_${k}");
		    my $dep  = { 'filterJob' => $cfg->{jobs}->{'02_filter'}->{$k}->{jobName} };
		    my $vars = {   'force'       =>  $FORCE_RERUN             ,
				   'verbose'     =>  $GLOBAL_VERBOSE          ,
				   'debug'       =>  $GLOBAL_DEBUG            ,
				   'outputDir'   =>  $cfg->{fastqcResultsDir} , # the PARENT output directory
				   'fastqc'      =>  $cfg->{bin}->{fastqc}    , # binary location
				   'inputFile'   =>  catfile($cfg->{filterResultsDir}, (noExtensions($filePrefix).".gz")) # <-- Must be the file's FULL PATH. The input file will have been FILTERED at this point, and is always re-compressed with ".gz", no matter what the initial input compression was. This should ALWAYS be .gz, assuming that the "02.filtering" step generates a .gz file.
		    };
		    my ($jn,$qs) = setJobInfo($jobName,$dep,$vars,$cfg->{monkeyPoo},"49.qc.pl");
		    $cfg->{jobs}->{'49_fastqc'}->{$k} = { "jobName"=>$jn, "qsub"=>$qs }; # <-- '49' so it gets alphabetized after all ALIGNMENT is done, guaranteed. Jobs are done in alphabetical order.
		    # ================================= DONE RUNNING FASTQC ON UNMAPPED (but filtered!) FASTQ FILES ====
		}
	    }
	}
	
	foreach my $sample (keys(%{$remember{$REM_ALIGN_BAM}})) {
	    my $prereqJob = $remember{$REM_ALIGN_BAM}{$sample}{"job"}; defined($prereqJob) or confess "failed to locate the required pre-requisite job variable!";
	    my $bamFile   = $remember{$REM_ALIGN_BAM}{$sample}{"bam"}; defined($bamFile)   or confess "failed to locate the required bam variable!";
	    my $dep = { 'prereqAlignJob' => $prereqJob };

	    # ================================= RUN FASTQC ON MAPPED BAM FILES INDIVIDUALLY ========================
	    my $fqcJobName = sanitizeJobName($cfg->{studyName} . "_s49.3_qc_on_bam_" . $sample);
	    my $fqcVars = {'force'=>$FORCE_RERUN, 'verbose'=>$GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG,
			   'outputDir'   =>  $cfg->{qcAfterMappingDir}, # the PARENT output directory
			   'fastqc'      =>  $cfg->{bin}->{fastqc}    , # binary location
			   'inputFile'   =>  $bamFile                  # <-- Must be the file's FULL PATH
	    };
	    my ($jobF,$qF) = setJobInfo($fqcJobName,$dep,$fqcVars,$cfg->{monkeyPoo},"49.qc.pl");
	    $cfg->{jobs}->{'49_qc_on_mapped'}->{$sample} = { "jobName"=>$jobF, "qsub"=>$qF }; # '49' so it gets alphabetized after all alignment is done, guaranteed.
	    # ================================= DONE RUNNING FASTQC ON MAPPED BAM FILES INDIVIDUALLY ==============
	    # ================================= RUN RSEQC ON MAPPED BAM FILES INDIVIDUALLY ========================
	    my $rJobName = sanitizeJobName($cfg->{studyName} . "_s47_rseqc_" . $sample);
	    my $rseqVars = {   'force'=>$FORCE_RERUN, 'verbose'=>$GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG,
			       'outputDir'       => $cfg->{rseqcQualityDir}        , # the output directory, which will be created
			       'inputBAM'        => $bamFile                       , # <-- Must be the input bam file's FULL PATH
			       'bedAnnot'        => $cfg->{bedAnnotForRSEQC}       , # location of bed-format annotation files for RSEQC
			       'python2_7_exe'   => $cfg->{bin}->{python2_7_exe}   , # binary location for python 2.7 or later
			       'scriptParentDir' => $cfg->{bin}->{rseqc_parent_dir}  # location where all the RSEQC executables are
	    };
	    my ($jobR,$qR) = setJobInfo($rJobName,$dep,$rseqVars,$cfg->{monkeyPoo},"47.rse.qc.pl"); # RSEQC is the name of the set of programs
	    $cfg->{jobs}->{'47_rseqc_on_bam'}->{$sample} = { "jobName"=>$jobR, "qsub"=>$qR }; # '47' so it gets alphabetized after all alignment is done, guaranteed.
	    # ================================= DONE RUNNING RSEQC ON MAPPED BAM FILES INDIVIDUALLY ===============
	}
	
	# ================================= RUN RSEQC ON ALL MAPPED BAM FILES **TOGETHER** ========================
	my $dependOnAllBams; # hash pointer
	foreach my $sample (keys(%{$remember{$REM_ALIGN_BAM}})) {
	    my $prereqJob = $remember{$REM_ALIGN_BAM}{$sample}{"job"}; defined($prereqJob) or confess "failed to locate the required pre-requisite job variable!";
	    $dependOnAllBams->{$prereqJob} = $prereqJob; # save this as a prerequisite!
	}
	my $multiJobName = sanitizeJobName($cfg->{studyName} . "_s46_qc_all_bam");
	my $multiVars = {  'force'=>$FORCE_RERUN, 'verbose'=>$GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG,
			   'outputDir'       => $cfg->{rseqcQualityDir}        , # the output directory, which will be created
			   'inputBamDir'     => $cfg->{mappingResultsDir}      , # we specify a DIRECTORY of bam files
			   'bedAnnot'        => $cfg->{bedAnnotForRSEQC}       , # location of bed-format annotation files for RSEQC
			   'python2_7_exe'   => $cfg->{bin}->{python2_7_exe}   , # binary location for python 2.7 or later
			   'scriptParentDir' => $cfg->{bin}->{rseqc_parent_dir}  # location where all the RSEQC executables are
	};
	my ($jobMulti,$qMulti) = setJobInfo($multiJobName,$dependOnAllBams,$multiVars,$cfg->{monkeyPoo},"46.qc.multiple.bam.pl"); # RSEQC is the name of the set of programs
	$cfg->{jobs}->{'46_qc_multi_bam'}->{"everything"} = { "jobName"=>$jobMulti, "qsub"=>$qMulti }; # '47' so it gets alphabetized after all alignment is done, gua
	# ================================= DONE RUNNING RSEQC ON ALL MAPPED BAM FILES TOGETHER ========================
    } # Done with FastQC quality control
    # put further analyses here...


    
    if ($cfg->{doExpression}) {
	# ================================= [START] RUN CUFFDIFF ON ALL MAPPED BAM FILES **TOGETHER** ========================
	# Depends on: *all* BAM files being aligned (in step 04a), and nothing else.
	my $bamsByGroup = ""; # <-- Will eventually look something like "exp1.bam,exp2.bam  ctrl1.bam,ctrl2.bam,ctrl3.bam"
	foreach my $category ( keys( %{$remember{$REM_EXPERIMENT_CATEGORY}} ) ) {
	    my @bamArray = ();
	    foreach my $replicate (keys(%{$remember{$REM_EXPERIMENT_CATEGORY}{$category}})) {
		my $bamPath = $remember{$REM_EXPERIMENT_CATEGORY}{$category}{$replicate}{"bam"};
		push(@bamArray, $bamPath);
	    }
	    $bamsByGroup .= ((length($bamsByGroup) == 0) ? '' : ' ') . join('<COMMA>', @bamArray); # <-- add a space between categoryes for every category EXCEPT the very first one!. The FINAL result will look something like "exp1.bam,exp2.bam  ctrl1.bam,ctrl2.bam,ctrl3.bam"
	}
	my $jname              = sanitizeJobName($cfg->{studyName} . "_s22_cuffdiff_all");
	my $dependOnAlignments = getHashPointerOfAllBamJobs();
	my $labelsCommaDelim = join("<COMMA>", keys(%{$remember{$REM_EXPERIMENT_CATEGORY}})); # Should look something like: "wildtype<COMMA>control<COMMA>drug<COMMA>whatever". Note that we can't use a literal comma (,), because qsub chews it up. So we use "<COMMA>" instead. The cuffdiff script knows how to deal with this.
	my $vars = { 'force'=>$FORCE_RERUN, 'verbose'=>$GLOBAL_VERBOSE, 'debug'=>$GLOBAL_DEBUG,
		     'cuffdiffExe'     => $cfg->{bin}->{cuffdiff}  ,
		     'outputDir'       => $cfg->{cuffdiffDir}      ,  # <-- the output directory, which will be created
		     'gtfFile'         => $cfg->{gtfFile}          ,
		     'genomeFasta'     => $cfg->{genomeFasta}      ,  # <-- it's OK (although not ideal) if the genomeFasta is undefined
		     'labels'          => $labelsCommaDelim        ,  # <-- labels is delimited by "<COMMA>" and spaces, NOT the actual comma ',', because of qsub weirdness
		     'bamFilesByGroup' => $bamsByGroup                # <-- bamFilesByGroup is delimited by "<COMMA>" and spaces, NOT the actual comma ',', because of qsub weirdness
	};
	my ($jjj,$qqq) = setJobInfo($jname,$dependOnAlignments,$vars,$cfg->{monkeyPoo},"22.cuffdiff.pl"); # RSEQC is the name of the set of programs
	$cfg->{jobs}->{'22_cuffdiff'}->{'everything'} = { "jobName"=>$jjj, "qsub"=>$qqq }; # '47' so it gets alphabetized after all alignment is done, gua
	# ================================= [DONE] RUN CUFFDIFF ON ALL MAPPED BAM FILES **TOGETHER** ========================
    }

}

#########################################################
### 2. generate writeup #################################
#########################################################

sub generateWriteup {
    # use job list to fill out a legible writeup
    my ($cfg) = @_;

    #printf STDERR "Printing writeup (not implemented yet)\n";

    return;
}

#########################################################
### 3. run jobs #########################################
#########################################################

sub runJobs { # Generate AND RUN the qsub "list of jobs to submit" file.
    my ($cfg)   = @_;
    my $jobs    = $cfg->{jobs};
    my $outfile = $cfg->{studyName} . "_jobs.sh"; # Name of the shell script that will we run below
    open (OF,">$outfile") or confess "Cannot write to the following 'qsub job list' file: \"$outfile\"";
    printf OF "#!/bin/bash\n";
    printf OF "set -e    # <-- abort if any 'qsub' command returns a non-zero value\n";
    printf OF "set -u    # <-- abort on undefined variables\n";
    foreach my $stepName (sort keys %$jobs) { # <-- it is CRITICALLY IMPORTANT that the jobs are added to the script in sorted order!
	foreach my $sampleName (sort keys %{$jobs->{$stepName}} ) {
	    print STDERR "Handling sample '$sampleName': step '$stepName'...\n";
	    printf OF  $jobs->{$stepName}->{$sampleName}->{qsub} . "\n"; # Each job has a 'qsub' command associated with it, which we finally print here.
	}
    }
    close(OF);
    verboseSystem("chmod 711 $outfile"); # Allow the script to be executed by this user...
    print STDERR "\n";
    if ($GLOBAL_DRY_RUN) { print STDERR "Dry-run: so we are not ACTUALLY submitting any jobs. Check the generated file \"$outfile\" to see what jobs would have been submitted.\n"; }
    else { print STDERR "Now submitting jobs in the file \"$outfile\" to the Torque job-management queue...\n"; }
    print STDERR "***********************************************************************\n";
    print STDERR "* Note: 1. Type 'qstat -s' to see your job status.\n";
    print STDERR "*       2. To cancel ALL your jobs--including other jobs--you you can use 'qdel all'---beware!\n";
    print STDERR "***********************************************************************\n";
    verboseSystem("bash ./$outfile"); # <-- This actually RUNS the jobs
    ($GLOBAL_VERBOSE && !$GLOBAL_DRY_RUN) && print STDERR "Here are the first 10 jobs (qstat -s):\n" . `qstat -s | head -n 10` . "\n";
}

1; # <-- Perl requires any module to end with a true value. Hence, this 1.

