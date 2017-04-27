#!/usr/bin/perl
# bananas.pm version 2 - parallelized for rigel
# Authors: Sean Thomas and Alex Williams, 2013
# This program is designed to automate many of the tasks that accompany next generation sequencing tasks. 

package bananas;

use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use Scalar::Util qw(looks_like_number);
use File::Spec::Functions qw(catfile); # For constructing filepaths ("catfile" is the function name)
use IO::Compress::Gzip qw(gzip $GzipError) ;
use List::Util qw(max);
use Carp; # backtrace on errors. Has the "confess" function. Use this instead of "die" if you want useful information!
require Exporter;

### 0.i Export functions
our (@ISA, @EXPORT, @EXPORT_OK);
@ISA = qw(Exporter);
@EXPORT_OK = qw(GLOBAL_VERBOSE GLOBAL_DRY_RUN);
@EXPORT = qw(loadConfig checkConfig buildJobSubmissionScript generateWriteup runJobSubmissionScript);

### 0.ii Option choices
our ($GLOBAL_VERBOSE) = 0; # ours to allow access to monkey, 0:quiet,1:verbose
our ($GLOBAL_DRY_RUN) = 0; # 0:fullRun,1:dryRun
our ($GLOBAL_COLOR) = 0; # 0 = only black and white messages. 1 = COLOR text to console
our ($FORCE_RERUN)    = 0; # 0:checkStatus,1:forceRerun

### 0.iii List recognized key names for the config file. Anything that is NOT in here will raise a warning.
my @OK_KEY_NAMES = qw(studyName sampleDir resultDir forceRerun minMapQ libraryAdapterFile bowtie2Index genomicBins genomicBinsID chromSizesFile geneWindows exonFile transcriptFile symbolXref genomeFasta repeatMask monkeyPoo aligner gtfFile tracksURL tracksSCP skipFastQC skipBrowser skipWindows);
my @OK_ALIGNERS  = qw(bowtie tophat);

#########################################################
### 0. Configuration File Functions #####################
#########################################################

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
    if (exists(${$cfgPointer}->{$key})) {
        die "\nDUPLICATE ASSIGNMENT OF '$key' DETECTED (filename: \"$configFileName\")\n>>>> The '$key' was defined a SECOND time on line $lineNum.\n";
    }
    if (not ("$key" ~~ @OK_KEY_NAMES)) { 
        warn "\n>>>> KEY WARNING: the key '$key' (on line $lineNum in the config file $configFileName) was NOT FOUND in the list of recognized keys."
    }
    ${$cfgPointer}->{$key} = $value;
}

# Same as regular "print" in most regards, except it ONLY actually prints if the $GLOBAL_VERBOSE flag is set. Only accepts one string argument. No arrays!
sub verbosePrint($;$) { # If it is a dry run, then it prepends [DRY RUN] to the statement.
    my ($text, $color) = @_;
    if ($GLOBAL_VERBOSE) {
	print STDERR "$text"; #printColorStderr("[DRY RUN]", "yellow on_red"); printColorStderr(" "); }
	#printColorStderr($text, $color);  # Prints to STDERR always!
    }
}

# Runs a system command, and prints it out if we are using "verbose" mode (i.e., monkey.pl --verbose)
# DRY RUN behavior: Does NOT actually run the system call if this is a dry run!
sub verboseSystem($) { 
    my ($cmd) = @_;
    verbosePrint(">>>>> System call >>>>> $cmd\n", "green");
    if ($GLOBAL_DRY_RUN) { return 0; } # Dry run ALWAYS "succeeds," (exit status 0) even though nothing is really run
    else {                 return(system($cmd)); } # <-- Mandatory that we RETURN the system call result!
}

sub loadConfig {
    # loads configuration file information into $cfg
    my ($file) = @_;
    my $cfg;
    # automatically detected skipping options
    $cfg->{skipPeaks}      = 1;      
    $cfg->{skipExpression} = 1;      
    $cfg->{skipExo}        = 1;    
    my $lineNum = 0;
    my $sampleNamesHash;
    open(CF,$file) or die "ERROR: INVALID CONFIG FILE SPECIFIED: Cannot open the configuration file \"$file\" !";
    while (my $line = <CF>) { 
	chomp($line);
	$lineNum++;
	next if ($line =~ /^#/);      
	next if ($line =~ /^(\s)*$/); 
	my @t = split(/\t/,$line);
	if (scalar @t == 1) {
	    my $item = $t[0];
	        if ($item =~ /([^=]+)=(.*)/) { # look for something with an EQUALS SIGN between it
		    my $key   = trimWhitespace($1);
		    my $value = trimWhitespace($2);
		    addKeyValuePairToConfig(\$cfg, $key, $value, $lineNum, $file);
		} else {
		    die "Error in specification file \"$file\" on line $lineNum: Currently, if a line only has ONE element on it, it must be of the form KEY = VALUE.\n";
		}
	} elsif (scalar(@t) == 2) { # look for any line with TWO ELEMENTS on it (i.e. tab-delimited key and value)
	    my $key = $t[0]; my $value = $t[1];
	    addKeyValuePairToConfig(\$cfg, $key, $value, $lineNum, $file);
	} elsif ((scalar @t == 5) or (scalar(@t) == 6)) {  
	    my $idName = $t[1];
	        if (defined($sampleNamesHash->{$idName})) { die "study design error: multiple samples with same id - $idName" } else { $sampleNamesHash->{$idName} = 1}
	        if ($idName =~ /.*[.].*[.].*/) {
		    die "INPUT FILE ID LINE ERROR in <$file> on line $lineNum: The sample ID CANNOT have more than one decimal point in it: <$idName>!\n";
		}
	        if (!($idName =~ /.*[.][0-9]+/)) {
		    die "INPUT FILENAME ERROR in file <$file> on line $lineNum: The name \"$idName\" DID NOT end in a period and then a replicate number. Examples of valid names: WT.1 WT.2 DRUG.1 DRUG.2 DRUG.3\n";
		}

	    my @id_replicate_split = split(/\./,$idName);
	    if (scalar(@id_replicate_split) != 2) {
		die "INPUT FILENAME ERROR in file <$file> on line $lineNum. The sample/input name must be of the format 'name.replicateNumber'. The offending name was \"$idName\" . An example of a good name is \"WILDTYPE.2\" or \"DRUG.1\". Please start counting replicates from ONE!\n";
	    }
	    $cfg->{details}->{$idName}->{name}      = $id_replicate_split[0];
	    $cfg->{details}->{$idName}->{replicate} = $id_replicate_split[1];
	        
	    # get data type <chip|rna|exo>
	    $cfg->{details}->{$idName}->{type} = $t[2];
	    if    (lc($t[2]) eq "chip") {   # chip-seq data
		$cfg->{skipPeaks} = 0;
	    } elsif (lc($t[2]) eq "rna") {  # rna-seq data
		$cfg->{skipExpression} = 0;
	    } elsif (lc($t[2]) eq "exo") {  # ChIP-exo data
		$cfg->{skipPeaks} = 0;
		$cfg->{skipExo}   = 0;
	    } else { die "INPUT SPECIFICATION ERROR in file <$file> on line $lineNum: the sample type (which must be either 'chip', 'exo', or 'rna') is unrecognized. The invalid type was: \"$t[2]\"\n";
	    }
	    
	    if (uc($t[3]) ne "NA" and $t[3] ne '') { # If this field isn't blank or NA...
		$cfg->{details}->{$idName}->{input} = $t[3];
	    }

	    my $isPairedEnd   = (scalar(@t) == 6); # if there are SIX items on this line, it means it's a paired end sample.
	    my $firstPairFile = $t[4];
	        
	    $cfg->{$t[0]}->{$idName}->{$firstPairFile} = 1; # indicate that we have a sample with this specific name
	        
	    if ($isPairedEnd) {
		my $secondPairFile = $t[5];
		$cfg->{$t[0]}->{$idName}->{$secondPairFile} = 1; # indicate that we have a paired end file!!!
	    }
	} else {
	    die "Configuration file must have 2,3, or 4 tab-delimited entries per line. Or 5 for paired-end data. Error on line <$lineNum> of file <$file>.\n";
	}
    }
    close(CF);

    $cfg->{writeUpFile} = catfile($cfg->{resultDir}, "00.writeup.txt");

    return($cfg);
}

sub isBooleanStringOrMissing {
    # Returns true if and only if the input is TRUE/true, FALSE/false, or undefined. Used for checking validity of inputs.
    my ($value) = @_;
    return (!defined($value) || (uc($value) eq "FALSE") || (uc($value) eq "TRUE"));
}

sub trimWhitespace { my $str=shift; $str =~ s/^\s+|\s+$//g; return($str) }; # trims leading and trailing whitespace

sub checkConfig {
    my ($cfg) = @_;

    (defined($cfg->{monkeyPoo})) or die "Configuration file does not contain a monkeyPoo entry\n";
    (-d $cfg->{monkeyPoo})       or die "monkeyPoo directory does not exist. The directory specified was: $cfg->{monkeyPoo}";

    (defined($cfg->{studyName}))        or die "Configuration file must contain a 'studyName' entry.";
    ($cfg->{studyName} !~ /[\/\\\s\t\-\,\.]/) or die "Study name cannot contain whitespaces/slashes/backslashes/hyphens/periods. You can use underscores. The offending name was: \"$cfg->{studyName}\"\n";
    ($cfg->{studyName} !~ /^\d/) or die "Study name cannot begin with a number: $cfg->{studyName}";

    (defined($cfg->{forceRerun})) or die "Config file does not contain a 'forceRerun' entry. You MUST specify either forceRerun=TRUE or forceRerun=FALSE.\n";
    if (lc($cfg->{forceRerun}) eq "false")   { $cfg->{forceRerun} = 0; }    # 0 is a "false" value for perl
    elsif (lc($cfg->{forceRerun}) eq "true") {  $cfg->{forceRerun} = 1; } # 1 is a "true" value for perl
    else {  die "forceRerun entry (which you specified as \"$cfg->{forceRerun}\") must be either TRUE or FALSE! You specified something else.\n"; }

    (defined($cfg->{sampleDir}))  or die "Configuration file must contain a 'sampleDir' entry!";
    ($cfg->{sampleDir} =~ /^\//)  or die "Sample directory must be a FULL PATHNAME (i.e., /path/to/place) and NOT a relative path! The offending name was: $cfg->{sampleDir}!"; # check to make sure it starts with a '/' (and is therefore probably a full path)
    ($cfg->{sampleDir} !~ /[\s]/) or die "Sample directory must NOT contain whitespaces. The offending name was: $cfg->{sampleDir}!";
    (-d $cfg->{sampleDir})        or die "Sample directory must ALREADY exist, but this directory appears not to: $cfg->{sampleDir}!";

    (defined($cfg->{resultDir}))  or die "Configuration file must contain a 'resultDir' entry!";
    ($cfg->{resultDir} =~ /^\//)  or die "Result directory must be a FULL PATHNAME (i.e., /path/to/place) and NOT a relative path! The offending name was: $cfg->{resultDir}!"; # check to make sure it starts with a '/' (and is therefore probably a full path)
    ($cfg->{resultDir} !~ /[\s]/) or die "Working directory cannot contain whitespaces. The offending name was: $cfg->{resultDir}";
    ((-d $cfg->{resultDir}) or mkdir($cfg->{resultDir})) or die "Unable to find or create working directory, check your permissions: $cfg->{resultDir}";

    $cfg->{filterResultsDir} = catfile($cfg->{resultDir},"02.sequenceFilter");
    $cfg->{fastqcResultsDir} = catfile($cfg->{resultDir},"03.fastqc");
    $cfg->{mappingResultsDir} = catfile($cfg->{resultDir},"04.mapping");
    $cfg->{tagResultsDir} = catfile($cfg->{resultDir},"05.tags");
    $cfg->{densityResultsDir} = catfile($cfg->{resultDir},"06.tagDensity");
    $cfg->{peakResultsDir} = catfile($cfg->{resultDir},"07.peaks");
    $cfg->{browserResultsDir} = catfile($cfg->{resultDir},"08.browserTracks");
    $cfg->{windowResultsDir} = catfile($cfg->{resultDir},"09.windowFigures");
    $cfg->{motifDiscDir} = catfile($cfg->{resultDir},"10.motifDiscovery");
    
    # all tools must include full hardcoded paths here, ideally in the binDir for simplifying this section's readability
    my $binDir = "/data/tools/bin";
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
    my $th = $cfg->{bin};
    foreach my $bk (keys %$th) {
	(-e ($th->{$bk})) or die "executable $th->{$bk} does not exist!";
    }
    
    (defined($cfg->{libraryAdapterFile}))              or die "MISSING CONFIG FILE OPTION: Configuration file must contain a 'libraryAdapterFile' entry!";
    (open(TF,$cfg->{libraryAdapterFile}) && close(TF)) or die "Cannot find the following libraryAdapterFile: $cfg->{libraryAdapterFile}!";

    (defined($cfg->{geneWindows}))              or die "MISSING CONFIG FILE OPTION: Configuration file must contain a 'geneWindows' entry!";
    (open(TF,$cfg->{geneWindows}) && close(TF)) or die "Cannot find the following geneWindows file: $cfg->{geneWindows}";

    (defined($cfg->{minMapQ}))       or die "MISSING CONFIG FILE OPTION: Configuration file must contain a 'minMapQ' entry! This is the minimum map quality (MAPQ) that a read must posses in order to NOT be disqualified after alignment. It should be between 0 (meaning 'keep everything') and 100. A standard value is 30.";
    ($cfg->{minMapQ} >= 0 && $cfg->{minMapQ} <= 100) or die "The 'minMapQ' option in the config file must be between 0 and 100, inclusive. The invalid map quality score specified in the file was: $cfg->{minMapQ}\n";

    # Added by Alex June 25, 2014
    isBooleanStringOrMissing($cfg->{skipFastQC})  or die "Configuration file problem: skipFastQC was defined as \"$cfg->{skipFastQC}\", but it must be the literal text TRUE, FALSE, or omitted entirely. Blank values are NOT ACCEPTABLE!";
    isBooleanStringOrMissing($cfg->{skipBrowser}) or die "Configuration file problem: skipBrowser was defined as \"$cfg->{skipBrowser}\", but it must be the literal text TRUE, FALSE, or omitted entirely. Blank values are NOT ACCEPTABLE!";
    isBooleanStringOrMissing($cfg->{skipWindows}) or die "Configuration file problem: skipWindows was defined as \"$cfg->{skipWindows}\", but it must be the literal text TRUE, FALSE, or omitted entirely. Blank values are NOT ACCEPTABLE!";

    # ======================== CHECK THE ALIGNER ====================
    (exists($cfg->{aligner})) or die "MISSING CONFIG FILE OPTION: The configuration file must specify an ALIGNER (the 'aligner') option, which was not specified!";
    $cfg->{aligner} = lc($cfg->{aligner}); # Convert the aligner name to LOWER CASE!
    ($cfg->{aligner} ~~ @OK_ALIGNERS) or die "The configuration file specified the following UNRECOGNIZED aligner: \"$cfg->{aligner}\". We expect the aligner to be something like 'bowtie' or 'tophat'.\n";
    if (($cfg->{aligner} eq "tophat")) { 
	(exists($cfg->{gtfFile})) or die "MISSING CONFIG FILE OPTION: If you specify 'tophat' as your aligner, you must ALSO specify a valid gtf file with the 'gtfFile' parameter!";
	(-e ($cfg->{gtfFile})) or die "The specified GTF-format annotation file (the 'gtfFile' setting in the config file) '$cfg->{gtfFile}' was NOT found on the filesystem. Tophat requires this file.";
    }
    # ======================== CHECK THE ALIGNER ====================
    
    # check samples, also make sure sample names not duplicated
    if (defined($cfg->{input})) {
        my $th = $cfg->{input};
        foreach my $k (keys %$th) {
	    my $th2 = $th->{$k};
	    foreach my $k2 (keys %$th2) { dieIfFileAccessFails(catfile($cfg->{sampleDir}, $k2), "INPUT FILE UNREADABLE"); }
        }
    }

    (defined($cfg->{sample})) or die "Configuration file does not contain any 'sample' lines. e.g.:\nsample <sampleLabel> <fullPathToSampleFile>!";

    my $thsamp = $cfg->{sample};
    foreach my $k (keys %$thsamp) {
	my $th2 = $thsamp->{$k};
	foreach my $sampkey (keys %$th2) { dieIfFileAccessFails(catfile($cfg->{sampleDir}, $sampkey), "SAMPLE FILE UNREADABLE"); }
	$th2 = $cfg->{details}->{$k};
	if (defined($th2->{input})) {
	    my $tch = $cfg->{input};
	    defined($tch->{$th2->{input}}) or die "Input for sample ($k) lists a matched input ($th2->{input}) that isn't specified in the configuration file.";
	}
    }

    if (defined($cfg->{bowtie2Index})) {
	$cfg->{genomeName} = basename($cfg->{bowtie2Index});
    } else {
        die "Configuration file does not contain a bowtie2Index entry, but that is MANDATORY! Even for tophat alignments!\n"
    }
    
    (open(SF,$cfg->{genomicBins})    && close(SF)) or die "Configuration file does not contain a valid genomicBins file!";
    (open(SF,$cfg->{genomicBinsID})  && close(SF)) or die "Configuration file does not contain a valid genomicBinsID file!";
    (open(SF,$cfg->{chromSizesFile}) && close(SF)) or die "Configuration file does not contain a valid chromSizesFile file!";
    (open(SF,$cfg->{geneWindows})    && close(SF)) or die "Configuration file does not contain a valid geneWindows file!";
    (open(SF,$cfg->{exonFile})       && close(SF)) or die "Configuration file does not contain a valid exonFile file!";
    (open(SF,$cfg->{transcriptFile}) && close(SF)) or die "Configuration file does not contain a valid transcriptFile file!";
    (open(SF,$cfg->{symbolXref})     && close(SF)) or die "Configuration file does not contain a valid symbolXref file, or that file could not be read.!";

    unless ($cfg->{skipPeaks}) {
	(defined($cfg->{genomeFasta}))             or die "Configuration file does not contain a valid 'genomeFasta' entry!";
	(open(SF,$cfg->{genomeFasta})&& close(SF)) or die "Configuration file does have a valid genomeFasta file: $cfg->{genomeFasta}!"; 
	(defined($cfg->{repeatMask}))              or die "Configuration file does not contain a valid 'repeatMask' entry!";
	(open(SF,$cfg->{repeatMask})&& close(SF))  or die "Configuration file does have a valid repeatMask bed file: $cfg->{repeatMask}!";
    }    

    #printf STDERR Dumper $cfg; die;

    # set skipping info
    if (defined($cfg->{skipBrowser}) && (uc($cfg->{skipBrowser}) eq "TRUE")) {
	$cfg->{skipBrowser} = 1; # no browser tracks for us!
    } else {
	$cfg->{skipBrowser} = 0;
	(defined($cfg->{tracksURL})) or die "MISSING CONFIG FILE OPTION: browser track URL basename ('tracksURL') needs to be specified in the config file. Or, if you don't want browser tracks, add 'skipBrowser = TRUE' to your config file.";
	(defined($cfg->{tracksSCP})) or die "MISSING CONFIG FILE OPTION: browser track scp location basename ('tracksSCP') needs to be specified in the config file.";
    }

    $cfg->{skipFastQC}  = ( defined($cfg->{skipFastQC})  && (uc($cfg->{skipFastQC})  eq "TRUE") ); # Set to 1 only if this value is DEFINED and also FALSE. Otherwise, 0.
    $cfg->{skipWindows} = ( defined($cfg->{skipWindows}) && (uc($cfg->{skipWindows}) eq "TRUE") ); # Set to 1 only if this value is DEFINED and also FALSE. Otherwise, 0.
}

#########################################################
### 1. Build Job Submission Function ####################
#########################################################

sub setJobInfo {
    # job name, dependencies, variables, script name
    my ($jobName,$dep,$vars,$binDir,$script) = @_;
    $jobName =~ tr/\./\_/;
    my $qsub = "$jobName=\`qsub -N $jobName ";
    if (scalar keys %$dep > 0) {
	$qsub = $qsub . " -W depend=afterok";
	foreach my $k (sort keys %$dep) {
	    my $v = $dep->{$k};
	    $qsub = $qsub . ":\$$dep->{$k}";
	}
    }
    if (scalar(keys(%$vars)) > 0) {
        $qsub = $qsub . " -v ";
        foreach my $k (sort keys %$vars) {
            $qsub = $qsub . "${k}='$vars->{$k}',";
        }
    }
    $qsub = $qsub . " " . catfile($binDir,$script) . "\`";
    verbosePrint("Verbose diagnostic info: Added a new qsub job:\n$qsub\n\n");
    return($jobName,$qsub);
}

sub buildJobSubmissionList {
    # this function organizes all jobs that will need to be run and creates a script
    # that will be run in order to process all of the jobs in the most efficient way possible

    # to add a new job, specify the base job name, dependencies, variables, script directory, and called script

    my ($cfg) = @_;
    my $sampleHash = $cfg->{sample};
    my $inputHash = "NA";
    if (defined($cfg->{input})) {$inputHash = $cfg->{input}};

    # build filtering jobs
    foreach my $h ($sampleHash,$inputHash) {
	if (uc($h) eq "NA") {next;}
	foreach my $k (keys %$h) {
	    my $tHash  = $h->{$k};
	    my @sKeys = sort keys %$tHash;
	    my $jobName = $cfg->{studyName} . "_step2_" . $k;
	    my $dep;
	    my $vars;
	    $vars->{force} = $FORCE_RERUN;
	    $vars->{fastqMcf} = $cfg->{bin}->{fastqMcf};
	    $vars->{fastqDir} = $cfg->{sampleDir};
	    $vars->{filterDir} = $cfg->{filterResultsDir};
	    $vars->{adapterFile} = $cfg->{libraryAdapterFile};
	    $vars->{inputFile1} = $sKeys[0];
	    if (scalar @sKeys == 2) { $vars->{inputFile2} = $sKeys[1] }
	    my ($jn,$qs) = setJobInfo($jobName,$dep,$vars,$cfg->{monkeyPoo},"02.filter.pl");
	    $cfg->{jobs}->{'02_filter'}->{$k}->{jobName} = $jn;
	    $cfg->{jobs}->{'02_filter'}->{$k}->{qsub} = $qs;
	}
    }

    # build qc jobs (conditional)
    unless ($cfg->{skipFastQC}) {
	foreach my $h ($sampleHash,$inputHash) {
	    if (uc($h) eq "NA") {next;} # skip anything where... whatever thi is is an "NA"
	    foreach my $k (keys %$h) {
                my $tHash  = $h->{$k};
                my @sKeys = sort keys %$tHash;
		my $jobName = $cfg->{studyName} . "_step3_" . $k;
		my $dep;
		$dep->{filterJob} = $cfg->{jobs}->{'02_filter'}->{$k}->{jobName};
		my $vars;
		$vars->{force} = $FORCE_RERUN;
		$vars->{fastqc} = $cfg->{bin}->{fastqc};
		$vars->{filterDir} = $cfg->{filterResultsDir};
		$vars->{fastqcDir} = $cfg->{fastqcResultsDir};
		$vars->{inputFile} = $sKeys[0];
		my ($jn,$qs) = setJobInfo($jobName,$dep,$vars,$cfg->{monkeyPoo},"03.fastqc.pl");
		$cfg->{jobs}->{'03_fastqc'}->{$k}->{jobName} = $jn;
		$cfg->{jobs}->{'03_fastqc'}->{$k}->{qsub} = $qs;
	    }
	}
    }

    # build mapping jobs (map to the genome with tophat or bowtie)
    foreach my $h ($sampleHash,$inputHash) {
	if (uc($h) eq "NA") {next;}
	foreach my $k (keys %$h) {
	    my $tHash   = $h->{$k};
	    my @sKeys   = sort keys %$tHash;
	    my $jobName = $cfg->{studyName} . "_step4_" . $k;
	    my $dep;
	    $dep->{filterJob} = $cfg->{jobs}->{'02_filter'}->{$k}->{jobName};
	    my $vars;
	    $vars->{force} = $FORCE_RERUN;
	    if ($cfg->{aligner} eq "tophat" && $cfg->{details}->{$k}->{type} eq "rna") {
		$vars->{aligner} = $cfg->{bin}->{tophat};
		$vars->{gtfFile} = $cfg->{gtfFile}; # not sure why we ONLY pass this along if it's tophat + rna (although it doesn't get used otherwise, of course)
	    } else {
		$vars->{aligner} = $cfg->{bin}->{bowtie2};
	    }
	    $vars->{bowtie2Index} = $cfg->{bowtie2Index};
	    $vars->{filterDir}    = $cfg->{filterResultsDir};
	    $vars->{mappingDir}   = $cfg->{mappingResultsDir};
	    $vars->{samtools}     = $cfg->{bin}->{samtools};
	    $vars->{minMapQ}      = $cfg->{minMapQ};
	    $vars->{sampleName}   = $k;
	    $vars->{inputFile1}   = $sKeys[0];
	    if (scalar @sKeys == 2) { $vars->{inputFile2} = $sKeys[1]; }
	    my ($jn,$qs) = setJobInfo($jobName,$dep,$vars,$cfg->{monkeyPoo},"04.mapping.pl");
	    $cfg->{jobs}->{'04_mapping'}->{$k}->{jobName} = $jn;
            $cfg->{jobs}->{'04_mapping'}->{$k}->{qsub} = $qs;
	}
    }

    # tags
    foreach my $h ($sampleHash,$inputHash) {
        if (uc($h) eq "NA") {next;}
	foreach my $k (keys %$h) {
	    my $inFile = catfile($cfg->{mappingResultsDir},"${k}_$cfg->{genomeName}_q$cfg->{minMapQ}.bam");
	    my $jobName = $cfg->{studyName} . "_step5_" . $k;
	    $jobName =~ tr/\./\_/;
	    my $dep;
	    $dep->{mappingJob} = $cfg->{jobs}->{'04_mapping'}->{$k}->{jobName};
	    my $vars;
	    $vars->{force} = $FORCE_RERUN;
	    $vars->{genome} = $cfg->{genomeName};
	    $vars->{seqType} = $cfg->{details}->{$k}->{type};
	    $vars->{bam2bed} = $cfg->{bin}->{bam2bed};
	    $vars->{sortBed} = $cfg->{bin}->{sortBed};
	    $vars->{tagsDir} = $cfg->{tagResultsDir};
	    $vars->{inputFile} = $inFile;
	    $vars->{sampleName} = $k;
	    my ($jn,$qs) = setJobInfo($jobName,$dep,$vars,$cfg->{monkeyPoo},"05.tags.pl");
	    $cfg->{jobs}->{'05_tags'}->{$k}->{jobName} = $jn;
            $cfg->{jobs}->{'05_tags'}->{$k}->{qsub} = $qs;
	}
    }

    # Tag mapping stats
    {
	my $jobName = $cfg->{studyName} . "_step5_xSummary";
	$jobName =~ tr/\./\_/;
	my $dep;
	my $th = $cfg->{jobs}->{'05_tags'};
	foreach my $k (sort keys %$th) {
	    $dep->{$k} = $cfg->{jobs}->{'05_tags'}->{$k}->{jobName};
	}
	my $vars;
	$vars->{force} = $FORCE_RERUN;
	$vars->{genome} = $cfg->{genomeName};
	$vars->{minMapQ} = $cfg->{minMapQ};
	$vars->{tagsDir} = $cfg->{tagResultsDir};
	$vars->{mappingDir} = $cfg->{mappingResultsDir};
	my ($jn,$qs) = setJobInfo($jobName,$dep,$vars,$cfg->{monkeyPoo},"05.xSummary.pl");
	$cfg->{jobs}->{'05_xSummary'}->{'xSummary'}->{jobName} = $jn;
	$cfg->{jobs}->{'05_xSummary'}->{'xSummary'}->{qsub} = $qs;
    }

    # Tag density
    foreach my $k (keys %$sampleHash) {
        my $jobName = $cfg->{studyName} . "_step6_" . $k;
	$jobName =~ tr/\./\_/;
        my $dep;
	$dep->{tagsJob} = $cfg->{jobs}->{'05_tags'}->{$k}->{jobName};
	my $vars;
	$vars->{force} = $FORCE_RERUN;
	$vars->{genome} = $cfg->{genomeName};
	$vars->{seqType} = $cfg->{details}->{$k}->{type};
	$vars->{bedmap} = $cfg->{bin}->{bedmap};
	$vars->{binsFile} = $cfg->{genomicBins};
	$vars->{binsFileID} = $cfg->{genomicBinsID};
	$vars->{tagsDir} = $cfg->{tagResultsDir};
	$vars->{densityDir} = $cfg->{densityResultsDir};
	$vars->{sampleName} = $k;
	if (defined($cfg->{details}->{$k}->{input})) {
	    $vars->{inputName} = $cfg->{details}->{$k}->{input};
	    $dep->{inputJob} = $cfg->{jobs}->{'05_tags'}->{$vars->{inputName}}->{jobName};
	}
	my ($jn,$qs) = setJobInfo($jobName,$dep,$vars,$cfg->{monkeyPoo},"06.density.pl");
	$cfg->{jobs}->{'06_density'}->{$k}->{jobName} = $jn;
	$cfg->{jobs}->{'06_density'}->{$k}->{qsub} = $qs;
    }

    ## ChIP/exo specific processing
    unless ($cfg->{skipPeaks}) {
	foreach my $k (keys %$sampleHash) {
	    # Do 'chip' and 'exo' specific processing here...
	    if ($cfg->{details}->{$k}->{type} eq "chip" || $cfg->{details}->{$k}->{type} eq "exo") {
		my $jobName = $cfg->{studyName} . "_step7_" . $k;
		$jobName =~ tr/\./\_/;
		my $dep;
		$dep->{densityJob} = $cfg->{jobs}->{'06_density'}->{$k}->{jobName};
		my $vars;
		$vars->{force}      = $FORCE_RERUN;
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
		$cfg->{jobs}->{'07_peaks'}->{$k}->{jobName} = $jn;
		$cfg->{jobs}->{'07_peaks'}->{$k}->{qsub}    = $qs;
	    }
	}

	# Now generate the peak SUMMARY stats... (this doesn't happen if there are no chip or exo samples...)
	my $jobName = $cfg->{studyName} . "_step7_xSummary";
	my $dep;
	my $th = $cfg->{jobs}->{'07_peaks'};
	foreach my $k (sort keys %$th) {
	    $dep->{$k} = $cfg->{jobs}->{'07_peaks'}->{$k}->{jobName};
	}
	my $vars;
	$vars->{force} = $FORCE_RERUN;
	$vars->{genome} = $cfg->{genomeName};
	$vars->{bedmap} = $cfg->{bin}->{bedmap};
	$vars->{tagsDir} = $cfg->{tagResultsDir};
	$vars->{peaksDir} = $cfg->{peakResultsDir};
	my ($jn,$qs) = setJobInfo($jobName,$dep,$vars,$cfg->{monkeyPoo},"07.xSummary.pl");
	$cfg->{jobs}->{'07_xSummary'}->{'xSummary'}->{jobName} = $jn;
	$cfg->{jobs}->{'07_xSummary'}->{'xSummary'}->{qsub} = $qs;
    }
    ## End ChIP/exo-specific processing

    ## This module generates UCSC-compatible browser tracks:
    unless ($cfg->{skipBrowser}) {
	foreach my $k (keys %$sampleHash) {
            my $jobName = $cfg->{studyName} . "_step8_" . $k;
	    my $dep;
	    $dep->{densityJob} = $cfg->{jobs}->{'06_density'}->{$k}->{jobName};
	    if ($cfg->{details}->{$k}->{type} eq "chip" || $cfg->{details}->{$k}->{type} eq "exo") {
		$dep->{peaksJob} = $cfg->{jobs}->{'07_peaks'}->{$k}->{jobName};
	    }
	    my $vars;
	    $vars->{force} = $FORCE_RERUN;
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
	    $cfg->{jobs}->{'08_browser'}->{$k}->{jobName} = $jn;
            $cfg->{jobs}->{'08_browser'}->{$k}->{qsub} = $qs;
	}

	# After the browser tracks have all been added to the queue, generate the SUMMARY...
        my $jobName = $cfg->{studyName} . "_step8_xSummary";
        my $dep; # hash
        my $th = $cfg->{jobs}->{'08_browser'};
        foreach my $k (sort keys %$th) {
	    $dep->{$k} = $cfg->{jobs}->{'08_browser'}->{$k}->{jobName};
        }
	my $vars;
	$vars->{force}      = $FORCE_RERUN;
	$vars->{scp}        = catfile($cfg->{tracksSCP},$cfg->{studyName});
	$vars->{browserDir} = $cfg->{browserResultsDir};
	my ($jn,$qs) = setJobInfo($jobName,$dep,$vars,$cfg->{monkeyPoo},"08.xSummary.pl");
	$cfg->{jobs}->{'08_xSummary'}->{'xSummary'}->{jobName} = $jn;
	$cfg->{jobs}->{'08_xSummary'}->{'xSummary'}->{qsub}    = $qs;
    }
    ## End browser track module

    ## Generate windows
    unless ($cfg->{skipWindows}) {
	foreach my $k (keys %$sampleHash) {
            my $jobName = $cfg->{studyName} . "_step9_" . $k;
	    my $analyze = catfile($cfg->{monkeyPoo},"09.xAnalyze.RScript");
            my $dep;
	    $dep->{densityJob} = $cfg->{jobs}->{'06_density'}->{$k}->{jobName};
	    my $vars;
	    $vars->{force} = $FORCE_RERUN;
	    $vars->{genome} = $cfg->{genomeName};
	    $vars->{bedmap} = $cfg->{bin}->{bedmap};
	    $vars->{sortBed} = $cfg->{bin}->{sortBed};
	    $vars->{densityDir} = $cfg->{densityResultsDir};
	    $vars->{windowDir} = $cfg->{windowResultsDir};
	    $vars->{sampleName} = $k;
	    $vars->{geneWindows} = $cfg->{geneWindows};
	    $vars->{analyze} = $analyze;
	    my ($jn,$qs) = setJobInfo($jobName,$dep,$vars,$cfg->{monkeyPoo},"09.windows.pl");
	    $cfg->{jobs}->{'09_window'}->{$k}->{jobName} = $jn;
            $cfg->{jobs}->{'09_window'}->{$k}->{qsub} = $qs;
	}
    }
    ## windows module

    ## Motif discovery module (only occurs if there are chip/exo data, and only on chip/exo data)
    unless ($cfg->{skipPeaks}) {
        foreach my $k (keys %$sampleHash) {
            if ($cfg->{details}->{$k}->{type} eq "chip" || $cfg->{details}->{$k}->{type} eq "exo") {
                my $jobName = $cfg->{studyName} . "_step10_" . $k;
                my $dep;
                $dep->{peaksJob} = $cfg->{jobs}->{'07_peaks'}->{$k}->{jobName};
                my $vars;
                $vars->{force} = $FORCE_RERUN;
                $vars->{genome} = $cfg->{genomeName};
                $vars->{bedops} = $cfg->{bin}->{bedops};
                $vars->{bedtools} = $cfg->{bin}->{bedtools};
                $vars->{peakMotifs} = $cfg->{bin}->{peakMotifs};
                $vars->{fasta} = $cfg->{genomeFasta};
                $vars->{repMask} = $cfg->{repeatMask};
                $vars->{infile} = catfile($cfg->{peakResultsDir},"${k}_$cfg->{genomeName}_peaks.bed");
                if ($cfg->{details}->{$k}->{type} eq "exo") {
                    $vars->{infile} = catfile($cfg->{peakResultsDir},"${k}_$cfg->{genomeName}_footprints.bed");
                }
                $vars->{motifsDir} = $cfg->{motifDiscDir};
                $vars->{sampleName} = $k;

                my ($jn,$qs) = setJobInfo($jobName,$dep,$vars,$cfg->{monkeyPoo},"10.motifDiscovery.pl");
                $cfg->{jobs}->{'10_motifDiscovery'}->{$k}->{jobName} = $jn;
                $cfg->{jobs}->{'10_motifDiscovery'}->{$k}->{qsub} = $qs;
            }
        }

        # summarize motif discovery, put all discovered motifs into a single transfac format file
        my $jobName = $cfg->{studyName} . "_step10_xSummary";
        my $dep;
        my $th = $cfg->{jobs}->{'10_motifDiscovery'};
        foreach my $k (sort keys %$th) {
            $dep->{$k} = $cfg->{jobs}->{'10_motifDiscovery'}->{$k}->{jobName};
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
        $cfg->{jobs}->{'10_motifDiscovery'}->{'xSummary'}->{jobName} = $jn;
        $cfg->{jobs}->{'10_motifDiscovery'}->{'xSummary'}->{qsub} = $qs;
    }
    ## end motif discovery module

    # put further analyses here...
    
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

sub runJobs {
    my ($cfg)   = @_;
    my $jobs    = $cfg->{jobs};
    my $outfile = $cfg->{studyName} . "_jobs.sh"; # Name of the shell script that will we run below
    open (OF,">$outfile") or confess "Cannot open the following 'jobs' file for writing: \"$outfile\"";
    printf OF "#!/bin/bash\n";
    #printf STDERR Dumper $jobs; die;
    foreach my $k (sort keys %$jobs) { 
	my $th = $jobs->{$k};
	foreach my $j (sort keys %$th) {
	    printf OF "$th->{$j}->{qsub}\n";
	}
    }
    close(OF);
    verboseSystem("chmod 711 $outfile"); # Allow the script to be executed by this user...
    print STDERR "\n";
    print STDERR "Now submitting jobs in the file \"$outfile\" to the Torque job-management queue...\n";
    print STDERR "***********************************************************************\n";
    print STDERR "* Note: 1. Type 'qstat -s' to see your job status.\n";
    print STDERR "*       2. To cancel ALL your jobs--including other jobs--you you can use 'qdel all'---beware!\n";
    print STDERR "***********************************************************************\n";
    verboseSystem("bash ./$outfile"); # <-- This actually RUNS the jobs
    ($GLOBAL_VERBOSE) && print STDERR "Here are the first 10 jobs (qstat -s):\n" . `qstat -s | head -n 10` . "\n";
}
