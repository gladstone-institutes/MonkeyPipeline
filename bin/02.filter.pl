#!/usr/bin/perl

my $syslog = "X.02.filter.syslog.txt";

use strict; use warnings; use Carp; # Carp = "confess" function
use File::Spec::Functions qw(catfile);
use Data::Dumper;
use lib "/wynton/group/gladstone/biocore/monkey"; use bananas_agw; # To make sure we find bananas_agw

my $force                = bananas_agw::envLooksTrue("force");
my $verbose              = bananas_agw::envLooksTrue("verbose");
my $debug                = bananas_agw::envLooksTrue("debug"); # extra-verbose
my $shouldFilterAdapters = bananas_agw::envLooksTrue("shouldFilterAdapters"); # If this is TRUE, then we filter. If it is FALSE, then we don't actually do any filtering (we just symlink input -> output VERBATIM, directly)

bananas_agw::requireEnvOrDie('fastqDir','filterDir','inputFile1', 'shouldFilterAdapters', 'downstreamBamFile'); # We always need these variables, even if we *aren't* filtering.

if (!$force and -e $ENV{'downstreamBamFile'}) { # Check to see if the EVENTUALLY CREATED bam file already exists for this set of FASTQs. If so, we will not bother to re-filter unless that final bam file is deleted. This way, the user doesn't need to keep around all these intermediate trimmed files.
	$verbose and print STDERR "[SKIPPING] skipping the generation of any filtered FASTQ files, as the downstream BAM that they generate ($ENV{'downstreamBamFile'}) already exists, so there is no need for these filtered files UNLESS the 'force' option is specified.";
	exit(0); # <-- EXIT EARLY
}

if ($shouldFilterAdapters) {
	bananas_agw::requireEnvOrDie('fastqMcf', 'adapterFile'); # We ALSO need these variables if we are filtering
	(-e $ENV{'fastqMcf'})    or confess "[ERROR] fastq-mcf executable location did not appear to actually exist (check that this is valid: $ENV{'fastqMcf'} ) ...";
	(-e $ENV{'adapterFile'}) or confess "[ERROR] Adapter file <$ENV{'adapterFile'}> appears not to be a valid file! ...";
}
(-d $ENV{'fastqDir'})    or confess "[ERROR] Error in filtering (02.filter.pl) -- the FASTQ input file directory <$ENV{'fastqDir'}> appears not to exist! ";

# check if reads are paired end
my $isPaired = ((not bananas_agw::isNA($ENV{'inputFile2'})) and ($ENV{'inputFile2'} ne "")) ? 1 : 0;
my $fq1      =               catfile($ENV{'fastqDir'}, $ENV{'inputFile1'});
my $fq2      = ($isPaired) ? catfile($ENV{'fastqDir'}, $ENV{'inputFile2'}) : "";


bananas_agw::dieIfFileAccessFails($fq1,                 "[ERROR] Failed to find the specified first-pair input file named <$fq1>. Check to make sure this file exists AND there is a full path to this file AND that the file does not have any special characters or spaces in the name.");
($isPaired) and bananas_agw::dieIfFileAccessFails($fq2, "[ERROR] Failed to find the specified second-pair input file named <$fq2>. Check to make sure this file exists AND there is a full path to this file AND that the file does not have any special characters or spaces in the name.");

# Check for bzipped files and warn the user if they are seen! fastq-mcf can only handle .gz and uncompressed files.
($fq1 !~ m/[.](bz2|bzip2)$/i and $fq2 !~ m/[.](bz2|bzip2)$/i) or confess "[ERROR]: Sorry, although fastq-mcf CAN handle gzipped (.gz) and uncompressed fasta files, but it CANNOT handle .bz2 files. Check <$fq1> and <$fq2> to make sure they are gzipped and not .bz2.";

my ($of1, $of2);   # of1 = output file 1 (the only one, if the inputs are single-end)
my ($tof1, $tof2); # Temporary files where filtering is written to BEFORE it is totally finished. Only FINISHED and properly-run jobs will get 'promoted' from $tof1 -> $of1 and $tof2 -> $of2
my $is_fq1_gzipped = ($fq1 =~ m/[.](gz|gzip)$/i);
my $is_fq2_gzipped = ($fq2 =~ m/[.](gz|gzip)$/i);
if ($shouldFilterAdapters) {
	# When we perform the filtering, we ALWAYS gzip the output!
	$of1 =               catfile($ENV{'filterDir'}, bananas_agw::noExtensions($ENV{'inputFile1'}).".gz");
	$of2 = ($isPaired) ? catfile($ENV{'filterDir'}, bananas_agw::noExtensions($ENV{'inputFile2'}).".gz") : "";

	$tof1 =               catfile($ENV{'filterDir'}, "temp_".bananas_agw::noExtensions($ENV{'inputFile1'}).".gz");
	$tof2 = ($isPaired) ? catfile($ENV{'filterDir'}, "temp_".bananas_agw::noExtensions($ENV{'inputFile2'}).".gz") : "";	
} else {
	# If we aren't actually performing filtering, then we just symlink the original input, which MAY possibly not be gzipped.
	# In that unlikely case, we will just symlink the original filename without the gz suffix.
	# Whatever the case may be, ".gzip" is always changed to ".gz".
	my $fq1_compression_suffix = ($is_fq1_gzipped) ? ".gz" : ""; # for the OUTPUT file
	my $fq2_compression_suffix = ($is_fq2_gzipped) ? ".gz" : ""; # for the OUTPUT file
	$of1 =               catfile($ENV{'filterDir'}, bananas_agw::noExtensions($ENV{'inputFile1'}.$fq1_compression_suffix));
	$of2 = ($isPaired) ? catfile($ENV{'filterDir'}, bananas_agw::noExtensions($ENV{'inputFile2'}.$fq2_compression_suffix)) : "";
	
	$tof1 =               catfile($ENV{'filterDir'}, "temp_".bananas_agw::noExtensions($ENV{'inputFile1'}.$fq1_compression_suffix));
	$tof2 = ($isPaired) ? catfile($ENV{'filterDir'}, "temp_".bananas_agw::noExtensions($ENV{'inputFile2'}.$fq2_compression_suffix)) : "";
}

# ================= EXIT EARLY IF THE OUTPUT FILES ALREADY EXIST ==================
my $resultFilesExistAlready = (-e "${of1}") and (!$isPaired or -e "${of2}");
if (not $force and $resultFilesExistAlready) { # Do the file(s) ALREADY exist AND we aren't forcing a re-run? Then skip it/them.
	(!$isPaired) && $verbose && print STDERR    "[OK] SKIPPING re-filtering of single-end input file <$fq1>, because the output file <$of1> already exists.\n";
	($isPaired)  && $verbose && print STDERR (qq{[OK] SKIPPING re-filtering of paired-end input files <$fq1> (left half)...\n} .
						  qq{                                              ...and <$fq2> (right half), because both output files (<$of1> and <$of2>) already exist.\n});
	exit(0);
}

# ================= IF WE GET TO THIS POINT, WE NEED TO GENERATE THE OUTPUT FILES ($of1 and $of2) =========
if ($is_fq1_gzipped) {
	my $code1 = bananas_agw::systemAndLog(qq{gzip --test "$fq1"}, $verbose); # Check that the gzipped file is valid!
	if (0 != $code1) {
		bananas_agw::systemAndLog(qq{touch "${of1}_FAILED--the_gzipped_input_file_was_corrupt.fix_me"}, $verbose); # Check that the gzipped file is valid!
		confess qq{[ERROR] in 02.filter.pl: The GZIPPED input file "$fq1" was NOT a valid gzip archive! It failed to verify when we ran 'gzip --test' on it. Cannot proceed. Fix this by fixing the source file! };
	}
}
if ($is_fq2_gzipped) {
	my $code2 = bananas_agw::systemAndLog(qq{gzip --test "$fq2"}, $verbose); # Check that the gzipped file is valid!
	if (0 != $code2) {
		bananas_agw::systemAndLog(qq{touch "${of2}_FAILED--the_gzipped_input_file_was_corrupt.fix_me"}, $verbose); # Check that the gzipped file is valid!
		confess qq{[ERROR] in 02.filter.pl: The GZIPPED input file "$fq2" was NOT a valid gzip archive! It failed to verify when we ran 'gzip --test' on it. Cannot proceed. Fix this by fixing the source file! };
	}
}

bananas_agw::mkdirOrDie($ENV{'filterDir'});

if ($shouldFilterAdapters) {
	# Note: fastq-mcf uses the DEFAULT quality filtering threshold if you don't give it one.
	#      *  For example, "-q 10" is the default if not overridden.
	# We have, as you can see here, overridden it with the "$phredQualThresh" variable to allow us to change this in the future.
        my $phredQualThresh = 10;
	my $logName = ($isPaired) ? "${of1}.paired.fastqmcf_log.txt" : "${of1}.fastqmcf_log.txt"; # Always the first file, even when paired-end
	my $cmd = "need to set this";
	if ($isPaired) { $cmd = qq{$ENV{'fastqMcf'} -q $phredQualThresh -o $tof1 -o $tof2 $ENV{'adapterFile'} $fq1 $fq2 > $logName}; }
	else {           $cmd = qq{$ENV{'fastqMcf'} -q $phredQualThresh -o $tof1          $ENV{'adapterFile'} $fq1      > $logName}; }

	my $exitCode = bananas_agw::systemAndLog($cmd, $verbose);
	($exitCode == 0) or confess("02.filter.pl: Filtering step failure -- exit code of fastq-mcf command must be zero! The command exited with code $exitCode, and the full command was: $cmd ");

	# Ok, since the files were apparently generated correctly above, NOW we can move them from "$tof" (temp out file) to "$of" (out file)
	bananas_agw::systemAndLog(                "/bin/mv -f $tof1 $of1", $verbose);
	($isPaired) and bananas_agw::systemAndLog("/bin/mv -f $tof2 $of2", $verbose);
} else {
	# If we skip the filtering step, we just make SYMLINKS to the original input fastq files.
	# This way, downstream analysis will still find the ""filtered"" (but not really) samples where it expects them to be.
	bananas_agw::systemAndLog(                "/bin/ln -s $fq1 $of1", $verbose); # Don't really filter, just symlink!
	($isPaired) and bananas_agw::systemAndLog("/bin/ln -s $fq2 $of2", $verbose); # Don't really filter, just symlink!
}

bananas_agw::dieIfFileAccessFails(                $of1, "[ERROR]: Failed to generate filtered 1st-end-of-the-pair output file <${of1}>");
(!$isPaired) or bananas_agw::dieIfFileAccessFails($of2, "[ERROR]: Failed to generate filtered 2nd-end-of-the-pair output file <${of2}>");

