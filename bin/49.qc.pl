#!/usr/bin/perl

# DESCRIPTION: This script runs FastQC for doing quality control on input files
#              (can be both .fq files and .bam files, depending on the input arguments).
# DEPENDS ON: The .fq or .bam files need to exist already.
# Warning: if you have a single quote on a PBS comment line above (e.g., a word like "don't"), it will trigger a "qsub: unmatched" error and exit code 1.

use strict; use warnings; use Carp; # Carp = "confess" function
use File::Spec::Functions qw(catfile);
#use Data::Dumper;
use File::Basename;
use lib "/wynton/group/gladstone/biocore/monkey"; # To make sure we find bananas_agw
use bananas_agw;

sub getOutputPrefixForFastQC($$) { # Alex's function for figuring out what the actual name of the output directory / output HTML file / output ZIP file that FastQC will generate.
	# Returns the SPECIFIC output prefix for this single input FASTQ file. The actual filename would then be something like "yourfasta.html"
	# **** Note that this is NOT THE ARGUMENT you pass to fastqc! Do NOT use this for --outdir! **********
	# **** (The actual --outdir argument should be the PARENT of this directory.)               **********
	# It's a little tricky to figure out exactly what the FastQC output file that we should look for will be named. Also, they keep CHANGING it. For the real answer, you have to look in the FastQC source, specifically in the "OfflineRunner.java" source file.
	# Annoyingly, fastqc will sometimes remove file extensions (for ".fastq" or ".fastq.gz" for example), but other times it won't (it keeps ".fq", even though you would expect it to remove that!)
	# This code below is probably STILL not totally correct.
	# Here is the code from FastQC (in the file "Analysis/OfflineRunner.java") that actually removes the file extensions:
	# String fileName = file.getFile().getName().replaceAll("\\.gz$","").replaceAll("\\.bz2$","").replaceAll("\\.txt$","").replaceAll("\\.fastq$", "").replaceAll("\\.fq$", "").replaceAll("\\.csfastq$", "").replaceAll("\\.sam$", "").replaceAll("\\.bam$", "")+"_fastqc.html";
	# Now you'd think it would remove ".fq" since it says so right there, but empirically that is not the case!
	my ($outDirPath, $inFilename) = @_;
	if ($inFilename =~ m/\//) { confess "What, the input file contains a '/'? It should NOT be a path! Just a filename."; }
	# String fileName = file.getFile().getName().replaceAll("\\.gz$","").replaceAll("\\.bz2$","").replaceAll("\\.txt$","").replaceAll("\\.fastq$", "").replaceAll("\\.fq$", "").replaceAll("\\.csfastq$", "").replaceAll("\\.sam$", "").replaceAll("\\.bam$", "")+"_fastqc.html";
	# remove specific extensions IN ORDER. Probably is CASE SENSITIVE even though it seems like it should not be!! Sorry, that's just how the FastQC code works!
	# *** We are just mimicking exactly what the Java code does in order, which is why this code below is so dumb looking! ***
	$inFilename =~ s/[.](gz)$//; $inFilename =~ s/[.](bz2)$//; $inFilename =~ s/[.](txt)$//; $inFilename =~ s/[.](fastq)$//; $inFilename =~ s/[.](fq)$//; $inFilename =~ s/[.](csfastq)$//; $inFilename =~ s/[.](sam)$//; $inFilename =~ s/[.](bam)$//;
	return catfile($outDirPath, "${inFilename}_fastqc"); # Note: returns the full path WITH THE DIRECTORY
}

bananas_agw::requireEnvOrDie("fastqc","outputDir","inputFile"); # <-- these variables MUST be defined in the %ENV hash
my $force              = bananas_agw::envLooksTrue("force");
my $verbose            = bananas_agw::envLooksTrue("verbose");
my $debug              = bananas_agw::envLooksTrue("debug");
my $inPath             = $ENV{'inputFile'};
my $fastqcOutDir       = $ENV{'outputDir'};
my $ncpu               = 1; # Don't run multithreaded
my $fastqcOutputRawWithExtension =                  catfile($fastqcOutDir, File::Basename::basename($inPath));
my $fastqcOutputBase             = getOutputPrefixForFastQC($fastqcOutDir, File::Basename::basename($inPath)); # Full path with the directory

#($verbose) && print STDERR "[OK] [FastQC: Checking to see if the following files exist already:\n";
#($verbose) && print STDERR "      * FastQC File status: ${fastqcOutputBase}.html " . ((-e "${fastqcOutputBase}.html") ? "EXISTS" : "nonexistent") . "\n";
#($verbose) && print STDERR "      * FastQC File status: ${fastqcOutputBase}.zip  " . ((-e "${fastqcOutputBase}.zip" ) ? "EXISTS" : "nonexistent") . "\n";
#($verbose) && print STDERR "      * FastQC File status: ${fastqcOutputRawWithExtension}.html " . ((-e "${fastqcOutputRawWithExtension}.html") ? "EXISTS" : "nonexistent") . "\n";
#($verbose) && print STDERR "      * FastQC File status: ${fastqcOutputRawWithExtension}.zip  " . ((-e "${fastqcOutputRawWithExtension}.zip" ) ? "EXISTS" : "nonexistent") . "\n";
if (!$force and (-e "${fastqcOutputBase}.html" or -e "${fastqcOutputBase}.zip" or -e "${fastqcOutputRawWithExtension}.html" or -e "${fastqcOutputRawWithExtension}.zip")) {
	# EXIT EARLY IF ANY of the candidate 'possible output file looks like this' already exists.
	# FastQC makes files with an specific method of determining filenames, which we *cannot* set, so we have to hope that one of the above files is correct
	# So we check to see if those files exist by looking for the possible options above.
	($verbose) && print STDERR "[OK] [Skipping FastQC re-run] since the QC output files already exist. (The name we looked for matched a pattern like like ${fastqcOutputBase}*).\n";
	exit(0); # <-- EXIT EARLY
}
($verbose) && print STDERR "[OK] [RUNNING FastQC] to generate the output base file path \"${fastqcOutputBase}\"...\n\n";
(-e $inPath) or confess "The input file <$inPath> (which should probably be a FULL PATH to a BAM or FASTQ file) could not be found! We were about to attempt to run FastQC on it! Check to make sure that this file exists and is a full path (not a relative path). ";
bananas_agw::mkdirOrDie($fastqcOutDir); # Has to come BEFORE symlinking and checking the input files
bananas_agw::systemAndLog("$ENV{'fastqc'}  --quiet --threads=$ncpu --outdir=${fastqcOutDir}  ${inPath}", $verbose); # Run FastQC
