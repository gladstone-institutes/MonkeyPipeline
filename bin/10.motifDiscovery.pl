#!/usr/bin/perl





use strict; use warnings; use Carp; # Carp has "confess"
use File::Spec::Functions qw(catfile);
use Data::Dumper;
use lib "/wynton/group/gladstone/biocore/monkey"; use bananas_agw; # To make sure we find bananas_agw

bananas_agw::requireEnvOrDie('genome','bedops','bedtools','peakMotifs','bedtools','fasta','repMask','infile','motifsDir','sampleName'); # <-- these vars MUST be defined in the %ENV hash
# Check that the ZLM_LICENSE environment variable points to a valid file.
(defined($ENV{'ZLM_LICENSE'})) or confess "ERROR: The environment variable 'ZLM_LICENSE'---which should point to a zlm license file (probably named 'vmatch.lic'---needs to be set. This environment variable will be used by peakMotifs.";
(-e $ENV{'ZLM_LICENSE'})       or confess "ERROR: the ZLM_LICENSE environment variable, which is needed by <peakMotifs> and points to a file (probably named 'vmatch.lic') on the filesystem, was set to the seemingly nonexistent (or non-readable) file <$ENV{'ZLM_LICENSE'}>. You should double check to make sure this file exists and is readable by the current user.";
(-r $ENV{'ZLM_LICENSE'})       or confess "ERROR: the ZLM_LICENSE file appears to exist, but is NOT READABLE by the current user! The file path is <$ENV{'ZLM_LICENSE'}>. Fix this with chmod. It is also possible that the ENCLOSING folder is not readable by this user.";
my $force       = bananas_agw::envLooksTrue("force");
my $verbose     = bananas_agw::envLooksTrue("verbose");

# make sure write directory exists or make it
bananas_agw::mkdirOrDie($ENV{'motifsDir'});
# check existence of input file
(-e $ENV{'infile'}) or confess "unable to find input bed file for motif finding. Offending input file name was: $ENV{'infile'}";

# set output file names
my $repMasked = catfile($ENV{'motifsDir'},"$ENV{'sampleName'}_$ENV{'genome'}_pks_repMask.bed");
my $fastaFile = catfile($ENV{'motifsDir'},"$ENV{'sampleName'}_$ENV{'genome'}_pks_repMask.fasta");

# generate fasta file if needed
if ($force or ! -e $fastaFile) {
	bananas_agw::systemAndLog(qq{$ENV{'bedops'}   -n -1 $ENV{'infile'} $ENV{'repMask'}   > $repMasked }, $verbose);
	bananas_agw::systemAndLog(qq{$ENV{'bedtools'} getfasta -fi $ENV{'fasta'} -bed $repMasked  >  $fastaFile}, $verbose);  # previously was "-fo OUTPUT" --- now it goes to STDOUT as of version 2.25
	unlink($repMasked);
}

# perform motif discovery if needed
my $thisMotifDir = catfile($ENV{'motifsDir'},"$ENV{'sampleName'}_pmResults");
if ($force or ! -d $thisMotifDir) {
	# Note: ZLM_LICENSE must be defined and point to a valid "vmatch.lic" file. Example: /home/sthomas/tools/vmatch-2.2.4-Linux_x86_64-64bit/vmatch.lic
	bananas_agw::mkdirOrDie($thisMotifDir);
	bananas_agw::systemAndLog("$ENV{'peakMotifs'} -i $fastaFile -prefix $ENV{'sampleName'} -outdir $thisMotifDir", $verbose);
}
