#!/usr/bin/perl





my $syslog = "X.07.NucleoAtac.syslog.txt";

#####my $pythonProfile = "/home/sthomas/envs/atac"; # <-------- hard coded... but not actually used anywhere!!!!!!!!!!!!

use strict; use warnings; use Carp; # Carp = "confess" function
use File::Spec::Functions qw(catfile);
use Data::Dumper;
use lib "/wynton/group/gladstone/biocore/MonkeyPipeline"; # To make sure we find bananas_agw
use bananas_agw;
use Cwd;

# check for correct qsub variables needed to perform density calculations
bananas_agw::requireEnvOrDie('sampleName','genome','seqType','bedops','nucleoatac','genomeFasta', 'atacDir','peaksDir','mappingDir','atacPython'); # <-- these vars MUST be defined in the %ENV hash
my $force       = bananas_agw::envLooksTrue("force");
my $verbose     = bananas_agw::envLooksTrue("verbose");
my $debug       = bananas_agw::envLooksTrue("debug");

my $k        = $ENV{'sampleName'};
my $genome   = $ENV{'genome'};
my $natac    = $ENV{'nucleoatac'};    bananas_agw::dieIfFileAccessFails($natac, "[ERROR in 07c.NucleoAtac.pl]: We were not able to find the required NucleoAtac executable.");
my $python   = $ENV{'atacPython'};    bananas_agw::dieIfFileAccessFails($python, "[ERROR in 07c.NucleoAtac.pl]: We were not able to find the required Python executable.");
my $fasta    = $ENV{'genomeFasta'};   bananas_agw::dieIfFileAccessFails($fasta, "[ERROR in 07c.NucleoAtac.pl]: We were not able to find the required input fasta file.");



my $inBed       = catfile($ENV{'peaksDir'},"${k}_${genome}_peaks.bed"); bananas_agw::dieIfFileAccessFails($inBed, "[ERROR in 07c.NucleoAtac.pl]: We were not able to find the input bed file.");
my $inBam       = catfile($ENV{'mappingDir'},"${k}_${genome}_q30.bam"); bananas_agw::dieIfFileAccessFails($inBam, "[ERROR in 07c.NucleoAtac.pl]: We were not able to find the input bam file.");

# make sure files and write directories exist
bananas_agw::mkdirOrDie($ENV{'atacDir'});
chdir("$ENV{'atacDir'}"); ($verbose) && print STDERR "[OK] [07c.NucleoAtac.pl]: Changing directory to \"$ENV{atacDir}\"...\n";

my $thisBase    = "${k}_${genome}";
my $outBGFile   = catfile($ENV{'atacDir'},"${k}_${genome}_nucOcc.bedGraph.gz");
my $outBedFile  = catfile($ENV{'atacDir'},"${k}_${genome}_nucOccPeaks.bed.gz");
my $opBGTmp     = catfile($ENV{'atacDir'},"${k}_${genome}.occ.bedgraph.gz");
my $opBedTmp    = catfile($ENV{'atacDir'},"${k}_${genome}.occpeaks.bed.gz");
my $padBed      = catfile($ENV{'atacDir'},"${k}_${genome}_paddedPeaks.bed");

if (-e $outBGFile and -e $outBedFile and !$force) { # Do the output files ALREADY exist, AND we aren't forcing a re-run? Then skip it.
	($verbose) && print STDOUT "[OK] [07c.NucleoAtac.pl: SKIPPING RE-RUN]: the output already exists at '$outBGFile' and '$outBedFile. No need to re-generate it.";
	exit(0); # Output ALREADY exists, so this is OK---exit without any errors!
}

bananas_agw::dieIfFileAccessFails($ENV{'bedops'});

my $BEDOPS_RANGE = 1500; # Alex does not actually know what this parameter is... (it was set by Sean earlier, so it's probably fine)
bananas_agw::systemAndLog("$ENV{'bedops'} -m --range $BEDOPS_RANGE $inBed > $padBed", $verbose, $syslog);
bananas_agw::systemAndLog("$python $natac run  --bed $padBed  --bam $inBam  --fasta $fasta  --out $thisBase", $verbose, $syslog);
bananas_agw::systemAndLog("/bin/mv -f  $opBGTmp   $outBGFile", $verbose, $syslog); # move the COMPLETED file. We can be pretty sure that the completed file is in one piece, and not partially completed, since we wrote to a temp file before renaming it all at once.
bananas_agw::systemAndLog("/bin/mv -f  $opBedTmp  $outBedFile", $verbose, $syslog); # move the COMPLETED file. We can be pretty sure that the completed file is in one piece, and not partially completed, since we wrote to a temp file before renaming it all at once.
bananas_agw::systemAndLog("/bin/rm $padBed", $verbose, $syslog);

bananas_agw::dieIfFileAccessFails($outBGFile,  "[ERROR] generating the expected final output BG file in 07c.NucleoAtac.pl");
bananas_agw::dieIfFileAccessFails($outBedFile, "[ERROR] generating the expected final output bed file in 07c.NucleoAtac.pl")

