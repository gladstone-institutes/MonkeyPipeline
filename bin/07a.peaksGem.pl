#!/usr/bin/perl

my $syslog = "X.07.gemPeaks.syslog.txt";

use strict; use warnings; use Carp; # Carp = "confess" function
use File::Spec::Functions qw(catfile);   use Data::Dumper;   use File::Basename;
use lib "/wynton/group/gladstone/biocore/MonkeyPipeline"; use bananas_agw; # To make sure we find bananas_agw
#use File::Path qw(make_path remove_tree);
use Cwd;


# ================ REPLICATES are currently not handled in any fashion, but GEM can handle them if desired. We'd have to add this to monkey somehow.
# Multi-condition v.s. Multi-replicates  (http://groups.csail.mit.edu/cgs/gem/#QA)
# GEM can analyze binding data from multiple conditions (time points) simultaneously. The user need to give them different names, for example, -–exptCond1 CTCF_cond1.bed -–exptCond2 CTCF_cond2.bed.
# For multiple replicates of same condition, you can specify multiple replicates as separate files, for example, -–exptCond1 CTCF_cond1_rep1.bed -–exptCond1 CTCF_cond1_rep2.bed (note that they need to have the same name). GEM will combine the replicates as one large dataset for analysis. 
#my $tagDensityThreshold = 100;

# check for correct qsub variables needed to perform density calculations
bananas_agw::requireEnvOrDie('sampleName','genome','seqType','gemJar','readDist','expBam','ctrlBam','peaksDir','genomeFastaPath','chrSizes'); # <-- these vars MUST be defined in the %ENV hash

# Note from Alex: ctrlBam can optionally not exist, and that's OK

my $force       = bananas_agw::envLooksTrue("force");
my $verbose     = bananas_agw::envLooksTrue("verbose");
my $debug       = bananas_agw::envLooksTrue("debug");

my $k      = $ENV{sampleName};
my $genome = $ENV{genome};
my $gem    = $ENV{gemJar};   (-e $gem) or confess "[07a.peaksGem.pl]: The required Java jar file for 'gem' (probably called 'gem.jar') does not exist at the following expected path: $gem";


my $finalOutDir = $ENV{peaksDir};
my $outPrefix  = "gem_${k}_${genome}"; # <-- this should be suitable to be a filename, so no spaces/quotes/dollar signs, etc!



my $gemOutParentDir  = catfile($finalOutDir, "${outPrefix}");
my $outPeakBed       = catfile($finalOutDir, "${k}_${genome}_gemPeaks.bed"); # <-- note that this lives OUTSIDE the gem directory, since the gem directory is DELETED at the bottom here!
my $gemOutSubdir     = catfile(              $gemOutParentDir,              "${outPrefix}_outputs");
my $outGpsEvents     = catfile(                                               $gemOutSubdir, "${outPrefix}_1_GEM_events.bed");  # xxx_1_GEM_events.txt: GPS Events predicted using xxx_1_Read_distribution. <--- this is the "less fancy" output
# Note: sometimes the "outGemEvents" file doesn't get generated, if there are "insignificant" results in the first part, with names like xxxx_1_GEM_insignificant.txt
my $outGemEvents     = catfile(                                               $gemOutSubdir, "${outPrefix}_2_GEM_events.bed");  # xxx_2_GEM_events.txt: GEM Events predicted using xxx_2_Read_distribution and xxx_1_KSM/PFMs motifs.



bananas_agw::mkdirOrDie($finalOutDir); # make sure write directories exist, or make them (this is the 'gem'-specific peaks dir)

chdir($finalOutDir); ($verbose) && print STDERR "[07a.peaksGem.pl]: Changing directory to the peaks directory, $finalOutDir...\n";



if (!$force and -e $outPeakBed) {
	($verbose) && print STDOUT "[OK] [SKIPPING RE-RUN]: the output GEM file '$outPeakBed' already exists, so we will not re-peak-call it.\n";
	exit(0); # This is OK, exit without problems
}

my $qVal = 2;	      # sets the q-value threshold (-log10, so 2 = 0.01)
my $kMin = 6;	      # minimum kmer length
my $kMax = 13;	      # maximum kmer length
my $threads = 1;      # Gem uses a TON of threads!!!! Even 1 means "1 times a bunch... number of chromosomes, maybe?" -- ALexW, March 11, 2016. max number of threads  # Note: it's possible to set $threads based on the user's "privileged"-ness status from bananas.pm. (i.e. admin/paid users could have 4 threads, regular users just 1 or something.)
my $gb = 25; # 10 GB crashes sometimes. AlexW increased it to 25 on March 11, 2016

# probably should be pbs_mem (without 'gb') and pbs_ncpus

# See GEM documentation for details about the output format: http://groups.csail.mit.edu/cgs/gem/#QA
#  xxx_0_Read_distribution.txt: The input read distribution (specified by --d).
#  xxx_0_GEM_events.txt: GPS Events used to re-estimate the read distribution of current dataset.
#  xxx_1_Read_distribution.txt: The read distribution estimated from xxx_0_GEM_events.
#  xxx_1_GEM_events.txt: GPS Events predicted using xxx_1_Read_distribution.
#  xxx_1_KSM/PFMs.txt: Motifs discovered using GPS events.
#  xxx_2_Read_distribution.txt: The read distribution estimated from GPS events.
#  xxx_2_GEM_events.txt: GEM Events predicted using xxx_2_Read_distribution and xxx_1_KSM/PFMs motifs.
#  xxx_2_KSM/PFMs.txt: Motifs discovered using GEM events.

my $expBam  = $ENV{expBam};   (-e $expBam) or confess "[07a.peaksGem.pl]: The required EXPERIMENTAL GROUP bam file appears not to exist at the following expected path: $expBam";
my $ctrlBam = $ENV{ctrlBam};  (bananas_agw::isNA($ctrlBam) || (-e $ctrlBam)) or confess "[07a.peaksGem.pl]: The expected CONTROL bam file (which was NOT 'NA', so it was specified as something, at least) appears not to exist at the following expected path: $ctrlBam";

my $javaGemOpt = "java -Xmx${gb}G  -jar $gem  --t $threads "
  . " --q $qVal  --d $ENV{'readDist'} --g $ENV{'chrSizes'}  --genome $ENV{'genomeFastaPath'} --expt $expBam --f SAM  --sl "
  . " --outBED --out $outPrefix ";

if ($ENV{seqType} =~ m/^(exo)$/i) { # if exo, use extra options and don't use ctrl sample
	bananas_agw::systemAndLog("$javaGemOpt --k_min $kMin --k_max $kMax --smooth 3 --mrc 20 --nrf", $verbose, $syslog);
} elsif ($ENV{seqType} =~ m/^(chip)$/i) {
	my $ctrlText =  bananas_agw::isNA($ctrlBam)  ?  " "  :  " --ctrl $ctrlBam "; # if it's NA, then don't add a control. Otherwise, add a control!
	bananas_agw::systemAndLog("$javaGemOpt  --k_min $kMin --k_max $kMax $ctrlText ", $verbose, $syslog);
} else {
	my $msg = "[Warning -- SERIOUS ERROR in 07a.peakGems.pl]! Your sequence type is \"$ENV{seqType}\", which is weird because we should really only be running this code if the sequence type is 'chip' or 'exo' (the next command will fail with a 'confess' message)...\n";
	print STDERR $msg; print STDOUT $msg;
}

(-e $outGpsEvents) or confess "07a.peaksGem.pl attempted to generate a GEM output file named \"$outGpsEvents\", but that file does not exist, so PROBABLY the call to gem failed. Possibly the sequence type was wrong (we expect it to be 'chip' or 'exo') or something went wrong with the system call.";
(-e $outGemEvents) or    warn "07a.peaksGem.pl attempted to generate a GEM output file named \"$outGemEvents\", but that file does not exist, so PROBABLY gem just found the results of part 1 to be 'insignificant' and did not generate the _2_ files.";

if (-e $outGpsEvents and !-e $outGemEvents) {
	bananas_agw::systemAndLog("cat $outGpsEvents | tail -n +2 > $outPeakBed", $verbose, $syslog); # Note: removes the *header* line from the final BED file, since it breaks all downstream tools
	bananas_agw::systemAndLog("touch $outPeakBed.WARNING_using_GEM_GPS_peaks_without_motifs.info.txt", $verbose, $syslog); # Note: remove the header line from the final BED file, since it breaks all downstream tools
	warn "07a.peaksGem.pl: Warning: using the perhaps-suboptimal peaks file '$outGpsEvents' since the '_2_' one was not available.";
} else {
	bananas_agw::systemAndLog("cat $outGemEvents | tail -n +2 > $outPeakBed", $verbose, $syslog); # Note: remove the header line from the final BED file, since it breaks all downstream tools
}

bananas_agw::dieIfFileAccessFails($outPeakBed, "[FAILURE in 07a.peaksGem.pl]: We were not able to generate the final output peaks file.");

# Finally, remove all the large files in the GEM directory if everything worked, apparently
if (-d $gemOutParentDir && !$debug) { # Unless we are in debug mode, remove the now-unnecessary gem directory
	($verbose) && print STDERR "Now removing the (unnecessary) 'gem' directory (\"$gemOutParentDir\")...\n";
	#File::Path::remove_tree($gemOutParentDir); # <-- remove EVERYTHING in this output prefix (gem directory)
	bananas_agw::systemAndLog("/bin/rm $gemOutSubdir/gem*_GEM_event* ", $verbose, $syslog);
	bananas_agw::systemAndLog("/bin/rm $gemOutSubdir/gem*_GEM_insignif* ", $verbose, $syslog);
	bananas_agw::systemAndLog("/bin/rm $gemOutSubdir/gem*_Read_dist* ", $verbose, $syslog);
	bananas_agw::systemAndLog("touch $gemOutSubdir/000_Note__we_automatically_deleted_the_large_files_here--you_can_prevent_this_with_monkeys_debug_option.touch ", $verbose, $syslog);
}

# If you run java -jar gem.jar with no arguments, you get these options:
# GEM command line options (see more options at our website)
#   Required parameters:
#   --d <read spatial distribution file>
#   --exptX <aligned read file for expt (X is condition name)>
# Required GEM motif discovery parameters, optional for GPS-only analysis:
#   --k <length of the k-mer for motif finding, use --k or (--kmin & --kmax)>
#   --k_min <min value of k, e.g. 6>
#   --k_max <max value of k, e.g. 13>
#   --seed <exact k-mer string to jump start k-mer set motif discovery>
#   --genome <the path to the genome sequence directory, for motif finding>
# Optional parameters:
#   --ctrlX <aligned reads file for ctrl (for each condition, ctrlX should match exptX)>
#   --g <genome chrom.sizes file with chr name/length pairs>
#   --f <read file format, BED/SAM/BOWTIE/ELAND/NOVO (default BED)>
#   --s <size of mappable genome in bp (default is estimated from genome chrom sizes)>
#   --a <minimum alpha value for sparse prior (default is esitmated from the whole dataset coverage)>
#   --q <significance level for q-value, specify as -log10(q-value) (default=2, q-value=0.01)>
#   --t <maximum number of threads to run GEM in paralell (default=#CPU)>
#   --out <output folder name and file name prefix>
#   --k_seqs <number of binding events to use for motif discovery (default=5000)>
# Optional flags:
# --fa use a fixed user-specified alpha value for all the regions
#         --help print this help information and exit
