#!/usr/bin/perl





my $syslog = "X.07b.bcpPeaks.syslog.txt";

use strict; use warnings; use Carp; # Carp = "confess" function
use File::Spec::Functions qw(catfile);
use Data::Dumper;
use lib "/data/work/Code/alexgw/monkey_agw"; use bananas_agw; # To make sure we find bananas_agw
use Cwd;

# check for correct qsub variables needed to perform density calculations
bananas_agw::requireEnvOrDie('sampleName','genome','seqType','bcp','expTagsBed','ctlTagsBed','peaksDir'); # <-- these vars MUST be defined in the %ENV hash

my $force = bananas_agw::envLooksTrue("force"); my $verbose = bananas_agw::envLooksTrue("verbose"); my $debug = bananas_agw::envLooksTrue("debug");

my $qVal   = 0.001;  # sets the q-value threshold
my $k      = $ENV{'sampleName'};
my $genome = $ENV{'genome'};
my $bcp    = $ENV{'bcp'};

# make sure write directories exists or make them
((-d $ENV{'peaksDir'}) or mkdir($ENV{'peaksDir'})) or confess "unable to find or create results directory: $ENV{'peaksDir'}!";
(-e $ENV{'expTagsBed'}) or confess "expTagsBed (experiment tags bed file, which is required) file does not exist: $ENV{'expTagsBed'}";
(-e $ENV{'ctlTagsBed'}) or confess "ctlTagsBed (control tags bed file, which is required) file does not exist: $ENV{'ctlTagsBed'}";
(-e $bcp) or confess "the required BCP executable does not exist: $bcp";

chdir("$ENV{'peaksDir'}");

my $outPeakFile = catfile($ENV{'peaksDir'},"${k}_${genome}_bcpPeaks.bed");
my $opTmp  = $outPeakFile . ".tmp";
my $opLog  = $outPeakFile . ".log";
my $expTmp = catfile($ENV{'peaksDir'},"${k}_${genome}_bcpExpTmp.bed");
my $ctlTmp = catfile($ENV{'peaksDir'},"${k}_${genome}_bcpCtlTmp.bed");

if ($force or ! -e $outPeakFile) {
	bananas_agw::systemAndLog("awk \'{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\t1\t\"\$5}\' $ENV{'expTagsBed'} > $expTmp", $verbose, $syslog);
	bananas_agw::systemAndLog("awk \'{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\t1\t\"\$5}\' $ENV{'ctlTagsBed'} > $ctlTmp", $verbose, $syslog);
	bananas_agw::systemAndLog("$bcp -1 $expTmp -2 $ctlTmp -f 200 -w 200 -p $qVal -3 $opTmp &> $opLog", $verbose, $syslog);
	bananas_agw::systemAndLog("sort-bed $opTmp > $outPeakFile", $verbose, $syslog);
	bananas_agw::systemAndLog("/bin/rm $opTmp $opLog $ctlTmp $expTmp", $verbose, $syslog);
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
