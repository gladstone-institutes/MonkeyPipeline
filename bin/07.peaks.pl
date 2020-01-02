#!/bin/perl

my $syslog = "X.07.peaks.syslog.txt";

use strict; use warnings; use Carp; # Carp = "confess" function
use File::Spec::Functions qw(catfile);
use Data::Dumper;
use lib "/wynton/group/gladstone/biocore/MonkeyPipeline"; # To make sure we find bananas_agw
use bananas_agw;

my $tagDensityThreshold = 100;

# check for correct qsub variables needed to perform density calculations
bananas_agw::requireEnvOrDie('genome','seqType','bedmap','bedops','sortBed','tagsDir','peaksDir','densityDir','sampleName');
my ($force, $verbose, $debug) = bananas_agw::envLooksTrue("force", "verbose", "debug");

bananas_agw::dieIfFileAccessFails($ENV{'bedops'}, "A required executable is missing.");
bananas_agw::dieIfFileAccessFails($ENV{'bedmap'}, "A required executable is missing.");
bananas_agw::dieIfFileAccessFails($ENV{'sortBed'}, "A required executable is missing.");

my $k          = $ENV{'sampleName'};
my $genome     = $ENV{'genome'};
bananas_agw::mkdirOrDie($ENV{'peaksDir'}); # make sure write directories exist, or make them

sub getBoundsPerPeak {
    my ($inFile,$outFile,$polarity) = @_;
    my $onp      = "";
    my $saved    = "";
    my $lastsave = "";
    my $maxv     = 0;
    my $fr       = 1;
    open(F0,$inFile) or confess "Cannot open tag count file: $inFile";
    open(OF,">$outFile");
    while(<F0>) {
	chomp;
	my @t = split(/\t/,$_);
	if ($fr) {
	    $fr=0;
	    $onp = $t[3];
	    $maxv = $t[4];
	    $saved = $_;
	}
	elsif ($onp ne $t[3]) {
	    if ($maxv>0) {
		printf OF "$saved\n";
	    }
	    else {
	    }
	    $onp = $t[3];
	    $maxv = $t[4];
	    $saved = $_;
	}
	elsif ($maxv < $t[4]) {
	    $maxv = $t[4];
	    $saved = $_;
	}
	elsif ($maxv == $t[4] && $polarity eq "-") {
            $saved = $_;
	}
    }
    close(F0);
    printf OF "$saved\n";
    close(OF);
}

sub getFootprints {
    my ($plus,$minus,$output) = @_;
    my $ms;
    my $ps;
    my $x;

    open(F0,$plus) or confess "cannot open plus strand tags file: '$plus'";
    while(<F0>) {
	chomp;
	my @t = split(/\t/,$_);
	$x->{$t[3]} = 1;
	$ps->{$t[3]}->{chr} = $t[0];
	$ps->{$t[3]}->{pos} = $t[1];
	$ps->{$t[3]}->{sc}  = $t[4];
    }
    close(F0);

    open(F0,$minus) or confess "cannot open minus strand tags file '$minus'";
    while(<F0>) {
	chomp;
	my @t = split(/\t/,$_);
	$x->{$t[3]} = 1;
	$ms->{$t[3]}->{chr} = $t[0];
	$ms->{$t[3]}->{pos} = $t[2];
	$ms->{$t[3]}->{sc}  = $t[4];
    }
    close(F0);

    my $fp = 1;
    open(OF,">$output") or confess "cannot open $output for writing";
    foreach my $k (keys %$x) {
	if (defined($ms->{$k})) {
	    if (defined($ps->{$k})) {
		my $score = $ps->{$k}->{sc} + $ms->{$k}->{sc};
		if ($ps->{$k}->{pos} < $ms->{$k}->{pos}) {
		    printf OF "$ms->{$k}->{chr}\t$ps->{$k}->{pos}\t$ms->{$k}->{pos}\tFP\_$fp\t$score\n";
		    $fp++;
		}
		elsif ($ps->{$k}->{pos} > $ms->{$k}->{pos}) {
		    printf OF "$ms->{$k}->{chr}\t$ms->{$k}->{pos}\t$ps->{$k}->{pos}\tFP\_$fp\t$score\n";
                    $fp++;
		}
		else {
		    my $sp = $ms->{$k}->{pos} - 1;
		    my $ep = $ms->{$k}->{pos} + 1;
		    printf OF "$ms->{$k}->{chr}\t$sp\t$ep\tFP\_$fp\t$score\n";
		}
		
	    }
	}
    }
    close(OF);
}

sub runExoAnalysis {
    my ($k,$genome,$peaks,$tagsDir,$peaksDir,$sortBed,$bedmap) = @_;
    
    # get + and - strand tags from tags folder
    my $pstf  = catfile($tagsDir, "${k}_${genome}_tags.pos.bed");
    my $ngtf  = catfile($tagsDir, "${k}_${genome}_tags.neg.bed");    
    (-e $pstf) or confess "pos strand tag file not found: $pstf";
    (-e $ngtf) or confess "pos strand tag file not found: $ngtf";

    # get 1bp bins for peak regions
    my $pk1bpf = catfile($peaksDir, "${k}_${genome}_peaks.1bp.bed");
    open (PF,$peaks) or confess "Cannot open peak file $peaks";
    open (OP,">$pk1bpf") or confess "Cannot open 1bp peak file for writing";
    while(<PF>) {
	chomp;
	my @t = split(/\t/,$_);
	my $stop = $t[1];
	for (my $i=0; $i < $t[2] - $t[1] - 1; $i++) {
	    my $start = $stop;
	    $stop++;
	    printf OP "$t[0]\t$start\t$stop\t$t[3]\n";
	}
    }
    close (OP);
    close (PF);

    # sort 1bp bins file
    my $tfn = catfile($peaksDir, "${k}_${genome}_peaks.1bp.bed.sorted");
    bananas_agw::systemAndLog("$sortBed $pk1bpf > $tfn", $verbose, $syslog);
    bananas_agw::systemAndLog("mv -f $tfn $pk1bpf", $verbose, $syslog);

    # do tag counts across peak regions
    my $posTagCountsFile = catfile($peaksDir, "${k}_${genome}_posTagDensity10bp_inPeaks_1bpRes.bed");
    my $negTagCountsFile = catfile($peaksDir, "${k}_${genome}_negTagDensity10bp_inPeaks_1bpRes.bed");
    my $tptcf = catfile($peaksDir, "${k}_${genome}_posTag.counts");
    my $tntcf = catfile($peaksDir, "${k}_${genome}_negTag.counts");
    bananas_agw::systemAndLog("$bedmap --count --range 20 $pk1bpf $pstf > $tptcf", $verbose, $syslog);
    bananas_agw::systemAndLog("paste $pk1bpf $tptcf > $posTagCountsFile", $verbose, $syslog);
    unlink($tptcf);
    bananas_agw::systemAndLog("$bedmap --count --range 20 $pk1bpf $ngtf > $tntcf", $verbose, $syslog);
    bananas_agw::systemAndLog("paste $pk1bpf $tntcf > $negTagCountsFile", $verbose, $syslog);
    unlink($tntcf);

    # get high-res strand-specific peaks
    getBoundsPerPeak($posTagCountsFile,$tptcf,"+");
    getBoundsPerPeak($negTagCountsFile,$tntcf,"-");

    # call footprints
    my $fpCalls = catfile($peaksDir, "${k}_${genome}_footprints.bed");
    getFootprints($tptcf,$tntcf,$fpCalls);
    bananas_agw::systemAndLog("$sortBed $fpCalls > $tptcf", $verbose, $syslog);
    bananas_agw::systemAndLog("mv -f $tptcf $fpCalls", $verbose, $syslog);
    unlink($tptcf);
    unlink($tntcf);
}

my $fn  = catfile($ENV{'densityDir'}, "${k}_${genome}_tagDensity.bed");

my $tmpA = catfile($ENV{'peaksDir'}, "${k}_${genome}.TEMP_A");
my $tmpB = catfile($ENV{'peaksDir'}, "${k}_${genome}.TEMP_B");
my $of1 = catfile($ENV{'peaksDir'}, "${k}_${genome}_tagDensity.bed.filter");
my $of2 = catfile($ENV{'peaksDir'}, "${k}_${genome}_tagDensity.bed.filter.merged");
my $of3 = catfile($ENV{'peaksDir'}, "${k}_${genome}_peaks.bed");

## ~~~~~~~~~~~~~~~~ The actual program execution occurs HERE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# check prior completion under various conditions. However, if "$force" is specified, then we ALWAYS go ahead with the code below!

my $fpCallsFilename = catfile($ENV{'peaksDir'}, "${k}_${genome}_footprints.bed"); # check for existence of final output file for exo-seq ONLY

if (!$force and $ENV{'seqType'} eq "exo" and -e $of3 and -e $fpCallsFilename) { # Exo-seq ONLY
	($verbose) && print STDOUT "[OK] [SKIPPING RE-RUN]: the exo-seq-specific output already exists at '$of3' and '$fpCallsFilename'. No need to re-generate it.";
	exit(0);
}

if (!$force and $ENV{'seqType'} ne "exo" and -e $of3) {	# Anything EXCEPT for exo-seq. If it's NOT 'exo'-seq, then we just need $of3 to exist
	($verbose) && print STDOUT "[OK] [SKIPPING RE-RUN]: the non-exo-seq output already exists at '$of3'. No need to re-generate it.";
	exit(0);
}



bananas_agw::dieIfFileAccessFails($fn, "[FAILURE in 07.peaks.pl]: We failed to find the 'fn' input.");

bananas_agw::systemAndLog("awk '{if (\$5 > $tagDensityThreshold) {print \$0}}' $fn > $tmpA", $verbose, $syslog); # generating temp file A
bananas_agw::systemAndLog("$ENV{'bedops'} -m --range 150 $tmpA > $tmpB", $verbose, $syslog); # generating temp file B
bananas_agw::systemAndLog("awk -v x=0 '{p1=\$2+150;p2=\$3-150;y = (p2 - p1); if (y>=20) {x=x+1; print \$1\"\\t\"p1\"\\t\"p2\"\\t$k\_pk\_\"x}}' $tmpB > $of1", $verbose, $syslog); # generating out file #1
bananas_agw::systemAndLog("$ENV{'bedmap'} --max $of1 $fn > $of2", $verbose, $syslog); # generating out file #2
bananas_agw::systemAndLog("paste $of1 $of2 | awk '{x=int(\$5); print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"x}'> $of3", $verbose, $syslog); # Generating out file #3
unlink($tmpA); # Remove temp file
unlink($tmpB); # Remove temp file
unlink($of1); # Remove temp file
unlink($of2); # Remove temp file

if ($ENV{'seqType'} eq "exo") {
	# Run some additional analysis for exo-seq
	runExoAnalysis($k, $genome, $of3, $ENV{'tagsDir'}, $ENV{'peaksDir'}, $ENV{'sortBed'}, $ENV{'bedmap'});
}
