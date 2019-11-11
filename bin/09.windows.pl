#!/usr/bin/perl

my $syslog = "X.09.windows.syslog.txt";

use strict; use warnings; use Carp; # Carp = "confess" function
use File::Spec::Functions qw(catfile);   use Data::Dumper;   use File::Basename;
use lib "/wynton/group/gladstone/biocore/monkey"; use bananas_agw; # To make sure we find bananas_agw

bananas_agw::requireEnvOrDie('force','genome','bedmap','sortBed','densityDir','windowDir','sampleName','geneWindows','analyze'); # <-- these vars MUST be defined in the %ENV hash
my ($force, $verbose, $debug) = bananas_agw::envLooksTrue("force", "verbose", "debug");
my $aDir   = $ENV{'windowDir'};
my $k      = $ENV{'sampleName'};
my $genome = $ENV{'genome'};

bananas_agw::dieIfFileAccessFails($ENV{'analyze'}); # analyze is probably something like "09.xAnalyze.RScript"

my $bf         = $aDir . "/" . $k . "_" . basename($ENV{'geneWindows'}) . ".bins";
my $wf         = $aDir . "/" . $k . "_geneWindows.bed";
my $mf3        = $aDir . "/" . $k . ".map";
my $mf4        = $aDir . "/" . $k . "_geneWindows.txt";
my $sumPdf     = $aDir . "/" . $k . "_geneWindows_sum.pdf";
my $logfile    = $aDir . "/" . $k . ".windows.log";






sub ourWriteToFile($$) {
	my ($file, $textToWrite) = @_;
	($file =~ '^>') or confess "$file must have '>' or '>>' prefixed\n";
	open(OUTFILE,"$file") or confess "error writing to file: $file\n";
	print OUTFILE $textToWrite;
	close(OUTFILE);
}

sub step0Filter {
	my ($infile,$pcFlank,$resolution,$ofname,$sortBed) = @_;
	my $numbins = int((100 + 2*$pcFlank)/$resolution);
	($resolution <= 20) or confess "Need a resolution lower than 20% of the gene: (your selected resolution, which was \"$resolution\") is too high--it cannot be above 20.\n";
	open(OF,">$infile.tmp") or confess "Cannot open outfile temp.out for writing\n";
	open(F1,$infile) or confess "Cannot open $infile\n";
	while (<F1>) {
		chomp;
		if ($_ =~ /^chr/) {
			my @temp = split(/\t/,$_);
			my $chr = $temp[0]; my $start = $temp[1]; my $stop = $temp[2]; my $name=$temp[3];
			my $length = $stop - $start;
			my $binSize = int($resolution/100*$length);
			if ($binSize < 20) {
				$binSize=20;
			}

			my $strand = $temp[4];
			my $tss = ($strand eq "+") ?   $start  : $stop;
			if ($strand eq "+") {
				my $thisBin;
				my $endBin = $tss - int($pcFlank*$length/100);
				unless ($endBin<100) {
					printf OF "$_\n";
				}
			} elsif ($strand eq "-") {
				my $thisBin;
				my $endBin = $tss + int($pcFlank*$length/100);
				unless ($start - int($pcFlank*$length/100) - 1 < 100) {
					printf OF "$_\n";
				}
			} else {
				confess qq{no strand selected--we expected a "+" or a "-", but instead we got this: \"$strand\" \n};
			}
		}
	}
	close(F1) or confess "failed to close file";
	bananas_agw::systemAndLog("$sortBed $infile.tmp > $ofname"); bananas_agw::dieIfFileAccessFails($ofname);
	unlink("$infile.tmp");
}

sub step1Bins {
	my ($infile,$pcFlank,$resolution,$ofname,$sortBed) = @_;
	my $numbins = int((100 + 2*$pcFlank)/$resolution);
	($resolution <= 20) or confess "Need a resolution lower than 20% of the gene: (your selected resolution, which was \"$resolution\") is too high--it cannot be above 20.\n";
	my $onbin = 0;
	my $online = 0;
	open(OF,">$infile.tmp") or confess "Cannot open outfile temp.out for writing\n";
	open(F1,$infile) or confess "Cannot open $infile\n";
	while (<F1>) {
		chomp;
		if ($_ =~ /^chr/) { # looks for this literal text ... ???
			my @temp = split(/\t/,$_);
			my $chr = $temp[0]; my $start = $temp[1]; my $stop = $temp[2]; my $name=$temp[3];
			my $length = $stop - $start;
			my $binSize = int($resolution/100*$length);
			if ($binSize==0) {
				$binSize=1;
			}

			my $strand = $temp[4];
			my $tss = ($strand eq "+") ?   $start  : $stop;
			if ($strand eq "+") {
				my $thisBin;
				my $endBin = $tss - int($pcFlank*$length/100);
				if ($endBin<1) {
					confess "This line is too close to the edge of the chromosome:\n$_\n";
				}
				for (my $i=0; $i<$numbins; $i++) {
					$thisBin = $endBin;
					$endBin = $thisBin + $binSize;
					if ($thisBin<1) {
						confess "bin is negative (+strand): $_\ni = $i, endBin = $endBin\n";
					}
					printf OF "$chr\t$thisBin\t$endBin\tid\t$onbin\n";
					$onbin++;
				}
			} elsif ($strand eq "-") {
				my $thisBin;
				my $endBin = $tss + int($pcFlank*$length/100);
				if ($start - int($pcFlank*$length/100) - 1 < 1) {
					confess "Line is too close to chromosome edge:\n$_\n";
				}
				for (my $i=0; $i<$numbins; $i++) {
					$thisBin = $endBin;
					$endBin = $thisBin - $binSize;
					if ($endBin<1) {
						confess "bin is negative: (-strand): $_\ni = $i, endBin = $endBin\n";
					}
					printf OF "$chr\t$endBin\t$thisBin\tid\t$onbin\n";
					$onbin++;
				}
			} else {
				confess "no strand selected\n";
			}
			$online++;
		}
	}
	close(F1) or die;
	bananas_agw::systemAndLog("$sortBed $infile.tmp \> $ofname");
	unlink("$infile.tmp");
}

# generates the mf3 file
sub step2Map {
	my ($k,$theGenome,$aDir,$dDir) = @_;
	my $input_tf = bananas_agw::catRequiredFile($dDir, "${k}_${theGenome}_tagDensity.bed"); # I guess this is from the density directory
	my $mf1_tmp  = $aDir . "/" . $k . ".map.tmp";
	my $mf2_tmp  = $aDir . "/" . $k . ".paste.tmp";
	bananas_agw::systemAndLog($ENV{'bedmap'} . " --max $bf $input_tf > $mf1_tmp", $verbose); bananas_agw::dieIfFileAccessFails($mf1_tmp);
	bananas_agw::systemAndLog("paste $bf $mf1_tmp    > $mf2_tmp", $verbose);                 bananas_agw::dieIfFileAccessFails($mf2_tmp);
	bananas_agw::systemAndLog("sort -n -k 5 $mf2_tmp > $mf3", $verbose);                     bananas_agw::dieIfFileAccessFails($mf3);
	unlink($mf1_tmp); unlink($mf2_tmp);
}

sub step3Split {
	my ($inbed,$mapfi,$nbins) = @_;
	my $b;
	my $onb=0;
	open(F1,$inbed) or confess "Cannot open $inbed\n";
	while (<F1>) {
		chomp;
		if ($_ =~ /^chr/) {
			my @temp=split(/\t/,$_);
			$b->[$onb]->{chr} = $temp[0];
			$b->[$onb]->{start} = $temp[1];
			$b->[$onb]->{stop} = $temp[2];
			$b->[$onb]->{name} = $temp[3];
			$b->[$onb]->{strand} = $temp[4];
			$onb++;
		}
	}
	close(F1);

	my $onitem=0;
	my $online=1;
	my $tmtx="";
	open (F2,$mapfi) or confess "Cannot open $mapfi\n";
	open (OF,">$mf4") or confess "Cannot open $mf4 for writing\n";
	while (<F2>) {
		chomp;
		if ($_ =~ /^chr/) {
			my @temp = split(/\t/,$_);
			$tmtx = $tmtx . "\t" . $temp[5];
			if ($online % $nbins == 0) {
				printf OF "$b->[$onitem]->{chr}\t$b->[$onitem]->{start}\t$b->[$onitem]->{stop}\t$b->[$onitem]->{name}\t$b->[$onitem]->{strand}"; # looks like it makes a BED line maybe?
				printf OF "$tmtx\n";
				$onitem++;
				$tmtx = "";
			}
			$online++;
		}
	}
	close(F2);
	close(OF);
}

# run through steps...
if (!$force and -e $sumPdf and -e $mf4) {
	($verbose) and print "[OK] Skipping re-generation of the existing files <$sumPdf> and <$mf4>.\n";
	exit; # we have ALREADY finished and generated the required output files
}

bananas_agw::mkdirOrDie($aDir); # Make the output directory

step0Filter($ENV{'geneWindows'}, 20, 1, $wf, $ENV{'sortBed'});
step1Bins($wf, 20, 1, $bf, $ENV{'sortBed'});
step2Map($k, $genome, $aDir, $ENV{'densityDir'});
step3Split($wf, $mf3, 140, $mf4); # generates 'mf4' file
bananas_agw::dieIfFileAccessFails($mf4);

# post-process and analyze
bananas_agw::systemAndLog("sed -i 's/NAN/0/g' $mf4"); # changes NAN to 0 in mf4. Note that it changes it IN PLACE
bananas_agw::systemAndLog("$ENV{'analyze'} $mf4 $sumPdf $aDir 1> $logfile 2> $logfile"); # the second argument to 'analyze' is for a PDF file that we will generate
bananas_agw::dieIfFileAccessFails($sumPdf);

# I guess we do not need any of these (possibly large) files anymore...
unlink($logfile); # Weird that we generate this file only to immediately delete it...
unlink($bf);
unlink($wf);
unlink($mf3); # unclear why we remove this file here...
