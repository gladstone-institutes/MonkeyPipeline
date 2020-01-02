#!/bin/perl

use strict; use warnings; use Carp; # Carp = "confess" function
use File::Spec::Functions qw(catfile);   use Data::Dumper;   use File::Basename;
use lib "/wynton/group/gladstone/biocore/MonkeyPipeline"; use bananas_agw; # To make sure we find bananas_agw

my $syslog = "X.06.density.syslog.txt";

#printf STDERR Dumper \%ENV; die;
bananas_agw::requireEnvOrDie('force','genome','seqType','sampleName','bedmap','binsFile','binsFileID','tagsDir','densityDir'); # <-- these vars MUST be defined in the %ENV hash
my $sampleFile = catfile($ENV{'tagsDir'},$ENV{'sampleName'}) . "_" . $ENV{'genome'} . "_tags.bed";
my $inputFile = "NA";
my $isMatched = 0;
my ($force, $verbose) = bananas_agw::envLooksTrue("force", "verbose");
if (defined($ENV{'inputName'}) and !bananas_agw::isNA($ENV{'inputName'})) {
	$inputFile = catfile($ENV{'tagsDir'},$ENV{'inputName'}) . "_" . $ENV{'genome'} . "_tags.bed";
	$isMatched = 1;
}

bananas_agw::mkdirOrDie($ENV{'densityDir'});
my $countsDir = catfile($ENV{'densityDir'},"counts"); # someDensityDir/counts/
bananas_agw::mkdirOrDie($countsDir);

my $HARD_CODED_TAG_SMOOTHING_DEFAULT = 75; # ??
my $HARD_CODED_TAG_SMOOTHING_FOR_RNA = 1;  # ??

sub normalizeCounts {
	my ($itags,$iCounts,$tags,$counts) = @_;
	my $c = bananas_agw::numLinesInFile($tags);
	my $tfile = "${counts}.normTmp";
	if (uc($itags) ne "NA" && uc($iCounts) ne "NA") {
		# matched input
		my $ic = bananas_agw::numLinesInFile($itags);
		my $ib = 0;
		my $itc;
		open(TF,$iCounts) or die "Cannot open counts file: $counts\n";
		while (<TF>) {
			$itc->[$ib] = $_;
			$ib++;
		}
		close(TF);

		my $i=0;
		open(OF,">$tfile") or die "Cannot open temporary file for normalization: $tfile\n";
		open(TF,$counts) or die "Cannot open count file for tag density normalization: $counts\n";
		while (<TF>) {
			my $x = ($_/$c) - ($itc->[$i]/$ic);
			my $v = int($x*$ib*10+0.49999)/10; # <-- 
			if ($v<0) {
				$v=0;
			}
			printf OF "$v\n";
			$i++;
		}
		close(TF);
		close(OF);
		system("mv -f $tfile $counts");
	} else {
		# no matched input
		my $b = bananas_agw::numLinesInFile($counts);
		open(OF,">$tfile") or die "Cannot open temporary file for normalization: $tfile\n";
		open(TF,$counts) or die "Cannot open count file for tag density normalization: $counts\n";
		while (<TF>) {
			my $v = int(  ($_ * $b / $c)   *100)/100;
			printf OF "$v\n"; 
		}
		close(TF);
		close(OF);

		bananas_agw::systemAndLog("mv -f $tfile $counts", $verbose, $syslog);
	}
}

my $tagSmooth = undef;   # sequence-type-dependent smoothing (RNA gets one value, everything else gets another)
if (lc($ENV{'seqType'}) eq "rna") {
	$tagSmooth = $HARD_CODED_TAG_SMOOTHING_FOR_RNA;
} else {
	$tagSmooth = $HARD_CODED_TAG_SMOOTHING_DEFAULT;
}

# get normalized tag density counts for data with matched/unmatched background tags
my $of = catfile($ENV{'densityDir'},$ENV{'sampleName'}) . "_" . $ENV{'genome'} . "_tagDensity.bed";

if ($force or ! -e $of) {
	if ($isMatched) {
		my $cof1 = catfile($ENV{'densityDir'},$ENV{'sampleName'}) . "_" . $ENV{'genome'} . "_tagDensity.bed.count";
		my $cof2 = catfile($ENV{'densityDir'},$ENV{'inputName'}) . "_" . $ENV{'sampleName'} . "_" . $ENV{'genome'} . "_tagDensity.bed.count";    
		system("$ENV{'bedmap'} --count --range $tagSmooth $ENV{'binsFile'} $sampleFile > $cof1");
		system("$ENV{'bedmap'} --count --range $tagSmooth $ENV{'binsFile'} $inputFile > $cof2");
		normalizeCounts($inputFile,$cof2,$sampleFile,$cof1);
		system("paste $ENV{'binsFileID'} $cof1 > $of");
		system("mv $cof1 $countsDir/");
		unlink($cof2);
	} else {
		my $cof1 = catfile($ENV{'densityDir'},$ENV{'sampleName'}) . "_" . $ENV{'genome'} . "_tagDensity.bed.count";
		system("$ENV{'bedmap'} --count --range $tagSmooth $ENV{'binsFile'} $sampleFile > $cof1");
		normalizeCounts("NA","NA",$sampleFile,$cof1);
		system("paste $ENV{'binsFileID'} $cof1 > $of");
		system("mv $cof1 $countsDir/");
	}
}

# commented out for testing
#my $ofNuc = catfile($ENV{'densityDir'},$ENV{'sampleName'}) . "_" . $ENV{'genome'} . "_nucTagDensity.bed";
#if ($ENV{'seqType'} eq "atac") {
#if ($force || ! -e $ofNuc) {  
#my $ofShtCounts = catfile($countsDir,$ENV{'sampleName'}) . "_" . $ENV{'genome'} . "_tags.bed.count";
#my $ofNucCounts = catfile($ENV{'densityDir'},$ENV{'sampleName'}) . "_" . $ENV{'genome'} . "_nucTagDensity.bed.count";
#my $ofNucTmp1 = catfile($ENV{'densityDir'},$ENV{'sampleName'}) . "_" . $ENV{'genome'} . "_nucTagDensity.bed.tmp1";
#my $ofNucTmp2 = catfile($ENV{'densityDir'},$ENV{'sampleName'}) . "_" . $ENV{'genome'} . "_nucTagDensity.bed.tmp2";
#my $longFile = catfile($ENV{'tagsDir'},$ENV{'sampleName'}) . "_" . $ENV{'genome'} . "_longFrags.bed";

#system("$ENV{'bedmap'} --count --range $tagSmooth $ENV{'binsFile'} $longFile  > $ofNucCounts");
#system("paste $ofShtCounts $ofNucCounts > $ofNucTmp1");
#system("cat $ofNucTmp1 | awk '{x = log((\$2+1)/(\$1+1))*abs(\$2-\$1) }' > $ofNucTmp2");
#system("paste $ENV{'binsFileID'} $ofNucTmp2 > $ofNuc");
#}
#}
