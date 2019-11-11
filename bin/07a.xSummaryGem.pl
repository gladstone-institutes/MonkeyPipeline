#!/usr/bin/perl





use strict; use warnings; use Carp; # Carp = "confess" function
use File::Spec::Functions qw(catfile);   use Data::Dumper;   use File::Basename;
use lib "/wynton/group/gladstone/biocore/monkey"; use bananas_agw; # To make sure we find bananas_agw

bananas_agw::requireEnvOrDie('force','genome','bedmap','tagsDir','peaksDir','sampleName'); # <-- these vars MUST be defined in the %ENV hash
my $k = $ENV{'sampleName'};
my $genome = $ENV{'genome'};

sub roundToTwoDecimalPlaces {
    my ($num) = @_; return(sprintf("%.2f", $num));
}

sub getSumOfNumbersInFile {
    my ($filename) = @_;
    my $sum = 0;
    open(SUMFILE,$filename) or die "cannot open $filename"; 
    while(<SUMFILE>) { chomp; $sum += $_; }
    close(SUMFILE);
    return($sum);
}

sub getPercentTagsInPeaks {
    my ($bedmap,$ofTemp,$peaks,$tags) = @_;
    system("$bedmap --count ${peaks} ${tags} > $ofTemp");
    my $totalCount  = getSumOfNumbersInFile($ofTemp);
    my $nTags       = bananas_agw::numLinesInFile($tags);         
    my $percent = (0 == $nTags) ? "DIVIDE_BY_ZERO" : roundToTwoDecimalPlaces(100*$totalCount/$nTags);
    unlink($ofTemp);
    return($percent);
}

sub ourWriteToFile {
    my ($file, $text) = @_;
    ($file =~ '^>') or die "unable to write to file: $file";
    open(OUTFILE,"$file") or die "unable to write to file: $file";
    print OUTFILE $text;
    close(OUTFILE);
}

my $ptipFile = catfile($ENV{'peaksDir'}, "0.percentTagsInPeaks.txt");
if ($ENV{'force'} || ! -e $ptipFile) {
    ourWriteToFile(">${ptipFile}", join("\t", ("Sample_name", "Percent Tags In Peaks"))."\n");

    # open tags directory and build list of samples and inputs
    my $samples;
    opendir(INDIR,$ENV{'peaksDir'}) or die "cannot open peaks directory: $ENV{'peaksDir'}";
    while (my $file = readdir(INDIR)) {
	if ($file =~ /^(.*?\.\d{1,})\_.*?\_peaks.bed$/) {
	    $samples->{$1} = $file;
	}
    }
    closedir(INDIR);

    foreach my $k (sort keys %$samples) {
	my $ofSorted  = catfile($ENV{'peaksDir'},$samples->{$k});
	my $ofTemp    = catfile($ENV{'peaksDir'}, "${k}.ptip");
	my $tags      = catfile($ENV{'tagsDir'},"${k}_${genome}_tags.bed");
	my $pc        = getPercentTagsInPeaks($ENV{'bedmap'},$ofTemp, $ofSorted,$tags);
	ourWriteToFile(">>${ptipFile}", "${k}\t${pc}\n");
    }
}
