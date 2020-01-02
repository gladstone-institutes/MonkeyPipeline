#!/bin/perl -w
# This runs after BCP had done its thing, and summarizes a bunch of files. (See "Percent tags in peaks" below)
use strict; use warnings; use Carp; # Carp = "confess" function
use File::Spec::Functions qw(catfile);   use Data::Dumper;   use File::Basename;
use lib "/wynton/group/gladstone/biocore/MonkeyPipeline"; use bananas_agw; # To make sure we find bananas_agw
bananas_agw::requireEnvOrDie('genome','bedmap','tagsDir','peaksDir','sampleName'); # <-- these vars MUST be defined in the %ENV ha
my ($force, $verbose) = bananas_agw::envLooksTrue("force", "verbose");
-e $ENV{bedmap}   or confess "Failed to find the required 'bedmap' executable,  which was supposed to be located here: $ENV{bedmap}";
-d $ENV{peaksDir} or confess "Failed to find the required 'peaksDir' directory, which was supposed to be located here: $ENV{peaksDir}";
-d $ENV{tagsDir}  or confess "Failed to find the required 'tagsDir' directory,  which was supposed to be located here: $ENV{tagsDir}";

sub roundToTwoDecimalPlaces {  my ($num) = @_; return(sprintf("%.2f", $num));  }

sub getSumOfNumbersInFile {
    my ($filename) = @_;
    my $sum = 0;
    open(SUMFILE,$filename) or die "cannot open $filename";
    while(<SUMFILE>) { chomp; $sum += $_; }
    close(SUMFILE);
    return($sum);
}

sub getPercentTagsInPeaks {
    my ($ofTemp,$peaks,$tags) = @_;
    system("$ENV{bedmap} --count ${peaks} ${tags} > $ofTemp");
    my $totalCount  = getSumOfNumbersInFile($ofTemp);
    my $nTags       = bananas_agw::numLinesInFile($tags);
    my $percent     = (0 == $nTags) ? "DIVIDE_BY_ZERO" : roundToTwoDecimalPlaces(100*$totalCount/$nTags);
    unlink($ofTemp); # remove temp file
    return($percent);
}

sub ourWriteToFile {
    my ($file, $text) = @_;
    ($file =~ '^>') or confess "the 'file' needs to start with either a single caret ('>' = write to new file, erasing anything that was there) or '>>' (append to existing file). You passed in this: $file";
    open(OUTFILE,"$file") or confess "Failed to write to file: $file"; print OUTFILE $text; close(OUTFILE);
}

my $ptipFile = catfile($ENV{'peaksDir'}, "0.percentTagsInPeaks.txt");
if ($force or (! -e $ptipFile)) {
    ourWriteToFile(">${ptipFile}", join("\t", ("Sample_name", "Percent Tags In Peaks"))."\n");
    # open tags directory and build list of samples and inputs
    my $samples; # hash pointer
    opendir(INDIR,$ENV{'peaksDir'}) or confess "Cannot open the following 'peaksDir' directory: $ENV{'peaksDir'}";
    while (my $file = readdir(INDIR)) {
	if ($file =~ /^(.*?\.\d{1,})\_.*?\_peaks.bed$/) {
	    $samples->{$1} = $file;
	}
    }
    closedir(INDIR);
    foreach my $k (sort keys %$samples) {
	my $ofTemp    = catfile($ENV{peaksDir}, "${k}.ptip");
	my $ofSorted  = catfile($ENV{peaksDir},$samples->{$k});
	my $tags      = catfile($ENV{tagsDir} , "${k}_$ENV{genome}_tags.bed");
	my $pc        = getPercentTagsInPeaks($ofTemp, $ofSorted, $tags);
	ourWriteToFile(">>${ptipFile}", "${k}\t${pc}\n");
    }
}
