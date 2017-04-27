#!/usr/bin/perl
use strict; use warnings; use Carp; # Carp = "confess" function
use File::Spec::Functions qw(catfile);   use Data::Dumper;   use File::Basename;
use lib "/data/work/Code/alexgw/monkey_agw"; use bananas_agw; # To make sure we find bananas_agw

bananas_agw::requireEnvOrDie('force','genome','minMapQ','tagsDir','mappingDir'); # <-- these vars MUST be defined in the %ENV hash
my $force  = bananas_agw::envLooksTrue("force");
my $genome                      = $ENV{"genome"};
my $q                           = $ENV{"minMapQ"};

sub roundToTwoDecimalPlaces($) { return(sprintf("%.2f", $_[0])); }

sub getMappingStats {
    my ($name, $tagsFilename, $logFilename) = @_;
    my $total  = 1; 
    my $mapped = "NA";
    if (-e $logFilename) {
	open(LFILE, "< $logFilename") or die "05.xSummary.pl: cannot find the required log file <$logFilename>: maybe the aligner was not bowtie2? (If that is the case, then this file will just plain not exist.) \n";
	while(<LFILE>) {
	    if ($_ =~ /(\d{1,})\sreads\;\sof\sthese\:/i) {
		$total = $1;
	    } elsif ($_ =~ /^(\S{1,})\% overall alignment rate/i) {
		$mapped = int($1/100 * $total);
	    }
	}
	close(LFILE);
    }
    my $finalLineCount = (-e $tagsFilename) ? bananas_agw::numLinesInFile($tagsFilename) : 0;
    my $percent        = roundToTwoDecimalPlaces(100*$finalLineCount/$total);
    if ($percent > 100 ) {
	$finalLineCount /= 2.0; # Alex's note: Maybe... maybe this means it was paired end or something (?)
	$percent        /= 2.0; #              for some reason, we divide everything by two in this case. (?)
    }
    return( join("\t", ($name, $total, $mapped, $finalLineCount, $percent)) );
}

sub ourWriteToFile { # Example usage: ourWriteToFile(">>${msFile}", "${tLine}\n"); (that would be an APPEND operation)
    my ($file, $text) = @_;
    ($file =~ '^>')       or die "unable to write to file: $file";
    open(OUTFILE,"$file") or die "unable to write to file: $file";
    print OUTFILE $text;
    close(OUTFILE);
}

my $msFile = catfile($ENV{'tagsDir'}, "0.tagStatistics.txt");
if ($force or ! -e $msFile) {
	ourWriteToFile(">${msFile}", join("\t", ("Sample_name", "Num_total_tags", "Num_mapped_tags", "Final_mapped_tags", "Percent_kept"))."\n");
	my $samples;
	opendir(INDIR,$ENV{'tagsDir'}) or die "cannot open tags directory: $ENV{'tagsDir'}";
	while (my $file = readdir(INDIR)) {
		if ($file =~ /^(.*?\.\d{1,})\_.*?\_tags.bed$/) {
			$samples->{$1} = $file;
		}
	}
	closedir(INDIR);
	foreach my $k (sort keys %$samples) {
		my $bowtieStatsLogFilename = undef; # Below: check these filenames in order to find the first on that exists. This is because we changed the name previously.
		# note: these are BOWTIE2 log files! So they will not be there if the aligner was something else!
		my @logFilenameGuesses = (catfile($ENV{'mappingDir'},      "${k}_${genome}_q${q}.bam.log")
					  , catfile($ENV{'mappingDir'},    "${k}_${genome}_q${q}.bam.bowtie2.log")
					  , catfile($ENV{'mappingDir'},    "${k}_${genome}_q${q}.bam.bowtie.log")
					  , glob(catfile($ENV{'mappingDir'},"${k}_${genome}*.log")) # Ok, look for ANYTHING that is 'close enough' if we don't find any of the hard-coded guesses above (glob is "wildcard expansion")
					 );
		foreach my $fff (@logFilenameGuesses) {
			if (-e $fff) { $bowtieStatsLogFilename = $fff; last; } # found the log file where it was supposed to be!
		}
		my $ofSorted = catfile($ENV{'tagsDir'},$samples->{$k});
		my $tLine    = getMappingStats($k, $ofSorted, $bowtieStatsLogFilename);
		ourWriteToFile(">>${msFile}", "${tLine}\n");
	}
}
