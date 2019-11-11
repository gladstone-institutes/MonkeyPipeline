#!/usr/bin/perl
my $syslog = "X.05.tags.syslog.txt";

use strict; use warnings; use Carp; # Carp = "confess" function
use File::Spec::Functions qw(catfile);   use Data::Dumper;   use File::Basename;
use lib "/wynton/group/gladstone/biocore/monkey"; use bananas_agw; # To make sure we find bananas_agw

bananas_agw::requireEnvOrDie('genome','seqType','samtools','bam2bed','sortBed','tagsDir','inputFile','sampleName','isPaired'); # <-- these vars MUST be defined in the %ENV hash
my ($force, $verbose) = bananas_agw::envLooksTrue("force", "verbose");

bananas_agw::dieIfFileAccessFails($ENV{'inputFile'}); # input file must exist (and be readable!)

# make sure write directory exists or make it
bananas_agw::mkdirOrDie($ENV{'tagsDir'});
# generate initial bed file

my $pathPre    = catfile($ENV{'tagsDir'},$ENV{'sampleName'}) . "_" . $ENV{'genome'}; # <-- The full path AND ALSO the first part of the file (sample name + genome)
my $oft        = "${pathPre}_tags.bed.tmp";# <-- Universal TEMP file for all methods
my $ofSorted   = "${pathPre}_tags_INTERMEDIATE_COMPUTATION_WILL_BE_DELETED_ON_SUCCESS.tmp.bed"; # <-- Universal TEMP file for all methods
my $ofFinal    = "${pathPre}_tags.bed";     # <-- Universal FINAL output file for all methods
my $pstf       = "${pathPre}_tags.pos.bed"; # EXO-seq only
my $nstf       = "${pathPre}_tags.neg.bed"; # EXO-seq only
my $ofLongTmp  = "${pathPre}_longFrags.bed.tmp"; # ATAC only
my $ofLongTmp2 = "${pathPre}_longFrags.bed.tmp2";# ATAC only
my $ofLongTmp3 = "${pathPre}_longFrags.bed.tmp3";# ATAC only
my $ofLong     = "${pathPre}_longFrags.bed";     # ATAC only

# =============================== CHECK TO SEE WHETHER OR NOT ALL THE EXPECTED OUTPUT FILES ALREADY EXIST (and we aren't forcing a re-run) =================
my @reqFinalFiles = ($ofFinal);
if ($ENV{'seqType'} eq "exo")  { push(@reqFinalFiles, ($pstf, $nstf)); } # EXO-seq requires that two additional files exist--pstf and nstf
if ($ENV{'seqType'} eq "atac") { push(@reqFinalFiles, $ofLong); }
my $allOutFilesExistAlready = 1; # <-- Initially hypothesize they DO all exist, then check for each one to possibly disprove this
foreach my $f (@reqFinalFiles) {
	if (not bananas_agw::fileIsReadableAndNonzero($f)) { $allOutFilesExistAlready = 0; } 	# Note: counts an EMPTY file as "not existing" for the purposes of re-running. This is the behavior we generally want.
}
if (not $force and $allOutFilesExistAlready) {
	($verbose) && print STDOUT "[OK] [SKIPPING RE-RUN]: the required final output files (specifically: " . join(", ", @reqFinalFiles) . ") all already exist.\n";
	exit;
}

# ======== If we get here, then at least one of the required output files does NOT already exist! ==========

# sequence-type-specific processing of tags
if ($ENV{'seqType'} eq "exo") { # process exo tags
    bananas_agw::systemAndLogDieNonzero("$ENV{'bam2bed'} --do-not-sort < $ENV{'inputFile'} | cut -f 1-3,4,6,8 > $oft", $verbose, $syslog);
    bananas_agw::systemAndLogDieNonzero("$ENV{'sortBed'} $oft | awk '{x=\$2;xu=\$2+1; if (\$5==\"-\") {x=\$3-1;xu=\$3}; print \$1\"\t\"x\"\t\"xu\"\t\"\$4\"\t\"\$5}' > $ofSorted", $verbose, $syslog);
    bananas_agw::dieIfFileAccessFails($ofSorted);
    bananas_agw::systemAndLogDieNonzero("awk '{if (\$5==\"+\") {print \$0}}' $ofSorted > $pstf", $verbose, $syslog);
    bananas_agw::systemAndLogDieNonzero("awk '{if (\$5==\"-\") {print \$0}}' $ofSorted > $nstf", $verbose, $syslog);
    bananas_agw::dieIfFileAccessFails($pstf);
    bananas_agw::dieIfFileAccessFails($nstf);
} elsif ($ENV{'seqType'} eq "rna") {
	bananas_agw::systemAndLogDieNonzero("$ENV{'bam2bed'} --do-not-sort < $ENV{'inputFile'} | cut -f 1-3,4,6,8 > $oft", $verbose, $syslog);
	bananas_agw::dieIfFileAccessFails($oft);
	# process rna tags
	my $rnaTmp = $oft . "rnaTmp";
	open(F0,$oft) or confess "[ERROR] 05.tags.pl: cannot open the file '$oft' for filtering spliced reads...";
	open(OF,">$rnaTmp");
	while (<F0>) {
		chomp;
		my @t = split(/\t/,$_);
		my $isGapped = 0;
		while ($t[5] =~ /(\d{1,})N/g) { # if... something has a number and then an N... meaning it's gapped? Field $t[5] might be the CIGAR string maybe?
			if ($1 >= 10) { # $1 is the number before the 'N' in column 5... maybe it's like "42N" and then $1 is 42.
				$isGapped = 1; # I guess we don't print out gapped reads? This is some way to accomplish that maybe. Not sure what $t[5] is specifically, or how it's possible to loop over it if it's just one element
			}
		}
		if (not $isGapped) {
			printf OF "$_\n"; # Looks like we just print the regular line as-is if it's not gapped. Gapped lines... just don't get printed ever maybe? Unclear.
		}
	}
	close(F0);
	close(OF);
	bananas_agw::systemAndLogDieNonzero("$ENV{'sortBed'} $rnaTmp | awk '{x=\$2;xu=\$3; print \$1\"\t\"x\"\t\"xu\"\t\"\$4\"\t\"\$5}' > $ofSorted", $verbose, $syslog);
	unlink($rnaTmp);
} elsif ($ENV{'seqType'} eq "atac") {
	($ENV{'isPaired'}) or confess "[ERROR] 05.tags.pl: ATAC-seq reads MUST be paired end---otherwise we cannot extract nucleosome positioning information.\n";
	# collapse paired reads into fragments
	bananas_agw::systemAndLogDieNonzero("$ENV{'samtools'} view $ENV{'inputFile'} | awk '{x = \$4-1; y = x + \$9; x=x+4; y=y-5; print \$3\"\t\"x\"\t\"y }'                       > $ofLongTmp ", $verbose, $syslog);
	bananas_agw::systemAndLogDieNonzero("cat                   $ofLongTmp        | awk '{if (\$3 > \$2 && \$2 > 0) { print \$0 }}'                                              > $ofLongTmp2", $verbose, $syslog);
	bananas_agw::systemAndLogDieNonzero("$ENV{'sortBed'}       $ofLongTmp2       | uniq                                                                                         > $ofLongTmp3", $verbose, $syslog);
	bananas_agw::systemAndLogDieNonzero("cat                   $ofLongTmp3       | awk '{x=\$2; xd=\$2+1; y=\$3-1; yd=\$3; print \$1\"\t\"x\"\t\"xd; print\$1\"\t\"y\"\t\"yd }' > $oft"       , $verbose, $syslog);
	bananas_agw::systemAndLogDieNonzero("$ENV{'sortBed'}       $oft                                                                                                             > $ofSorted"  , $verbose, $syslog);
	bananas_agw::systemAndLogDieNonzero("cat                   $ofLongTmp3       | awk '{if (\$3-\$2 >=150) {print \$0}}'                                                       > $ofLong"    , $verbose, $syslog);
	bananas_agw::dieIfFileAccessFails($ofLong);
	unlink($ofLongTmp,$ofLongTmp2,$ofLongTmp3);
} elsif ($ENV{'seqType'} eq "chip") {
	my $shiftDistance = 100;
	bananas_agw::systemAndLogDieNonzero("$ENV{'bam2bed'} --do-not-sort < $ENV{'inputFile'} | cut -f 1-3,4,6,8 > $oft", $verbose, $syslog);
	bananas_agw::systemAndLogDieNonzero("$ENV{'sortBed'} $oft | awk '{x=\$2 + $shiftDistance; xu=x+1; if (\$5==\"-\") {x=\$3-$shiftDistance-1;xu=\$3-$shiftDistance}; print \$1\"\t\"x\"\t\"xu\"\tid\t\"\$5}' | uniq > $ofSorted", $verbose, $syslog);
} elsif ($ENV{'seqType'} eq "other") {
	bananas_agw::systemAndLogDieNonzero("$ENV{'bam2bed'} --do-not-sort < $ENV{'inputFile'} | cut -f 1-3,4,6,8 > $oft", $verbose, $syslog);
	bananas_agw::systemAndLogDieNonzero("$ENV{'sortBed'} $oft | awk '{x=\$2;xu=\$2+1; if (\$5==\"-\") {x=\$3-1;xu=\$3}; print \$1\"\t\"x\"\t\"xu\"\t\"\$4\"\t\"\$5}' > $ofSorted", $verbose, $syslog);
} else {
	confess qq{[ERROR] 05.tags.pl: The sequence type "$ENV{'seqType'}" is unrecognized / not currently supported. This is a PROGRAMMING ERROR and not a config file error! Fix it in 'bananas_agw.pm' or in 'bin/05.tags.pl'. };
}

bananas_agw::dieIfFileAccessFails($ofSorted);
bananas_agw::systemAndLogDieNonzero("/bin/mv -f $ofSorted $ofFinal", $verbose, $syslog); # Move the 'almost done' file to the final destination now that we know it is theoretically correct. This avoids having a good-looking file that is ACTUALLY broken.
bananas_agw::dieIfFileAccessFails($ofFinal);

if (-e $oft) { unlink($oft); } # clean up and remove a temp file, assuming everything else went well


