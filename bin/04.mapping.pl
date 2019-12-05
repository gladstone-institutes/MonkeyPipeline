#!/usr/bin/perl -w
my $syslog = "X.04.mapping.syslog.txt";

use strict; use warnings; use Carp; # Carp = "confess" function
use File::Spec::Functions qw(catfile);   use Data::Dumper;   use File::Basename;
use lib "/wynton/group/gladstone/biocore/MonkeyPipeline"; use bananas_agw; # To make sure we find bananas_agw

bananas_agw::requireEnvOrDie('alignerPath','alignerExeWithoutPath','samtools','minMapQ','filterDir','mappingDir','sampleName','inputFile1','finalBam','sortMethod');

sub getSampleFileName($$) { my ($rDir, $file) = @_; return catfile($rDir,(bananas_agw::noExtensions($file) . ".gz")); }
my ($force, $verbose) = bananas_agw::envLooksTrue("force", "verbose");
my $finalBam    = $ENV{'finalBam'}; #catfile($ENV{'mappingDir'},"$sampleName_${genome}_q${mapq}.bam"); # example of what this looks like: mysample_mm9_q30.bam
my $isPaired    = (defined($ENV{'inputFile2'}) && (length($ENV{'inputFile2'}) > 0) && (uc($ENV{'inputFile2'}) ne "NA")); # if inputFile2 is defined, then we assume we DO have a paired sample!
my $fqFile1     =               getSampleFileName($ENV{'filterDir'},$ENV{'inputFile1'});
my $fqFile2     = ($isPaired) ? getSampleFileName($ENV{'filterDir'},$ENV{'inputFile2'})  :  ""; # make sure it is EMPTY if the file is not paired!
my $mapq        = $ENV{'minMapQ'}; # should be a number from 0 to ... well, it's different for each aligner. STAR: 3 means "unique" probably. TOPHAT, set to 30+ if you want unique.
my $sampleName  = $ENV{'sampleName'};
my $gtf         = exists($ENV{'gtfFile'}) ? $ENV{'gtfFile'} : undef; # <-- A GTF annotation file is used by Tophat and STAR, but not by Bowtie! Can be "undef" or NA!

# should get the 'pbs_mem' variable and use that to set the 'limitGenomeGenerateRAM' variable as well
my $numThreads  = exists($ENV{'pbs_ncpus'}) ? $ENV{'pbs_ncpus'} : 1; # Default number of threads is 1 if the 'pbs_ncpus' parameter wasn't explicitly specified
$numThreads =~ m/^(\d)+$/ or confess "[ERROR] 04.mapping.pl. numThreads (pbs_ncpus) must be an integer, whereas you passed in this value: $ENV{'pbs_ncpus'} <-- non-integer value";

my $alignerOutCode;
my $outSam;

sub samtoolsIndex($) {
	my ($inBam) = @_;
	(-e $inBam) or confess "[ERROR] 04.mapping.pl: input file to samtools index <$inBam> could not be found! Cannot index a nonexistent BAM file.";
	bananas_agw::systemAndLog("$ENV{samtools} index ${inBam}", $verbose, $syslog); # <-- Generate an index (bam.bai) file, which is OCCASIONALLY needed (for example, for RSEQC and for IGV)
	(-e "${inBam}.bai") or confess "[ERROR] 04.mapping.pl: in samtools indexing: failed to created the expected '.bai' index file named <${inBam}.bai>. Maybe there was a problem finding samtools?";
}

# ========== Check prior completion and exit if the files exist ========
if (-e $finalBam and !$force) { # Does the file ALREADY exist AND we aren't forcing a re-run? Then skip it.
	($verbose) && print STDOUT "[OK] [SKIPPING RE-RUN]: the already-aligned file '$finalBam' (for sample '$sampleName') already exists, so we will not re-align it.\n";
	(-e "${finalBam}.bai") or samtoolsIndex(${finalBam}); # generate .bai file if it doesn't exist
	exit(0); # This is OK, exit without an error!
}

my $BYCOORD = "BY_COORD";
my $BYNAME  = "BY_NAME";

# =========== Check a bunch of error conditions ===========
(uc($ENV{sortMethod}) eq $BYNAME or uc($ENV{sortMethod}) eq $BYCOORD) or confess "[ERROR] 04.mapping.pl: The 'sortMethod' parameter must be either the literal text '$BYNAME' or '$BYCOORD'. It cannot be ANYTHING else. You specified '$ENV{sortMethod}', however. Fix this!";
(-e $fqFile1)                                               or confess "[ERROR] 04.mapping.pl: Cannot find the input fastq file 1 ('fqFile1') named '$fqFile1'! Check to make sure this file really exists on the filesystem. Also make sure it is a FULL path, not a relative path.";
($mapq >= 0 and $mapq <= 999 and ($mapq == int($mapq)))     or confess "[ERROR] 04.mapping.pl: mapq was set to the weird value of '$mapq'--but it should really be an INTEGER in the range from 0 and ... maybe 50 at most? (Probably lower, like 30).";
bananas_agw::mkdirOrDie($ENV{'mappingDir'});

sub samToFilteredAndSortedBam($$$$) {
	# Takes the '$inSam' file, sorts it by the required method, and generates '$outBam'
	my ($inSam, $outBam, $sortBy, $bamLog) = @_;
	(-e $inSam) or confess "[ERROR] 04.mapping.pl: an earlier step, before sorting, somehow FAILED to generate the file '$inSam'. Something went wrong. Fix the mapping script!";
	my $bamTmp        = "${outBam}_TEMP_DELETE_ME_SOON_bam_conversion_in_progress_tmp.bam";
	my $almostDoneBam = "${outBam}_TEMP_DELETE_ME_SOON_final_conversion_in_progress.tmp.bam"; # the next-to-last step in the conversion. This way if we see the final file actually made, we KNOW it was the full file, and not a partial write.
	if (uc($sortBy) eq uc($BYCOORD)) { 		# <-- POSITIONALLY SORTED (COORDINATE SORTED output bam)
		bananas_agw::systemAndLog("$ENV{samtools} view -q ${mapq} -S -h -u -b ${inSam} > ${bamTmp}       2> ${bamLog};    $ENV{samtools} sort ${bamTmp} -o ${almostDoneBam}  2>>  ${bamLog}", $verbose, $syslog);
	} elsif (uc($sortBy) eq uc($BYNAME) ) {
		bananas_agw::systemAndLog("$ENV{samtools} view -q ${mapq} -S -h    -b ${inSam} > ${almostDoneBam} 2> ${bamLog}", $verbose, $syslog); print STDERR "WARNING: lexographic sorting is usually a bad idea and will probably break things downstream! It will for sure break RSeQC, and most likely other tools!\n";
	}			# <-- LEXOGRAPHICALLY SORTED output bam
	bananas_agw::systemAndLog("/bin/mv -f ${almostDoneBam} ${outBam}", $verbose, $syslog); # move it to the final destination file
	if (-e $outBam) {
		(-e $bamTmp) && unlink($bamTmp);
	} else {
		bananas_agw::systemAndLog("touch ${outBam}_FAILED_FOR_SOME_REASON_CHECK_LOGS", $verbose, $syslog);
		confess "[ERROR] 04.mapping.pl: Failed to create a final BAM file ('$outBam')! Sometimes this happens if there is something wrong in the samtools sorting, so maybe check that first.";
	}
}

my $alignerLog = "${finalBam}.$ENV{alignerExeWithoutPath}.log"; # WARNING: Do not change this name without changing the code in "05.xSummary.pl" to look for this updated filename! There is a variable named 'bowtieStatsLogFilename' that MUST have the same name. And it REQUIRES bowtie!


if ($ENV{alignerExeWithoutPath} =~ m/^(tophat|tophat2)/i) { # ~~~~~~~~ TOPHAT (actually Tophat2) ~~~~~~~~
	bananas_agw::requireEnvOrDie('bowtie2Index');
	($ENV{sortMethod} eq $BYCOORD) or confess "With tophat, we actually only support '$BYCOORD', not '$BYNAME' sorting.";
	(defined($gtf)) or confess "FATAL ERROR in 4.mapping.pl: Tophat needs a GTF annotation file. However, the 'gtfFile' line in your monkey config file did NOT specify one! If you do not have such a file, set it to 'gtfFile = NA' to not use a GTF file at all.";
	((uc($gtf) eq "NONE") or (uc($gtf) eq "NA") or (-e $gtf)) or confess "[FATAL ERROR] 04.mapping.pl: The following GTF annotation file, which was specified in the config file, and is MANDATORY for Tophat, doesn't exist: \"$gtf\"";
	my $dirThatContainsBowtie = dirname($ENV{alignerPath}); # the "dirname" function is located in File::Basename. This is IMPORTANT to include on the tophat command, otherwise Tophat will not locate Bowtie properly if it isn't in the default path (which it might not be!).
	my $subdirName            = "${sampleName}.tophat.dir";	# NOT the path
	my $topOutDir             = catfile($ENV{'mappingDir'}, "${subdirName}");
	my $tophat_GTF_string     = ((uc($gtf) eq "NONE") or (uc($gtf) eq "NA")) ? " " : " --GTF=$gtf ";
	my $pairedMateDist        = 50; # <-- this is MANUALLY specified here for now... turns out basically not to matter
	my $tophatArgs            = " ${tophat_GTF_string} --min-anchor=5 --segment-length=25 --no-coverage-search --segment-mismatches=2 --splice-mismatches=2 --num-threads=$numThreads " . (($isPaired) ? " --mate-inner-dist=$pairedMateDist " : " ");
	my $topOutBamInSubdir     = catfile($topOutDir, "accepted_hits.bam"); # <-- tophat makes a file with this specific name, so DO NOT CHANGE IT!
	my $filteredBamNameOnly   = "accepted_filtered_q${mapq}.bam";
	my $filteredBam           = catfile($topOutDir, $filteredBamNameOnly); # <-- this is the final output file of interest with only the reads that met the 'mapq' threshold
	my $mapqLog               = catfile($topOutDir, "mapq.samtools.filtering.txt.log");
	bananas_agw::systemAndLog("PATH=\$PATH:${dirThatContainsBowtie} &&  $ENV{alignerPath} -o $topOutDir $tophatArgs $ENV{'bowtie2Index'} $fqFile1 $fqFile2 2> ${alignerLog}", $verbose, $syslog); # <-- Note that the PATH is manually specified in order to add the directory with Bowtie in it, to prevent Tophat from failing to find Bowtie under certain circumstances.
	(-e $topOutBamInSubdir) or confess "[ERROR] 04.mapping.pl: the 'top out bam in subdir' file failed to be generated to the expected location of: $topOutBamInSubdir";
	$alignerOutCode = bananas_agw::systemAndLog("$ENV{samtools} view -q $mapq -h -b $topOutBamInSubdir > $filteredBam 2> $mapqLog", $verbose, $syslog); # filter for the minMapQ
	my $filteredBamLocalPath = catfile($subdirName, $filteredBamNameOnly); # NOT an absolute path!!
	bananas_agw::systemAndLog("/bin/ln -f -s $filteredBamLocalPath $finalBam", $verbose, $syslog); # <-- note: in the tophat case, the 'finalBam' is a symlink (this is OK). The 'real' non-symlink file is $filteredBam in the tophat directory.

} elsif ($ENV{alignerExeWithoutPath} =~ m/^(bowtie|bowtie2)/i) { # ~~~~~~~~ BOWTIE (actually Bowtie2) ~~~~~~~~
	# ==================================================== BOWTIE2 ===============================
	bananas_agw::requireEnvOrDie('bowtie2Index');
	$outSam  = catfile($ENV{'mappingDir'}, "${sampleName}.sam"); # example of what this looks like: mysample.sam.
	my $FILE_PAIR_INPUT_TEXT = ($isPaired) ? " -1 $fqFile1 -2 $fqFile2 " : " -U $fqFile1 ";
	my $samTmp  = "${outSam}.temp_mapping_in_progress.sam";
	$alignerOutCode = bananas_agw::systemAndLog("$ENV{alignerPath} --threads=$numThreads -x $ENV{'bowtie2Index'} $FILE_PAIR_INPUT_TEXT -S $samTmp   2>  ${alignerLog}", $verbose, $syslog);
	bananas_agw::systemAndLog("/bin/mv -f $samTmp $outSam", $verbose, $syslog); # Move the temp file to the final destination now that we know it is theoretically correct. This avoids having a good-looking sam file that is actually broken.
	# ==================================================== END BOWTIE2 ===========================
} elsif ($ENV{alignerExeWithoutPath} =~ m/^(bwa)/i) {
	# ==================================================== BWA ==============================
	bananas_agw::requireEnvOrDie('bwaIndex');
	$outSam = catfile($ENV{'mappingDir'}, "${sampleName}.sam"); # example of what this looks like: mysample.sam.
	# Note below: BWA is aligned with the 'mem' command and not the 'aln' command as you might expect. See details: https://www.biostars.org/p/117225/
	# From that site: "In a test (recovery of alignments):
	# bowtie2: 30% / bwa aln: 25% / bwa mem: 85%  " <-- so use this 85% one, of course! Here's the full paper: http://arxiv.org/abs/1303.3997
	my $unalignedReadsBam = "${finalBam}"; $unalignedReadsBam =~ s/[.]bam$/_unaligned.bam/i;
	my $samTmp            = "${outSam}.temp_mapping_in_progress.sam";
	my $alignedAndUnSam   = "${samTmp}.both_aligned_and_not--temp.sam";
	$alignerOutCode = bananas_agw::systemAndLog("$ENV{alignerPath} mem   -t $numThreads   -R \'\@RG\\tID:${sampleName}\\tSM:${sampleName}\\tLB:${sampleName}\\tPL:ILLUMINA\' $ENV{bwaIndex}  $fqFile1 $fqFile2 > $alignedAndUnSam 2> ${alignerLog}", $verbose, $syslog);
	(-e $alignedAndUnSam) or confess "[ERROR] 04.mapping.pl: the 'aligned and unaligned sam' file failed to be generated (by bowtie) to the expected location of: $alignedAndUnSam";
	bananas_agw::systemAndLog("$ENV{samtools} view -S -h -F 0x4    $alignedAndUnSam > $samTmp"); #   aligned reads ONLY. It's a SAM file, **not** a bam file! We convert it to BAM below
	bananas_agw::systemAndLog("$ENV{samtools} view -S -h -f 0x4 -b $alignedAndUnSam > $unalignedReadsBam"); # unaligned reads ONLY. This is a BAM file!
	(-e $samTmp and -e $unalignedReadsBam) or confess "[ERROR] 04.mapping.pl: after the aligner $ENV{alignerExeWithoutPath} generated the aligned file '$alignedAndUnSam', samtools (at $ENV{samtools}) somehow was unable to split the file into the aligned and unaligned reads! This means that at least one of the expected files were NOT present at the following locations: either\n$samTmp\nor possibly:\n$unalignedReadsBam";
	unlink($alignedAndUnSam); # <-- this is the only temporary file. Do NOT get rid of $unalignedReadsBam!
	bananas_agw::systemAndLog("/bin/mv -f $samTmp $outSam", $verbose, $syslog); # Move the temp file to the final destination now that we know it is theoretically correct. This avoids having a good-looking sam file that is actually broken.
	# ==================================================== END BWA ===========================
} elsif ($ENV{alignerExeWithoutPath} =~ m/^(star)/i) {
	# ==================================================== STAR ==============================
	# STAR can actually directly generate output BAM files now! (This is a new-ish feature.): see https://github.com/alexdobin/STAR/releases   http://seqanswers.com/forums/showthread.php?t=46592
	bananas_agw::requireEnvOrDie('starIndexDir');
	(defined($gtf)) or confess "FATAL ERROR in 4.mapping.pl: STAR needs either a GTF annotation file, or the GTF file to be explicitly specified as 'NA'. However, the 'gtfFile' line in your monkey config file did NOT specify one! If you do not have such a file, set it to 'gtfFile = NA' to not use a GTF file at all.";
	my $gtfCmd = bananas_agw::fileIsReadableAndNonzero($gtf) ? " --sjdbGTFfile $gtf " : " "; # either supply a gtf --- or there is no GTF
	my $starOutDir = catfile($ENV{'mappingDir'}, "${sampleName}.star.dir"); # full path! Creates a new directory.
	$outSam        = catfile($starOutDir, "Aligned.out.sam"); # <-- STAR has hard-coded filenames, and this is one of them. We will convert this to a BAM file and then DELETE the original SAM file, so beware!
	my $starOutSortedBam = catfile($starOutDir, "Aligned.sortedByCoord.out.bam"); # <-- STAR has hard-coded filenames
	bananas_agw::mkdirOrDie($starOutDir);
	$alignerOutCode = bananas_agw::systemAndLog(qq{$ENV{alignerPath}  --runThreadN $numThreads  --genomeDir $ENV{'starIndexDir'} }
						    . qq{ --readFilesCommand gunzip -c } # <-- incredibly, this works somehow, despite the lack of a quotation marks, and a space (gunzip -c is totally valid with no quotes)
						    . qq{ --readFilesIn $fqFile1 $fqFile2 }
						    . qq{ --outFileNamePrefix $starOutDir/ }
						    . qq{ --outSAMtype BAM SortedByCoordinate } # outputs a file named "Aligned.sortedByCoord.out.bam"
						    # limitGenomeGenerateRAM = 31000000000  # default is 31 GB!!!  1000000000 * num_gigabytes
						    # outTmpKeep = None
						    # outReadsUnmapped = None
						    . qq{ $gtfCmd }
						    . qq{ --alignSJDBoverhangMin 3 } # --alignSJDBoverhangMin is used at the mapping step to define the minimum allowed overhang over splice junctions. For example, the default value of 3 would prohibit overhangs of 1b or 2b.
						    , $verbose, $syslog);
	(0 == $alignerOutCode) or confess "[ERROR] 04.mapping.pl: the aligner $ENV{alignerExeWithoutPath} returned the following non-zero exit code: $alignerOutCode . Something went wrong. Fix the mapping script!";
	(-e $starOutSortedBam) or confess "[ERROR] 04.mapping.pl: the 'sorted bam' file failed to be generated (by STAR) to the expected location of: $starOutSortedBam";
	bananas_agw::systemAndLog("/bin/ln -f -s $starOutSortedBam $finalBam", $verbose, $syslog); # <-- note: in the tophat case, the 'finalBam' is a symlink (this is OK). The 'real' non-symlink file is $filteredBam in the tophat directory.
	# ==================================================== END STAR ==========================
} else {
	confess "[ERROR] 04.mapping.pl: Error in MAPPING READS: although the only supported aligners are 'tophat', 'bowtie', 'bwa', and 'star', but yours was <$ENV{alignerExeWithoutPath}>, which was allegedly at the following full path: <$ENV{alignerPath}>. Check to make sure that is IN FACT exactly the correct aligner and that the full path is right! Watch out for capitalization, permissions, and misspellings!";
}

(0 == $alignerOutCode) or confess "[ERROR] 04.mapping.pl: the aligner $ENV{alignerExeWithoutPath} returned the following non-zero exit code: $alignerOutCode . Something went wrong. Fix the mapping script!";

if (-e ${finalBam}) {
	# I guess the final bam file ALREADY exists, so no need to try to convert SAM -> BAM here
} else {
	# If there's no final BAM file, then we expect a *SAM* file to exist, which we will convert to e a BAM file with 'samToFilteredAndSortedBam'!
	(-f $outSam)           or confess "[ERROR] 04.mapping.pl: somehow, the aligner $ENV{alignerExeWithoutPath} failed to generate the expected output SAM file (not BAM file!) named <$outSam>.";
	samToFilteredAndSortedBam($outSam, $finalBam, $ENV{sortMethod}, "${finalBam}.conversion_log.txt"); # Generates the $finalBam BAM file AND filters by minimum map quality from the input $outSam SAM file
}

if (-e ${finalBam} && -e ${outSam}) {
	unlink($outSam); ($verbose) && print STDERR "[OK] 04.mapping.pl: since mapping with $ENV{alignerExeWithoutPath} finished successfully, we DELETED the redundant SAM file <$outSam>.\n"; # <-- delete the output sam file after it gets converted to a bam file!
}

# By now, the FINAL bam file should exist!
if (not -e $finalBam) {
	bananas_agw::systemAndLog("touch ${finalBam}_ALIGNMENT_WITH_$ENV{alignerExeWithoutPath}_FAILED_FOR_SOME_REASON_CHECK_LOGS", $verbose, $syslog); # make a 'look, it's super broken!' file
	confess "[ERROR] 04.mapping.pl: Failed to create the final BAM file ('$finalBam')! Sometimes this happens if there is something wrong in the samtools sorting, so maybe check that first.";
}

samtoolsIndex(${finalBam}); # generate .bai file and check to make sure it exists




# This isn't used, but it's totally possible that we could create indexes from an input genome that was not already indexed.
#  sub makeBwaIndex() {
#  	# ==================== not used at all / yet!! ============
#  	# check that bwa index exists / make it / or fail if unable to make it
#  	# NOT TESTED
#  	# Example of how to call it:
#  	# my $indexPrefix = $toolDir . "/bwaIndex_" . $cfg->{expID};
#  	# checkBwaIndex($indexPrefix, $cfg->{bwa}, $cfg->{genomeFasta});
#  	my $bPre = "what, this is undefined";
#  	my ($indexPrefix, $bwa_exe, $genomeFasta) = @_;
#  	my $bp1  = "${bPre}.amb";
#  	my $bp2  = "${bPre}.ann";
#  	my $bp3  = "${bPre}.bwt";
#  	my $bp4  = "${bPre}.pac";
#  	my $bp5  = "${bPre}.sa";
#  	if (! -e $bp1 || ! -e $bp2 || ! -e $bp3 || ! -e $bp4 || ! -e $bp5) {
#  		# if index doesn't exist, make it
#  		system("${bwa_exe} index -p ${indexPrefix} ${genomeFasta} "); #>& /dev/null ");
#  		# if we try to make it and it still doesn't exist, then die...
#  		if (! -e $bp1 || ! -e $bp2 || ! -e $bp3 || ! -e $bp4 || ! -e $bp5) {
#  			confess "cannot generate bwa index ($indexPrefix) with given permissions";
#  		}
#  	}
#  }
