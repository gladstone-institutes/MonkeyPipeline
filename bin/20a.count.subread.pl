#!/usr/bin/perl

use strict; use warnings; use Carp; # Carp = "confess" function
use File::Spec::Functions qw(catfile);   use Data::Dumper;   use File::Basename;
use lib "/wynton/group/gladstone/biocore/MonkeyPipeline"; use bananas_agw; # To make sure we find bananas_agw

bananas_agw::requireEnvOrDie('subreadFeatureCountsExe', 'subreadCountsFilename', 'outputDir','gtf','singleEndBamFilesSpaceDelim','pairedEndBamFilesSpaceDelim'); # <-- these vars MUST be defined in the %ENV hash
my ($force, $verbose, $debug) = bananas_agw::envLooksTrue("force", "verbose", "debug");

my $EXON_STR    = "exon";
my $GENE_ID_STR = "gene_id";
my $NUM_THREADS = (exists($ENV{pbs_ncpus}) && defined($ENV{pbs_ncpus})) ? $ENV{pbs_ncpus} : 1; # see if we got extra cpus specified

($NUM_THREADS =~ m/[1-9][0-9]*/) or confess "[PROGRAMMING ERROR]: We got a really WEIRD argument for 'pbs_ncpus'---it should be just a number, but we got this: $NUM_THREADS.";

($ENV{'subreadCountsFilename'} =~ m/[-_.A-Za-z.]/) or confess "[ERROR]: The output subread filename contained some weird characters besides A-Z, 0-9, ., _, and -. This is NOT ALLOWED. The offending filename was: <$ENV{'subreadCountsFilename'}>.";


# Note: subread featureCounts makes a whole bunch of ~20 megabyte temp files with names like:
#   temp-core-072447-00E0ED5D1ADC.sam-TH03-BK000004.tmp
#   temp-core-072447-00E0ED5D1ADC.sam-TH04-BK000000.tmp
#   temp-core-072447-00E0ED5D1ADC.sam-TH04-BK000001.tmp (total size of all of them for a large run was ~2 gigabytes)
# These are dumped into the WORKING directory, so the first thing we do for SUBREAD is to 'cd' into the final directory so that these don't litter the user's home directory.

my $featureCountsExe   = $ENV{'subreadFeatureCountsExe'}; (-e $featureCountsExe)  or confess "[ERROR]: Could not find the specified executable in <$featureCountsExe>! Double check that it exists.";
my $outDir             = $ENV{'outputDir'};                 # The output directory for differential expression results. Annoyingly does NOT provide counts!
my $gtf                = $ENV{'gtf'}; (-e $gtf)               or confess "[ERROR]: The input reference GTF file <$gtf> (which should be a FULL PATH to a GTF file) could not be found! Double check that it exists.";
my $finalOut           = catfile($outDir, $ENV{'subreadCountsFilename'});
my $logPath            = catfile($outDir, "subread.count.log");

my $omitAmbigFinalOut           = catfile($outDir, "omit.ambiguous." . $ENV{'subreadCountsFilename'});
my $omitAmbigLogPath            = catfile($outDir, "omit.ambiguous.subread.count.log");
my $allBamStr          = "$ENV{'singleEndBamFilesSpaceDelim'} $ENV{'pairedEndBamFilesSpaceDelim'}"; # single end and paired end files together.
$allBamStr =~ s/^\s+|\s+$//g; # remove leading and trailing whitespace

foreach my $f (split(/\s+/, "$allBamStr")) { (!defined($f) or $f eq '') or bananas_agw::dieIfFileAccessFails($f, "The .bam file '$f' was specified, but it cannot be found at that path. The list of ALL bam files is: $allBamStr\nCheck the inputs in bananas.pm"); } # Sanity check: the specified BAM files must exist, if they aren't blank

if (bananas_agw::fileIsReadableAndNonzero("$finalOut") and !$force) { # Output files exist already, AND we aren't forcing a re-run!
    ($verbose) && print STDERR "[OK] [SKIPPING] Skipping the re-run of SubRead's FeatureCounts program, because the output file <$finalOut> already exists, so we assume that featureCounts ran correctly.\n";
    exit(0);
}


# =========== [START] CHECK FOR tss_id and p_id in the GTF file ========
# "Note: subread requires that transcripts in the input GTF be annotated with certain attributes in order to look for changes in primary transcript expression, splicing, coding output, and promoter use. These attributes are:"
my $num_exon    = `head -n 1000 $gtf | grep -c "$EXON_STR"`; # Just check the first 1000 lines of the GTF for this field.
my $num_gene_id = `head -n 1000 $gtf | grep -c "$GENE_ID_STR"`;
($num_exon    > 0) or confess("[FAILURE: SERIOUS LIKELY ERROR in GTF input to SUBREAD]! The input GTF file ($gtf) did NOT have even a single '$EXON_STR' field (usually 'exon').    This field is REQUIRED in order to properly handle the GTF file.\n");
($num_gene_id > 0) or confess("[FAILURE: SERIOUS LIKELY ERROR in GTF input to SUBREAD]! The input GTF file ($gtf) did NOT have even a single '$GENE_ID_STR' field (usually 'gene_id'). This field is REQUIRED in order to properly handle the GTF file. n");
# =========== [DONE] CHECK FOR tss_id and p_id in the GTF file =========
bananas_agw::mkdirOrDie($outDir);
bananas_agw::systemAndLog(qq{touch "${logPath}"}, $verbose);
# We are NOT using the '-f' option: -f: count at "feature" (normally this means EXON) level instead of gene level
# -O: assign reads to ALL overlapping features
# -M: assign multi-mapping reads to all features as well
# -Q: required quality minimum (0 is our default)

# --primary: only include primary hits
bananas_agw::systemAndLog(qq{cd "$outDir"}, $verbose); # <-- Causes the temp files to be written in $outDir instead of the user's home directory (which is the default)
bananas_agw::systemAndLog("${featureCountsExe} -p -O -M -T ${NUM_THREADS}   -t $EXON_STR -g $GENE_ID_STR -a $gtf  -o ${finalOut}          $allBamStr  2>&1 | tee --append ${logPath} "         , $verbose);
# Below: the original way we ran this command before: (Note that it's ok if '-p' is included here, because it just doesn't have any effect on single-end reads.)
bananas_agw::systemAndLog("${featureCountsExe} -p       -T ${NUM_THREADS}   -t $EXON_STR -g $GENE_ID_STR -a $gtf  -o ${omitAmbigFinalOut} $allBamStr  2>&1 | tee --append ${omitAmbigLogPath} ", $verbose);

(-e $finalOut) or confess "Failed to generate output counts file <$finalOut>...";
($verbose) && print STDERR "[OK] We have now generated the final subreads 'counts' file in the file <${finalOut}>.\n";

