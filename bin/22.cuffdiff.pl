#!/bin/perl

my $syslog = "X.22.cuffdiff.syslog.txt";

use strict; use warnings; use Carp; # Carp = "confess" function
use File::Spec::Functions qw(catfile);   use Data::Dumper;   use File::Basename;
use lib "/wynton/group/gladstone/biocore/MonkeyPipeline"; use bananas_agw; # To make sure we find bananas_agw

# Procedure for this script:
# Run CuffDiff (note: not Cufflinks, even though you'd expect it!)
# Cuffdiff INTERNALLY runs cufflinks, and then:
#   If you tell it to do differential expression, it gives you all pair-wise P-values. BUT NOT THE COUNTS!
#   If you tell it to do only counts, it gives you FPKM, but not P-values, and not raw counts.
#   So you have to run it twice, and then run another ad-hoc script that Alex wrote to join the files together.
#   Downside: no raw counts provided (?). Maybe. Run htseq-counts to get those perhaps.

bananas_agw::requireEnvOrDie('cuffdiffExe', 'outputDir','gtfFile','labels','bamFilesByGroup');
my ($force, $verbose, $debug) = bananas_agw::envLooksTrue("force", "verbose", "debug");

my $COMMA_PLACEHOLDER_FROM_QSUB = "<COMMA>"; # <-- qsub HATES literal commas in arguments. We use this placeholder and manually replace it with a ',' below.
my $numThreads = (exists($ENV{pbs_ncpus}) && defined($ENV{pbs_ncpus})) ? $ENV{pbs_ncpus} : 1; # see if we got extra cpus specified
($numThreads =~ m/[1-9][0-9]*/) or confess "[PROGRAMMING ERROR]: We got a really WEIRD argument for 'pbs_ncpus'---it should be just a number, but we got this: $numThreads.";

my $cuffdiffExe        = $ENV{'cuffdiffExe'};
my $outDir             = $ENV{'outputDir'};                 # The output directory for differential expression results. Annoyingly does NOT provide counts!
my $outDirNoDiffExpr   = catfile(${outDir}, "1_individual_sample_counts_dir"); # The output directory for the cuffdiff run WITHOUT differential expression (only counts)
my $gtfFile            = $ENV{'gtfFile'}; # <-- required!
my $labels             = $ENV{'labels'};  # <-- Cuffdiff argument: must look like this: "DRUG,CONTROL,OTHERDRUG"
$labels                =~ s/${COMMA_PLACEHOLDER_FROM_QSUB}/,/g;          # <-- labels should be delimited by the special text "<COMMA>" since qsub CANNOT accept a literal comma for some reason, even when it's escaped.
my $bamFilesByGroup    = $ENV{'bamFilesByGroup'}; # <-- Cuffdiff argument: must look like this: "experiment.rep1.bam,experiment.rep2.bam  ctrlRep1.bam,ctrlRep2.bam,ctrlRep3.bam  otherReplicate1.bam,otherReplicate2.bam" (replicates are comma-delimited)
$bamFilesByGroup       =~ s/${COMMA_PLACEHOLDER_FROM_QSUB}/,/g; # <-- bamFilesByGroup should be delimited by the special text "<COMMA>" since qsub CANNOT accept a literal comma.
my $optionalGenomeFasta = $ENV{'genomeFasta'}; # <-- optional!
(-e $gtfFile)          or confess "[ERROR]: The input reference GTF file <$gtfFile> (which should be a FULL PATH to a GTF file) could not be found! Double check that it exists. This GTF-format file contains annotation data that is REQUIRED for Cuffdiff to generate usable output data in our pipeline.";
(-e $cuffdiffExe)      or confess "[ERROR]: Could not find the specified cuffdiff executable in <$cuffdiffExe>! Double check that it exists.";

if (defined($optionalGenomeFasta)) { (-e $optionalGenomeFasta) or confess "[ERROR]: The input 'genomeFasta' file was DEFINED (it was '$optionalGenomeFasta'), but does NOT EXIST! Double check this!";
} else {  print STDERR "[Warning]: Running Cuffdiff *without* doing --frag-bias-correct, since no 'genomeFasta' file was specified. Presumably there was no 'genomeFasta' entry in the monkey configuration file.\n"; }

# ========== SANITY CHECK OF SOME CUFFDIFF-SPECIFIC VARIABLES ========
($labels =~ /^[-_,.A-Za-z0-9]+$/) or confess "[ERROR]: The 'labels' argument ($labels) contained non-alphanumeric characters! It CANNOT contain any non-alphanumeric characters except for the dot, hyphen, or underscore! It should look something like this: TREATMENT,CTRL,DRUGX,OTHERDRUG";
my @labelGroupsArr = split(/,/, $labels);
my @bamGroupsArr = split(/ /, $bamFilesByGroup);
(scalar(@bamGroupsArr) == scalar(@labelGroupsArr)) or confess "[ERROR]: The 'labels' argument and 'bamFilesByGroup' argument must have the SAME number of groups, but they don't!\nlabels was: $labels\nand bamFilesByGroup was: $bamFilesByGroup\nThis is probably a coding error somewhere in bananas.pm. Error in cuffdiff was triggered";
# ========== SANITY CHECK OF SOME CUFFDIFF-SPECIFIC VARIABLES ========

my $logPath           = catfile(${outDir}, "0.log.txt"); # See if this file already exists in order to assess if Cuffdiff has been run already
my $logPathNoDiffExpr = catfile(${outDirNoDiffExpr}, "0.log.nodiff.txt"); # See if this file already exists in order to assess if Cuffdiff has been run already
if (-e "$outDir/gene_exp.diff" and -e "${outDirNoDiffExpr}/gene_exp.diff" and !$force) { # Output files exist already, AND we aren't forcing a re-run!
    ($verbose) && print STDERR "[OK -- SKIPPING] Skipping the re-run of Cuffdiff, because the output files <$outDir/gene_exp.diff> and <$outDirNoDiffExpr/gene_exp.diff> already exist, so we assume that Cuffdiff hopefully ran correctly.\n";
} else {
    # =========== [START] CHECK FOR tss_id and p_id in the GTF file ==========
    # "Note: Cuffdiff requires that transcripts in the input GTF be annotated with certain attributes in order to look for changes in primary transcript expression, splicing, coding output, and promoter use. These attributes are:"
    my $num_tss_id = `head -n 1000 $gtfFile | grep -c "tss_id"`; # Just check the first 1000 lines of the GTF for some tss_id field. That should be plenty to find at least one.
    my $num_p_id   = `head -n 1000 $gtfFile | grep -c "p_id"`;
    ($num_tss_id > 0) or print STDERR "[WARNING: SERIOUS POSSIBLE ERROR in GTF input to CUFFDIFF]! The input GTF file to cuffdiff ($gtfFile) did NOT have even a single 'tss_id' field. This field is REQUIRED in order to properly handle the GTF file. You can generate this field using cuffcompare, if necessary. See the docs at http://cufflinks.cbcb.umd.edu/manual.html#cuffdiff_input\n";
    ($num_p_id > 0)   or print STDERR "[WARNING: SERIOUS POSSIBLE ERROR in GTF input to CUFFDIFF]! The input GTF file to cuffdiff ($gtfFile) did NOT have even a single 'p_id' field. This field is REQUIRED in order to properly handle the GTF file. You can generate this field using cuffcompare, if necessary. See the docs at http://cufflinks.cbcb.umd.edu/manual.html#cuffdiff_input\n";
    # =========== [DONE] CHECK FOR tss_id and p_id in the GTF file ==========
    
    my $fragBiasCorrectStr = (defined($optionalGenomeFasta)) ? qq{ --frag-bias-correct="${optionalGenomeFasta}" } : qq{ };
    my $verboseStr = ($verbose) ? " --verbose " : " --quiet ";
    my $fdr = 0.05;
    
    my $cuffdiffParams = " $verboseStr "
	. " --num-threads=${numThreads} "
	. " --max-bundle-frags=2147483000 " # This prevents "too many reads" genes from being reported as ZERO due to 'HIDATA'
	. " --FDR=${fdr} "
	. " --upper-quartile-norm "
	. " --compatible-hits-norm "
	. " --multi-read-correct "
	. " ${fragBiasCorrectStr} ";

    # ====================== PART 1 of 2 for cuffdiff: calculate differential expression between groups ==========================
    bananas_agw::mkdirOrDie($outDir);
    bananas_agw::systemAndLog("touch ${logPath}", $verbose, $syslog);
    bananas_agw::systemAndLog("${cuffdiffExe} " . " --output-dir=${outDir} $cuffdiffParams --labels=${labels} " . " $gtfFile $bamFilesByGroup  2>&1 | tee --append ${logPath} ", $verbose, $syslog);
    #               EXECUTABLE             OPTIONS IN ANY ORDER                                          ARGUMENTS THAT *MUST* COME LAST
    # ====================== END OF PART 1 of 2 for cuffdiff =====================================================

    # ====================== PART 2 of 2 for cuffdiff: calculate per-file FPKM / counts ==========================
    # Now calculate PER FILE FPKM, and NOT differential expression at all!
    # Note the by default, cuffdiff only provides *average* FPKM values for EXPERIMENTAL GROUPS. It does not provide per-file data.
    # If you want that, you have to run it with a "fake" "--no-diff" comparison where each file is its own group. That's what's going on below
    my @bamPathsArray       = split(/[,\s]/, $bamFilesByGroup);    # Split on comma OR space! something like ("longpath/ctrl1.bam" "longpath/ctrl2.bam" "longpath/expr1.bam" "longpath/expr2.bam") (in an array)
    my @bamOnlyBasenames    = map{basename($_)} @bamPathsArray; # something like ("ctrl1.bam" "ctrl2.bam" "expr1.bam" "expr2.bam") (in an array)
    my $basenameLabelString = join(',', @bamOnlyBasenames); # something like "ctrl1.bam,ctrl2.bam,expr1.bam,expr2.bam" (a string)
    my $bamIndividualString = join(' ', @bamPathsArray);    # something like "longpath/ctrl1.bam longpath/ctrl2.bam longpath/expr1.bam longpath/expr2.bam" (a string)
    bananas_agw::mkdirOrDie($outDirNoDiffExpr);
    bananas_agw::systemAndLog("touch ${logPathNoDiffExpr}", $verbose, $syslog);
    bananas_agw::systemAndLog("${cuffdiffExe} " . " --no-diff " . " --output-dir=${outDirNoDiffExpr} $cuffdiffParams --labels=${basenameLabelString} " . " $gtfFile $bamIndividualString 2>&1 | tee --append ${logPathNoDiffExpr} ", $verbose, $syslog);
    #               EXECUTABLE                           OPTIONS IN ANY ORDER                                                                    ARGUMENTS THAT *MUST* COME LAST
    # ====================== END OF PART 2 of 2 for cuffdiff ==========================
}

