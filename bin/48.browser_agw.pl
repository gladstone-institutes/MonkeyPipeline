#!/usr/bin/perl

# Note: this file generates the non-binned (traditional) wiggle/bam tracks only.
# If you are looking for SEAN THOMAS's binned browser tracks (primarily intended for ChIP-seq), then check "08.browser.pl"

use strict; use warnings; use Carp; # Carp = "confess" function
use File::Spec::Functions qw(catfile);   use Data::Dumper;   use File::Basename;
use lib "/data/work/Code/alexgw/monkey_agw"; use bananas_agw; # To make sure we find bananas_agw

bananas_agw::requireEnvOrDie('binDir', 'outputDir','convertToGenomeBrowserExe','inputMultiFileString', 'studyName', 'serverBaseURL', 'rsyncDestinationBase');

my $force    = bananas_agw::envLooksTrue("force");
my $verbose  = bananas_agw::envLooksTrue("verbose");
my $dbug     = bananas_agw::envLooksTrue("debug");

my $studyName       = $ENV{'studyName'};
my $localOutDirPath = bananas_agw::removeTrailingSlash($ENV{'outputDir'}); # This is the staging area for the "08.browser.whatever" output directory, but NOT the output directory for the web server
my $successFilePath = catfile($localOutDirPath, "done.touch.txt");

my $trackInfoFilename = "all.txt"; # The filename of the 'info' file with all the track data
my $trackInfoLocalPath = catfile($localOutDirPath, $trackInfoFilename);

my $serverBaseURL      = bananas_agw::removeTrailingSlash($ENV{'serverBaseURL'}); # Note: should probably be https and not http!!!
my $projectSpecificURL = $serverBaseURL . "/" . $studyName;  # <-- Do NOT use "catfile" here, because it messes up the two slashes in "http://"
my $trackInfoURL       = $serverBaseURL . "/" . $studyName . "/" . $trackInfoFilename; # <-- Do NOT use "catfile" here, because it messes up the two slashes in "http://"

my $alsoCopyBamFiles   = bananas_agw::envLooksTrue("tracksShowEveryRead");
#print STDERR "Base URL is $serverBaseURL\n";
#print STDERR "project URL is $projectSpecificURL\n";

my $rsyncDestBase = bananas_agw::removeTrailingSlash($ENV{'rsyncDestinationBase'});
if (!$force and -d $localOutDirPath and -e $successFilePath) {
    ($verbose) && print STDERR "[OK] [Skipping re-run]: Looks like this browser directory ($localOutDirPath) ALREADY exists, so we are not re-running the browser track genereation.\n";
    exit(0);
}

(-d $ENV{'binDir'})                               or confess "The browser track 'bin directory' (where we will look for 'wigToBigWig' and 'genomeCoverageBed' was not found): directory was not specified or does not exist! Here is the (seemingly wrong) directory name: <$ENV{'binDir'}> (assuming it exists).";
(-e $ENV{'convertToGenomeBrowserExe'})            or confess "Cannot find the executable for converting to genome browser tracks. We expected it to have the following name and full path: '$ENV{'convertToGenomeBrowserExe'}'...";
(-f catfile($ENV{'binDir'}, "wigToBigWig"))       or confess "Cannot find 'wigToBigWig' in the binDir ($ENV{'binDir'}). Check to make sure it is there and can be executed!";
(-f catfile($ENV{'binDir'}, "genomeCoverageBed")) or confess "Cannot find 'genomeCoverageBed' in the binDir ($ENV{'binDir'}). Check to make sure it is there and can be executed!";

bananas_agw::mkdirOrDie($localOutDirPath); # Has to come BEFORE symlinking and checking the input files
bananas_agw::systemAndLog(qq{chmod a+rx $localOutDirPath}, $verbose); # Make sure everyone can read and 'execute' (list) this directory
# =============== MAKE SYMLINKS, AND CHECK TO MAKE SURE THE ELIGIBLE INPUT FILES ARE REALLY .bam or .sam files ============
my @eligibles = split(/\s/, $ENV{'inputMultiFileString'}); # Check each input bam/sam file. Input files are SPACE_DELIMITED!
for my $f (@eligibles) {
    # The file must be either a bam/sam file, OR possibly an empty line (^\s*$)
    (-e $f && ($f =~ /\.(sam|bam)/i) or $f =~ /^\s*$/) or confess "[ERROR] FILE MISSING. Error in generating browser tracks: could not find the specified file '$f' on the filesystem!... exiting now.";
    ($f =~ /^[\/]/) or confess "[ERROR] FILE PATH INVALID (NOT AN ABSOLUTE PATH). The input bam/sam files at this point should be FULL PATHS to the filenames--in other words, they must start with a '/' character! We CANNOT deal with relative paths! We also CANNOT deal with files with *spaces* of any sort in them! Fix this. The offending filename was: '$f'.";
    bananas_agw::systemAndLog(qq{ln -s ${f} $localOutDirPath/}, $verbose); # Create a symlink in the new browser directory. This is needed for the conversion script to work!
}
# ======================== DONE CHECKING ===================

my $BAMS_ALREADY_SORTED = 1; # Assume this is true. It is if we use tophat, but NOT bowtie, so beware...
my $sort_arg     = ($BAMS_ALREADY_SORTED) ? "--nosort" : "--sort"; # note: right now we ALWAYS assume the files are sorted!
my $keep_bam_arg = ($alsoCopyBamFiles)    ? "--bam"    : "--nobam"; # If we do NOT want the bam files to be retained, then specify "--nobam".

# Note 'nosort' above---we assume the input files are ALREADY SORTED! This is true for Tophat but not for Bowtie!

my $BROWSER_CMD = qq{PATH=$ENV{'binDir'}:\$PATH && cd $localOutDirPath && $ENV{'convertToGenomeBrowserExe'}  --url=\"$projectSpecificURL\" $sort_arg $keep_bam_arg  *.[sbSB][aA][mM] }; # note how we ADD the 'binDir' to the path to ensure that certain executables (wigToBigWig and genomeCoverageBed) can be found by the script! This is REQUIRED!. The weird [sBSB...] thing at the end finds all SAM and BAM files.
bananas_agw::systemAndLog($BROWSER_CMD, $verbose);

($verbose) && print STDERR "[OK] Note: assuming that the input BAM files for the browser tracks are SORTED. Which is true for Tophat-aligned but not Bowtie-aligned files.\n";
($verbose) && print STDERR "[OK] Now copying the files for project URL '$projectSpecificURL' to the output directory.\n";

my $destNoSlash  = $rsyncDestBase . "/" . $studyName; # <-- DO NOT USE 'catfile' here --- The specific destination directory for this particular project. Do NOT use 'catfile' for this!

my $shouldCopyFilesToBrowserDestination = not bananas_agw::isNA($rsyncDestBase); # We are rsyncing the files UNLESS the rsync destination is 'NA'
my $looksLikeLocalCopying               = ($rsyncDestBase !~ '[:@]'); # if there's no ':' or '@' in the rsync command, then this is PROBABLY a local copy operation, so we can make a symlink to that directory
if ($looksLikeLocalCopying && $shouldCopyFilesToBrowserDestination) { bananas_agw::systemAndLog("ln -s $destNoSlash ${localOutDirPath}.server", $verbose); } # Make a symlink to the served-files directory

bananas_agw::systemAndLog(qq{   echo "## Paste this link into the Genome Browser custom track box: $trackInfoURL"                                              > $trackInfoLocalPath ; }
			  . qq{ echo "## Remember to set the default track location to an area with data (e.g. 'Gapdh' for human, 'Actn1' for mouse)"         >> $trackInfoLocalPath ; }
#			  . ($looksLikeLocalCopying ? qq{ echo "## Tracks were copied to this filesystem location: $destNoSlash"         >> $trackInfoLocalPath ; } : " ") # If it's a local copy, then we can safely say where the files were copied to. If it's remote, we shouldn't reveal the username.
			  . qq{ echo "## Finally, in the Genome Browser, go to My Data --> Session --> Save Session and save it as (for example): $studyName" >> $trackInfoLocalPath ; }, $verbose);

if ($shouldCopyFilesToBrowserDestination and not($trackInfoURL =~ /^https:[\/][\/]/i)) {
    warn("WARNING: the URL for the 'info' file (which was $trackInfoURL) for the genome browser DID NOT start with 'https://', but we expected it to. Note: HTTPS and not just HTTP. This is crucial for our system since it has issues serving large files over http.");
    bananas_agw::systemAndLog(qq{echo "## WARNING: The URL above ($trackInfoURL) *DID NOT* start with 'https://'--that is probably a mistake." >> $trackInfoLocalPath });
}

bananas_agw::systemAndLog(qq{echo "## Note: that link should probably be an 'https' link. Paste this link into the Genome Browser: $trackInfoURL" >> $trackInfoLocalPath});
bananas_agw::systemAndLog(qq{cat ${localOutDirPath}/Browser*.txt >> $trackInfoLocalPath }, $verbose); # Create the track info file!

(-e $trackInfoLocalPath) or confess "Failure to create the local track info file, which we tried (and failed) to write to $trackInfoLocalPath";

bananas_agw::systemAndLog("chmod a+r $trackInfoLocalPath $localOutDirPath $localOutDirPath/Browser*", $verbose); # Make sure everyone can READ the browser files

if ($shouldCopyFilesToBrowserDestination) {
    bananas_agw::systemAndLog("rsync --copy-links --chmod=u+rx,g+rx,o+rx --perms --dirs $localOutDirPath                              ${destNoSlash}" , $verbose); # Make a new (possibly remote) directory.
    bananas_agw::systemAndLog("rsync --copy-links --chmod=u+r,g+r,o+r    --perms        $trackInfoLocalPath $localOutDirPath/Browser* ${destNoSlash}/", $verbose); # Make a new (possibly remote) directory for this project. You have to use '-r' or it omits the directory.
} else {
    print STDERR "[NOTE]: Not rsync-ing, because the rsync destination was set to \"" . $bananas_agw::NA_VALUE . "\", which means that we are not copying anything to a destination directory.\n";
}

# This must be at the VERY END: we create a 'success' file!
bananas_agw::systemAndLog(qq{echo 'This file is used to indicate that the browser generation ran correctly and does not need to be re-run.' > $successFilePath}, $verbose);
