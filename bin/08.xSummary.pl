#!/bin/perl
eval 'exec /bin/perl -w -S $0 ${1+"$@"}'
if 0; # not running under some shell

use strict; use warnings; use Carp; # Carp = "confess" function
use File::Spec::Functions qw(catfile);   use Data::Dumper;   use File::Basename;
use lib "/wynton/group/gladstone/biocore/MonkeyPipeline"; use bananas_agw; # To make sure we find bananas_agw

bananas_agw::requireEnvOrDie('rsyncDir','browserDir','trackFileSuffix'); # <-- these vars MUST be defined in the %ENV hash
my $force    = bananas_agw::envLooksTrue("force");
my $verbose  = bananas_agw::envLooksTrue("verbose");
my $dbug     = bananas_agw::envLooksTrue("debug");

my $trackSuff                   = $ENV{'trackFileSuffix'}; # something like ".tracks.txt" probably
my $tempTrackFile               = catfile($ENV{'browserDir'}, "_______TEMP_SOON_MOVE_THIS_FILE______"); # this must NOT contain the text above, or it will be removed!
my $finalSummaryTracksFile      = catfile($ENV{'browserDir'}, "0.allTracksUCSC.txt"); # this must NOT contain the text above, or it will be removed!

(-d $ENV{'browserDir'}) or confess "cannot open browser directory (which should already exist!) which was expected to be named this:  $ENV{'browserDir'}";
bananas_agw::systemAndLog("cat $ENV{'browserDir'}/*${trackSuff} > $tempTrackFile", $verbose);
bananas_agw::dieIfFileAccessFails($tempTrackFile);
bananas_agw::systemAndLog("/bin/rm $ENV{'browserDir'}/*${trackSuff}", $verbose); # remove all the useless input track files...
bananas_agw::systemAndLog("/bin/mv -f $tempTrackFile $finalSummaryTracksFile", $verbose);
bananas_agw::dieIfFileAccessFails($finalSummaryTracksFile);

#this section simply prints out the command to a file instead of copying the files directly, deprecated...
#my $scpFile = catfile($ENV{'browserDir'},"copyTracks.sh");
#open(OF,">$scpFile") or die "cannot open scp script file for writing\n";
#printf OF "cp $ENV{'browserDir'}/* $ENV{'scp'}\n";
#close(OF);

bananas_agw::mkdirOrDie($ENV{'rsyncDir'});
bananas_agw::systemAndLog("rsync --copy-links --chmod=u+r,g+r,o+r --perms $ENV{'browserDir'}/*   $ENV{'rsyncDir'}\n", $verbose);

