#!/usr/bin/perl -w

# https://metacpan.org/pod/distribution/Test-Simple/lib/Test/Tutorial.pod

use Test::More tests => 3;
use File::Temp;
use Carp;


my $SHOULD_GENERATE_GOLD_STANDARD_FILES = 0;
my $GOLD_PARENT_DIR = "/data/applications/monkey/test_gold";

my $TESTDIR = "/data/applications/monkey/test_gold"; #File::Temp->newdir(CLEANUP=>0); # CLEANUP=>0 so we don't auto-delete delete the temp directory

#`mkdir -p $TESTDIR`;

#print qq{Testing in this test dir: "$TESTDIR" ...};


#chdir($TESTDIR) or confess qq{Failed to change to test directory.};

my $forceGold = 0;
my $goldInBackground = 1;

sub notify($) {
	my ($msg) = @_;
	chomp($msg);
	print STDERR $msg . "\n";
}

sub make_gold_standard($$$;$$) {
	my ($subfolderName, $goldCmd, $forceRerun, $runInBackground) = @_;
	if (!defined($forceRerun)) { $forceRerun = 0; }
	if (!defined($runInBackground)) { $runInBackground = 0; }

	my $fullPath = "${GOLD_PARENT_DIR}/$subfolderName";
	my $expectedOutputFile = "${fullPath}/99.gold.standard.ok.touch";
	if (-e $expectedOutputFile and !$forceRerun) {
		notify(qq{[OK] [SKIPPING] Output file '$expectedOutputFile' already exists, so not re-running this command: $goldCmd});
		return;
	}
	if ($forceRerun) {
		warning("Note: force rerun is NOT currently deleting the existing files, which it probably should.");
		warning("Probably should be deleting the existing output files here.");
	}
	my $outDirOption = " --out=$fullPath ";

	my $completeCmd = "$goldCmd $outDirOption";
	(system($goldCmd) == 0) or confess "[ERROR] Gold-standard-making in '$subfolderName' failed with non-zero exit code! Command was: $goldCmd";
	(system(qq{touch "$expectedOutputFile"})) or confess qq{[ERROR] Couldn't make the 'we have finished successfully' file ($expectedOutputFile).};
}

make_gold_standard("test01_bowtie", "monkey_agw --test1 --no-scheduler", $forceGold, $goldInBackground);




#ok( 1 + 1 == 2 );
#ok( 1 + 1 == 3 );
#ok( 2 + 2 == 4 , " yeah math works");
#is( 2, 2, " yep sure is");
