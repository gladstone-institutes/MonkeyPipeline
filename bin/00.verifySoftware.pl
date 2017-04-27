#!/usr/bin/perl -w

use strict; use warnings; use Carp; # Carp = "confess" function
use File::Spec::Functions qw(catfile);   use Data::Dumper;   use File::Basename;
use lib "/data/work/Code/alexgw/monkey_agw"; use bananas_agw; # To make sure we find bananas_agw

$ENV{'req_exes'} = "bash perl python fsadf";
$ENV{'req_r_libs'} = "edgeR";

$ENV{'req_python2_libs'} = "matplotlib requests Tkinter" # <-- these 3 are for AltAnalyze Lineage Profiler. But we can mostly live without them (?)
  . "";

$ENV{'req_perl_libs'} = "XML";


# ALTANALYZE REQUIRES that the "requests" package be installed, and also numpy
#      sudo pip install --upgrade requests
#      sudo pip install --upgrade numpy

#bananas_agw::requireEnvOrDie('req_exes','req_r_libs','req_python2_libs','req_perl_libs');
my ($force, $verbose, $debug) = bananas_agw::envLooksTrue("force", "verbose", "debug");

# Checks everything at "exe_paths" to make sure they are executable... maybe
# Delim is a ":"

my @req_exes         = split(/\s+/, $ENV{'req_exes'});
my @req_r_libs       = split(/\s+/, $ENV{'req_r_libs'});
my @req_python2_libs = split(/\s+/, $ENV{'req_python2_libs'});
my @req_perl_libs    = split(/\s+/, $ENV{'req_perl_libs'});

sub explode($) {
	my ($msg) = @_;
	confess "[ERROR] Cannot find a required component on the COMPUTE SERVER!\nSpecifically:\n      $msg\nNote that this is a problem on the COMPUTE SERVER, not the login node! Check your permissions on the COMPUTE SERVER and make sure this file/program/library is present and correct.\nAlso make sure that any executable programs give the FULL ABSOLTUE PATH, not just a program name. Again: this is A PROBLEM ON THE COMPUTE SERVER, so do not try to debug it on the login node! Also be careful about different versions of programs (i.e. there may be more than one installation of 'python' or 'wget' or other program).";
}

foreach my $x (@req_exes) {
	my $exitCode = system(qq{which '$x' > /dev/null});
	if (0 == $exitCode) {
		# looks like we found the required exectuable! No problem, let's go on to the next one.
	} else {
		$x =~ m/^[\/]/ or warn("The required exectuable <$x> seemed to not start with a '/', and also we could not find it using 'which'. This probably means you did not give the full path to it. This probably-incorrect (or absent) path will almost certainly make it impossible for us to properly run that executable.");
		-e $x or explode("The required exectuable <$x> could not be found!");
		-x $x or explode("The required exectuable <$x> was found, but was not EXECUTABLE. Someone will need to 'chmod +x' it.");
	}
 	# Check to make sure each of these exist!
}

foreach my $r_lib (@req_r_libs) {
	# require('...')
}

foreach my $py_lib (@req_python2_libs) {
	bananas_agw::dieIfMissingPythonLib($pythonExe, $py_lib);
}
#foreach my $recommended_py (@recommended_python2_libs) {
#	if (not bananas_agw::pythonHasLib($pythonExe, $py_lib)) {
#		my $msg = "[WARNING]: Continuing... not a fatal problem, but something is going to probably break because we are missing the python library '$recommended_py' on the computer that is actually running the computation. Check 1) your version of python and 2) that the PYTHONPATH is set correctly!";
#		print STDERR $msg;
#		print STDOUT $msg;
#		warning($msg);
#	}
#}

foreach my $perl_lib (@req_perl_libs) {
	# use ...;
}

# R things to look for:
# edgeR
# DESeq2
# gapmap
# genefilter

# cummeRbund
#data.table
#pheatmap
#ggplot2
#S4Vectors
#methods
#stats4
#BiocGenerics
#parallel
