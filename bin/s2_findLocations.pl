#!/bin/perl
eval 'exec /bin/perl -w -S $0 ${1+"$@"}'
if 0; # not running under some shell

use strict;
use warnings;
use Data::Dumper;

unless (scalar @ARGV==2) { die "requires an input file of regexes from step1 and a fasta file with headers like >chr1:3002825-3003025" }
# program assumes that line breaks are not found within sequence stretches...

sub revComp {
    my ($s) = @_;
    my $rc = reverse($s);
    $rc =~ tr/ACTGactg/TGACtgac/;
    return ($rc);
}

my $m;
open(F0,$ARGV[0]) or die "Cannot open $ARGV[0]";
while(<F0>) {
    chomp;
    my @t = split(/\t/,$_);
    $m->{$t[0]}->{fwd} = $t[2];
    $m->{$t[0]}->{fwdseq} = $t[1];
    $m->{$t[0]}->{rev} = $t[3];
    $m->{$t[0]}->{revseq} = revComp($t[1]);
}
close(F0);

my $onChr = "";
my $start = 0;
my $stop = 0;
open(F0,$ARGV[1]) or die "Cannot open $ARGV[0]";
while(<F0>) {
    chomp;
    if ($_ =~ /^\>(\S{1,})\:(\d{1,})-(\d{1,})/) {
	$onChr = $1;
	$start = $2 - 1;
	$stop = $3;
    }
    elsif ($onChr ne "" && $_ =~ /[ACGTUacgtu]/) {
	foreach my $k (sort keys %$m) {
	    #printf STDERR "looking for: $m->{$k}->{fwd}";
	    while ($_ =~ /$m->{$k}->{fwd}/g) {
		my $tstart = $start + $-[0] + 1;
		my $tstop = $start + $+[0] + 1;
		printf "$onChr\t$tstart\t$tstop\t$k\t1\t\+\t$m->{$k}->{fwdseq}\n";
	    }
	    while ($_ =~ /$m->{$k}->{rev}/g) {
                my $tstart = $start + $-[0] + 1;
                my $tstop = $start + $+[0] + 1;
                printf "$onChr\t$tstart\t$tstop\t$k\t1\t\-\t$m->{$k}->{revseq}\n";
            }
	}
    }
}
close(F0);
