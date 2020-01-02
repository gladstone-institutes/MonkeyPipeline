#!/bin/perl
use strict;
use warnings;
use Data::Dumper;

sub checkCase {
    if ($_[0] =~ /^[[:upper:]]/) {
        return 1;
    }
    else {
        return 0;
    }
}

unless (scalar @ARGV==1) {die "requires an input transfac matrix file output from RSAT peak-motif"}
my $m;
my $on = "";
open(F0,$ARGV[0]) or die "cannot open $ARGV[0]";
while (<F0>) {
    chomp;
    if ($_ =~ /^ID\s{1,}(.*)$/) {
	$on = $1;
    }
    elsif ($_ =~ /^DE\s{1,}(.*)$/) {
	my $iup = $1;
	my $p1 = 0;
	my $p2 = 0;
	
	for (my $i=0; $i<length($iup); $i++) {
	    my $t = substr($iup,$i,1);
	    if (checkCase($t)) {
		if ($p1 == 0) { $p1 = $i }
		$p2 = $i;
	    }
	}
	$m->{$on}->{iupac} = $iup;
	my $me = $p2 - $p1 + 1;
	$m->{$on}->{core} = substr($iup,$p1,$me);
	$m->{$on}->{p1} = $p1;
	$m->{$on}->{p2} = $p2;
    }
    elsif ($_ =~ /^CC\s{1,}consensus.regexp:\s{1,}(.*)$/) {
	my $reg = $1;
	my $u = "";
	my $building = 0;
	my $on2 = 0;
	for (my $i=0; $i<length($reg); $i++) {
	    my $t = substr($reg,$i,1);
	    if ($t eq "\[") {
		$building = 1;
		$u = "\[";
	    }
	    elsif ($t eq "\]") {
		$building = 0;
		$m->{$on}->{fwd}->[$on2] = $u . "\]";
		$u = "";
		$on2++;
	    }
	    elsif ($building) {
		$u .= $t;
	    }
	    else {
		$m->{$on}->{fwd}->[$on2] = $t;
		$on2++;
	    }
	}
    }
    elsif ($_ =~ /^CC\s{1,}consensus.regexp.rc:\s{1,}(.*)$/) {
	my $reg = $1;
        my $u = "";
        my $building = 0;
        my $on2 = 0;
        for (my $i=0; $i<length($reg); $i++) {
            my $t = substr($reg,$i,1);
            if ($t eq "\[") {
                $building = 1;
		$u = "\[";
            }
            elsif ($t eq "\]") {
                $building = 0;
                $m->{$on}->{rev}->[$on2] = $u . "\]";
                $u = "";
                $on2++;
            }
            elsif ($building) {
                $u .= $t;
            }
            else {
                $m->{$on}->{rev}->[$on2] = $t;
                $on2++;
            }
        }
    }
}
close(F0);

foreach my $k (sort keys %$m) {
    my $s = $m->{$k}->{iupac};
    my $fwd = "";
    my $rev = "";
    my $p1r = length($s) - 1 - $m->{$k}->{p2};
    my $p2r = length($s) - 1 - $m->{$k}->{p1};
    
    for (my $i=0; $i<length($s); $i++) {
	if ($i >= $m->{$k}->{p1} && $i <= $m->{$k}->{p2}) {
	    $fwd = $fwd . $m->{$k}->{fwd}->[$i];
	}
	if ($i >= $p1r && $i <= $p2r) {
            $rev = $rev . $m->{$k}->{rev}->[$i];
        }
    }

    printf "$k\t$m->{$k}->{core}\t$fwd\t$rev\n";
}
