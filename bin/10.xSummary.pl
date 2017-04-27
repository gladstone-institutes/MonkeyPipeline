#!/usr/bin/perl





use strict; use warnings; use Carp; # Carp has "confess"
use File::Spec::Functions qw(catfile);
use Data::Dumper;
use lib "/data/work/Code/alexgw/monkey_agw"; use bananas_agw; # To make sure we find bananas_agw

bananas_agw::requireEnvOrDie('studyName','motifDir','peakDir','bedops','bedmap','bedtools','genomeFasta','step1','step2','sortBed'); # <-- these vars MUST be defined in the %ENV hash
my $force       = bananas_agw::envLooksTrue("force");
my $verbose     = bananas_agw::envLooksTrue("verbose");
my $allMotifs = catfile($ENV{'motifDir'},"0.$ENV{'studyName'}_allDiscoveredMotifs.tf");
my $allRegex  = catfile($ENV{'motifDir'},"0.$ENV{'studyName'}_allDiscoveredMotifs.regex");
my $profile   = catfile($ENV{'motifDir'},"0.$ENV{'studyName'}_allDiscoveredMotifs.prf");

(-d $ENV{'motifDir'}) or die "cannot open motif discovery directory $ENV{'motifDir'}";
(-d $ENV{'peakDir'}) or die "cannot open motif discovery directory $ENV{'peakDir'}";

if ($force || ! -e $allMotifs || ! -e $profile || ! -e $allRegex) {
    # get samples with motifs
    my $samples;
    opendir(INDIR,$ENV{'motifDir'}) or die "cannot open motif discovery directory: $ENV{'motifDir'}";
    while (my $sdir = readdir(INDIR)) {
      if ($sdir =~ /^(.*?\.\d{1,})\_pmResults$/) {
	my $sName = $1;
	my $sBase = catfile($ENV{'motifDir'},$sdir,"results","discovered_motifs","${sName}_motifs_discovered");
	if (-e "${sBase}.tf" && -e "${sBase}_table.tab") {
	    $samples->{$sName} = $sBase;
	    #printf STDERR "$sName\t$samples->{$sName}\n";
	}
      }
    }
    closedir(INDIR);

    # for each sample, load table, substitute motif names, and print summary file
    open(AMO,">$allMotifs") or die "Cannot open $allMotifs";
    printf AMO "VV  Custom TRANSFAC MATRIX TABLE\nXX\n\/\/\n";
    close(AMO);

    open(OF,">$profile") or die "Cannot open $profile";
    printf OF "$ENV{'studyName'} motifs\n0.$ENV{'studyName'}_allDiscoveredMotifs.prf\n MIN_LENGTH  300\n0.0\n";
    foreach my $s (sort keys %$samples) {
      my $table   = $samples->{$s} . "_table.tab";
      my $tf      = $samples->{$s} . ".tf";
      my $tmpFile = "${tf}.tmp";
      bananas_agw::systemAndLog("cp $tf $tmpFile", $verbose);
      if (-e $table && -e $tf) {
	my $fr=1;
	#printf STDERR "$table";
	open(F0,$table) or die "Cannot open $table";
	while(<F0>) {
	    chomp;
	    if ($fr) { $fr=0 }
	    elsif ($_ =~ /NO MOTIF FOUND/) { }
	    else {
		my @t = split(/\t/,$_);
		#printf STDERR "$_\n";
		my $oldName = $t[1];
		my $newName = $s . "_motif_" . $t[0];
		bananas_agw::systemAndLog("sed -i 's/$oldName/$newName/' $tmpFile", $verbose);
		printf OF " 1.000 0.700 0.980 $newName $newName\n";
	    }
	}
	close(F0);
	bananas_agw::systemAndLog("cat $tmpFile >> $allMotifs", $verbose);
	unlink($tmpFile);
      }
    }
    printf OF "\/\/\n"; # //[\n]
    close(OF);

    my $tfTmp = $allMotifs . ".tmp";
    bananas_agw::systemAndLog("$ENV{'step1'} $allMotifs > $allRegex", $verbose);
}

# get peak fasta sequences
my $pkBed = catfile($ENV{'motifDir'},"1.$ENV{'studyName'}_allPeaksMerged.bed");
my $pkBedTmp = catfile($ENV{'motifDir'},"1.$ENV{'studyName'}_allPeaksMerged.bed.tmp");
my $pkFa  = catfile($ENV{'motifDir'},"1.$ENV{'studyName'}_allPeaksMerged.fasta");
if ($force || ! -e $pkBed || ! -e $pkFa) {
    bananas_agw::systemAndLog("$ENV{'bedops'} -m $ENV{'peakDir'}/*peaks.bed > $pkBedTmp", $verbose);
    my $on = 1;
    open(F0,$pkBedTmp) or die;
    open(OF,">$pkBed") or die;
    while(<F0>) {
	chomp;
	my @t = split(/\t/,$_);
	if (scalar @t == 3) {
	    printf OF "$t[0]\t$t[1]\t$t[2]\tpk_$on\n";
	    $on++;
	}
    }
    close(F0);
    close(OF);
    unlink($pkBedTmp);
    bananas_agw::systemAndLog("$ENV{'bedtools'} getfasta -fi $ENV{'genomeFasta'} -bed $pkBed  >  $pkFa", $verbose); # previously was "-fo OUTPUT" --- now it goes to STDOUT as of version 2.25
}

# perform motif scan
my $motifMoa = catfile($ENV{'motifDir'},"2.$ENV{'studyName'}_allMotifs.moa");
if ($force || ! -e $motifMoa) {
    bananas_agw::systemAndLog("$ENV{'step2'} $allRegex $pkFa | $ENV{'sortBed'} - > $motifMoa", $verbose);
}

# get peak level summary
my $pkSummary = catfile($ENV{'motifDir'},"3.$ENV{'studyName'}_peakSummary.txt");
my $pkoMatrix = catfile($ENV{'motifDir'},"3.$ENV{'studyName'}_pkOverlap.mtx");
my $fpoMatrix = catfile($ENV{'motifDir'},"3.$ENV{'studyName'}_fpOverlap.mtx");
my $pkSamples;
my $fpSamples;
if ( $force || ! -e $pkSummary) {
    # get peak and footprint file names
    if ($force || ! -e $pkoMatrix || ! -e $fpoMatrix) {
	opendir(INDIR,$ENV{'peakDir'}) or die "cannot open peak directory: $ENV{'peakDir'}";
	while (my $file = readdir(INDIR)) {
	    if ($file =~ /^(.*?\.\d{1,})\_.*peaks.bed$/) {
		my $sName = $1;
		$pkSamples->{$sName} = $file;
		my $opFile = catfile($ENV{'motifDir'},"${sName}_pkOverlaps.map");
		bananas_agw::systemAndLog("$ENV{'bedmap'} --count $pkBed $ENV{'peakDir'}/$file > $opFile", $verbose);
	    }
	    elsif ($file =~ /^(.*?\.\d{1,})\_.*footprints.bed$/) {
		my $sName = $1;
		$fpSamples->{$sName} = $file;
		my $opFile = catfile($ENV{'motifDir'},"${sName}_fpOverlaps.map");
		bananas_agw::systemAndLog("$ENV{'bedmap'} --count $pkBed $ENV{'peakDir'}/$file > $opFile", $verbose);
	    }
	}
	closedir(INDIR);
    }
    
    # generate pk matrix
    if ($force || ! -e $pkoMatrix) {
	my $pkoMatrixTmp = catfile($ENV{'motifDir'},"3.$ENV{'studyName'}_pkOverlap.mtx.tmp");
	my $pkoHead = catfile($ENV{'motifDir'},"3.$ENV{'studyName'}_pkOverlap.mtx.head");
	bananas_agw::systemAndLog("ls -1 $ENV{'motifDir'}/*pkOverlaps.map | awk '{gsub(/\\/.*\\//,\"\",\$1); print}' | sed 's/pkOverlaps.map/pkO/' | tr '\\n' '\\t' | sed 's/\\t\$//' | awk '{print \$0}' > $pkoHead", $verbose);
	bananas_agw::systemAndLog("paste $ENV{'motifDir'}/*pkOverlaps.map > $pkoMatrixTmp", $verbose);
	bananas_agw::systemAndLog("cat $pkoHead $pkoMatrixTmp > $pkoMatrix", $verbose);
	unlink($pkoMatrixTmp);
	unlink($pkoHead);
    }

    # generate fp matrix
    if ($force || ! -e $fpoMatrix) {
	my $fpoMatrixTmp = catfile($ENV{'motifDir'},"3.$ENV{'studyName'}_fpOverlap.mtx.tmp");
	my $fpoHead = catfile($ENV{'motifDir'},"3.$ENV{'studyName'}_fpOverlap.mtx.head");
	bananas_agw::systemAndLog("ls -1 $ENV{'motifDir'}/*fpOverlaps.map | awk '{gsub(/\\/.*\\//,\"\",\$1); print}' | sed 's/fpOverlaps.map/fpO/' | tr '\\n' '\\t' | sed 's/\\t\$//' | awk '{print \$0}' > $fpoHead", $verbose);
	bananas_agw::systemAndLog("paste $ENV{'motifDir'}/*fpOverlaps.map > $fpoMatrixTmp", $verbose);
	bananas_agw::systemAndLog("cat $fpoHead $fpoMatrixTmp > $fpoMatrix", $verbose);
	unlink($fpoMatrixTmp);
	unlink($fpoHead);
    }

    # generate motif matrix
    my $motifOverlaps = catfile($ENV{'motifDir'},"3.$ENV{'studyName'}_motifOverlaps.map");
    if ($force || ! -e $motifOverlaps) {
	bananas_agw::systemAndLog("$ENV{'bedmap'} --echo-map-id $pkBed $motifMoa > $motifOverlaps", $verbose);
    }

    # load motifs
    my $mtH;
    my $mtA;
    my $on = 0;
    my $motifMtx = catfile($ENV{'motifDir'},"3.$ENV{'studyName'}_motif.mtx");
    open(F0,$allRegex) or die "cannot open motif regex file $allRegex";
    open(OF,">$motifMtx") or die "cannot open $motifMtx for writing";
    while(<F0>) {
	chomp;
	my @t = split(/\t/,$_);
	$mtH->{$t[0]}->{seq} = $t[1];
	$mtH->{$t[0]}->{on}  = $on;
	$mtA->[$on] = $t[0];
	if ($on>0) {printf OF "\t"}
	printf OF "$t[1]";
	$on++
    }
    close(F0);
    printf OF "\n";

    # load map file and 
    open(F0,$motifOverlaps) or die "cannot open motif overlaps file $motifOverlaps";
    while(<F0>) {
        chomp;
	my @t = split(/\;/,$_);
	my $d;
	foreach my $k (@t) {
	    $d->{$k} = 1;
	}
	for (my $i=0; $i<$on; $i++) {
	    if ($i>0) {printf OF "\t";}
	    my $onmo = $mtA->[$i];
	    if (defined($d->{$onmo})) { printf OF "1"}
	    else { printf OF "0" }
	}
	printf OF "\n";
    }
    close(F0);
    close(OF);

    # assemble full matrix
    my $tmpPk   = catfile($ENV{'motifDir'},"3.$ENV{'studyName'}_pkFile.tmp");
    my $tmpHead = catfile($ENV{'motifDir'},"3.$ENV{'studyName'}_pkFile.head");
    bananas_agw::systemAndLog("echo \"chr\tstart\tstop\tID\"                > $tmpHead", $verbose);
    bananas_agw::systemAndLog("cat   $tmpHead $pkBed                        > $tmpPk", $verbose);
    bananas_agw::systemAndLog("paste $tmpPk $pkoMatrix $fpoMatrix $motifMtx > $pkSummary", $verbose);
    bananas_agw::systemAndLog("/bin/rm $ENV{'motifDir'}/*.map", $verbose);
    unlink($tmpPk);
    unlink($tmpHead);
    unlink($pkoMatrix);
    unlink($fpoMatrix);
    unlink($motifMtx);
}
