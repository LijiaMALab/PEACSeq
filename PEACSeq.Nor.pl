#!/usr/bin/perl
#Informatic Biology departments of Beijing Genomics Institute (BGI) 
use strict;
use Getopt::Long;
my %opts;

GetOptions(\%opts, "n:i","bait:s","help!");
if (@ARGV != 4) {
	print "usage: .pl ../293TAll.T3.OWG.NotI5.F.PEbwa.bam.log 293TAll.T3.OWG.NotI19.R.PEbwa.bam.log T1.NotI10F.NotI19R.PEACSeq out
		-n	1000000 (total reads number normalized to, by default 1M)
		-f	filetype
		\n";
	exit;
}

my $NorT=1000000;
$NorT=$opts{"n"} if($opts{'n'} ne "");

my $infile1 = shift;
$infile1 ="zcat $infile1|" if($infile1 =~/\.gz$/);
open(IN1, $infile1) or die $!;
my $infile2 = shift;
$infile2 ="zcat $infile2|" if($infile2 =~/\.gz$/);
my $infile3 = shift;
$infile3 ="zcat $infile3|" if($infile3 =~/\.gz$/);
open(IN3, $infile3) or die $!;
my $outfile = shift;
if($outfile =~/\.gz$/){
	$outfile ="|gzip >$outfile";
}else{
	$outfile =">$outfile";
}
open(OUT, "$outfile") or die $!;


my $NumReads=0;
while (<IN1>) {
	next if (/^#/);
	chomp;
	if(/read2/){
		if(/^(\d+) \+ 0 read2/){
			$NumReads += $1;
			print "ReadsNum:$NumReads\t$_\t$infile1\n";
		}else{
			print "Error, cant match ReadsNum $_\n";
		}
	}
}
close IN1;
while (<IN1>) {
	next if (/^#/);
	chomp;
	if(/read2/){
		if(/(\d+) \+ 0 read2/){
			$NumReads += $1;
			print "ReadsNum:$NumReads\t$_\n";
		}else{
			print "Error, cant match ReadsNum $_\n";
		}
	}
}
close IN1;

$"="\t";
my @col=(10,11,12,14,15,16,18,19,20,21,22,23,37,38,39,40,42,43,44,45);
my $Flag1=0;
my $Flag2=0;
while (<IN3>) {
	chomp;
	my @info = split/\t/;
	if (/^#/){
		$info[47]="Flag2";
		print OUT "@info\n";
		next;
	}
	if(($info[32]eq"+"&&($info[37]>1||$info[43]>1))||($info[32]eq"-"&&($info[38]>1||$info[42]>1))){
		$Flag1=1;
		$info[46]=$Flag1;
#		if($Flag1 != $info[46]){
#			print "Error, $_ not match\n";
#		}
	}
	if(($info[32]eq"+"&&$info[37]>1&&$info[43]>1)||($info[32]eq"-"&&$info[38]>1&&$info[42]>1)){
		$Flag2=1;
		$info[47]=$Flag2;
	}
	foreach my$c(@col){
		$info[$c-1] = sprintf("%.6f",$info[$c-1]*$NorT/$NumReads);
	}
	print OUT "@info\n";
}
close IN3;
close OUT;
