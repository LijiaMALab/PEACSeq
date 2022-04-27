#!/usr/bin/perl
#Informatic Biology departments of Beijing Genomics Institute (BGI) 
use strict;
use Getopt::Long;
my %opts;

GetOptions(\%opts, "f:s","r:s","q:i","b!","help!");
if (@ARGV != 2) {
	print "usage: .pl T3.15F.18R.PEACSeq.nor out
		-q	min mapq qulity, by defalut 0
		-b	quialiti is apply to both strands, by default not setted.
		-l	minimal length,by default 30nt
		-f	T6.15F.nodup.sort.bam
		-r	T6.18R.nodup.sort.bam
		\n";
	exit;
}

my $minQ=0;
$minQ=$opts{"q"} if($opts{'q'} ne "");
my $FlagBoth=0;
$FlagBoth=1 if(exists($opts{'b'}));
my $F="";
$F=$opts{"f"} if($opts{'f'} ne "");
$F ="samtools view $F|" if($F =~/\.bam$/);

my $infile2 = shift;
$infile2 ="zcat $infile2|" if($infile2 =~/\.gz$/);
open(IN2, $infile2) or die $!;
my $outfile = shift;
if($outfile =~/\.gz$/){
	$outfile ="|gzip >$outfile";
}else{
	$outfile =">$outfile";
}
open(OUT, "$outfile") or die $!;

my %Pos;
my %Sites;
while (<IN2>) {
	next if (/^#/);
	chomp;
	my @info = split/\t/;
	my $strand = $info[32];
	my $chr = $info[27];
	$chr =~s/^chr//i;
	my $p = $info[29]-3;
	$p = $info[28]+3 if($strand eq "-");
	for(my $i=($p-30);$i<($p-10);$i++){
		$Pos{$chr}{$i}{'U'} .="$info[30]\t";
	}
	for(my $i=($p+10);$i<($p+30);$i++){
		$Pos{$chr}{$i}{'D'} .="$info[30]\t";
	}
	if($Sites{$info[30]} ne ""){
		print "Error, $info[30] is not uniq\n";
	}
	$Sites{$info[30]}{'detail'} = $_;
	$Sites{$info[30]}{'strand'} = $strand;
	$Sites{$info[30]}{'chr'} = $chr;
}
close IN2;
open(IN2, $infile2) or die $!;

open(IN1F, $F) or die $!;
my $c=0;
my %FlagLQ;
while (<IN1F>) {
	last if($FlagBoth==0);
	next if (/^#/);
	chomp;
	my @info = split/\t/;
	if($info[4]<$minQ){
		$FlagLQ{$info[0]}=1;
	}
	if($c%1000000==0){
		print "reads $c reads, $F\n";
	}
	$c++;
}
close IN1F;
##################
open(IN1F, $F) or die $!;
my $c=0;
while (<IN1F>) {
	next if (/^#/);
	chomp;
	my @info = split/\t/;
	my $chr = $info[2];
	my $strand;
	next if($info[4]<$minQ);
	next if($FlagBoth==1 && $FlagLQ{$info[0]}==1);
	my $matchLen=0;
	my $TotalLen = length($info[9]);
	while($info[5]=~/(\d+)M/g){
		$matchLen+=$1;
	}
	next if($matchLen<30 && $matchLen/$TotalLen<0.5);
	$strand = "+" if(($info[1]&32)==32);
	if(($info[1]&128)==128){
		$strand = "-" if(($info[1]&16)==16);
		my ($s,$e);
		if($strand eq "+"){
			$s = $info[3];
			$e = $info[3]+$info[8];
		}elsif($strand eq "-"){
			$s = $info[7];
			$e = $info[7]+abs($info[8]);
		}
		for(my$i=$s;$i<=$e;$i++){
			if($Pos{$chr}{$i}{'U'} ne ""){
				my @S = split/\t/,$Pos{$chr}{$i}{'U'};
				foreach my$s(@S){
					$Sites{$s}{'FU'}{$strand}{$info[0]}=1;
				}
			}
			if($Pos{$chr}{$i}{'D'} ne ""){
				my @S = split/\t/,$Pos{$chr}{$i}{'D'};
				foreach my$s(@S){
					$Sites{$s}{'FD'}{$strand}{$info[0]}=1;
				}
			}
		}
	}
	if($c%1000000==0){
		print "reads $c reads, $F\n";
	}
	$c++;
}
close IN1F;
undef %FlagLQ;

my $R="";
$R=$opts{"r"} if($opts{'r'} ne "");
$R ="samtools view $R|" if($R =~/\.bam$/);
if($opts{'r'} ne ""){
	open(IN1R, $R) or die $!;
	my $c=0;
	my %FlagRQ;
	while (<IN1R>) {
		last if($FlagBoth==0);
		next if (/^#/);
		chomp;
		my @info = split/\t/;
		if($info[4]<$minQ){
			$FlagRQ{$info[0]}=1;
		}
		if($c%1000000==0){
			print "reads $c reads, $R\n";
		}
		$c++;
	}
	close IN1R;
###############
	open(IN1R, $R) or die $!;
	my $c=0;
	while (<IN1R>) {
		next if (/^#/);
		chomp;
		my @info = split/\t/;
		next if($info[4]<$minQ);
		next if($FlagBoth==1 && $FlagRQ{$info[0]}==1);
		my $matchLen=0;
		my $TotalLen = length($info[9]);
		while($info[5]=~/(\d+)M/g){
			$matchLen+=$1;
		}
		next if($matchLen<30 && $matchLen/$TotalLen<0.5);
		my $chr = $info[2];
		my $strand;
		$strand = "+" if(($info[1]&32)==32);
		if(($info[1]&128)==128){
			$strand = "-" if(($info[1]&16)==16);
			my ($s,$e);
			if($strand eq "+"){
				$s = $info[3];
				$e = $info[3]+$info[8];
			}elsif($strand eq "-"){
				$s = $info[7];
				$e = $info[7]+abs($info[8]);
			}
			for(my$i=$s;$i<=$e;$i++){
				if($Pos{$chr}{$i}{'U'} ne ""){
					my @S = split/\t/,$Pos{$chr}{$i}{'U'};
					foreach my$s(@S){
						$Sites{$s}{'RU'}{$strand}{$info[0]}=1;
					}
				}
				if($Pos{$chr}{$i}{'D'} ne ""){
					my @S = split/\t/,$Pos{$chr}{$i}{'D'};
					foreach my$s(@S){
						$Sites{$s}{'RD'}{$strand}{$info[0]}=1;
					}
				}
			}
		}
		if($c%1000000==0){
			print "reads $c reads, $R\n";
		}
		$c++;
	}
	close IN1R;
	undef %FlagRQ;
}
undef %Pos;

while (<IN2>) {
	chomp;
	$_=~s/\t$//;
	if (/^#/){
		print OUT "$_\tFUPnum\tFDPnum\tFUNum\tFDNnum\tFDSBRatio\tRUPnum\tRDPnum\tRUNnum\tRDNnum\tRDSBRatio\tFlag1\tFlag10\n";
		next;
	}
	my @info = split/\t/;
	my ($FUPnum,$FUNnum,$RUPnum,$RUNnum,$FDPnum,$FDNnum,$RDPnum,$RDNnum);
	$FUPnum = scalar(keys %{$Sites{$info[30]}{'FU'}{"+"}});
	$FDPnum = scalar(keys %{$Sites{$info[30]}{'FD'}{"+"}});
	$FUNnum = scalar(keys %{$Sites{$info[30]}{'FU'}{"-"}});
	$FDNnum = scalar(keys %{$Sites{$info[30]}{'FD'}{"-"}});
	my ($FRatio,$RRatio);
	$FRatio = sprintf("%.6f",($FUNnum+0)/($FDPnum+1+$FUNnum)) if($Sites{$info[30]}{'strand'} eq "+");
	$FRatio = sprintf("%.6f",($FDPnum+0)/($FUNnum+1+$FDPnum)) if($Sites{$info[30]}{'strand'} eq "-");
	$RUPnum = scalar(keys %{$Sites{$info[30]}{'RU'}{"+"}});
	$RDPnum = scalar(keys %{$Sites{$info[30]}{'RD'}{"+"}});
	$RUNnum = scalar(keys %{$Sites{$info[30]}{'RU'}{"-"}});
	$RDNnum = scalar(keys %{$Sites{$info[30]}{'RD'}{"-"}});
	$RRatio = sprintf("%.6f",($RDPnum+0)/($RUNnum+1+$RDPnum)) if($Sites{$info[30]}{'strand'} eq "+");
	$RRatio = sprintf("%.6f",($RUNnum+0)/($RDPnum+1+$RUNnum)) if($Sites{$info[30]}{'strand'} eq "-");
	my $Flag1=0;
	my $Flag10=0;
	if(($Sites{$info[30]}{'strand'}eq"+" && $FDPnum>0 && $RUNnum>0)||($Sites{$info[30]}{'strand'}eq"-" && $FUNnum>0 && $RDPnum>0)){
		$Flag1=1;
	}
	if(($Sites{$info[30]}{'strand'}eq"+" && $FDPnum>=10 && $RUNnum>=10)||($Sites{$info[30]}{'strand'}eq"-" && $FUNnum>=10 && $RDPnum>=10)){
		$Flag10=10;
	}
	print OUT "$_\t$FUPnum\t$FDPnum\t$FUNnum\t$FDNnum\t$FRatio\t$RUPnum\t$RDPnum\t$RUNnum\t$RDNnum\t$RRatio\t$Flag1\t$Flag10\n";
}
close IN2;
close OUT;
