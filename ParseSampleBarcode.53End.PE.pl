#!/usr/bin/perl
use strict;
use Getopt::Long;
my %opts;

GetOptions(\%opts, "L5:s","R5:s","L3:s","R3:s","v1!","V2!","BR5fq!","gz!","N!","help!");
if (@ARGV != 4) {
	print "usage: .pl file1-1[,file2-1](Reads1) file1-2[,file2-2](Reads2) out1 out2
		-L5	Read1 UMI length and suffix, by default(0,0,)
		-L3	Read1 UMI length and prefix, by default(0,0,)
		-R5	Read2 UMI length and suffix, by default(0,0,)
		-R3	Read2 UMI length and prefix, by default(0,0,)
		-v1	by default not set, by default discard reads without barcode,
			setted: output reads without barcode.
		-V2	by default not set.
			setted: only report barcode detail,discard all reads.
		-BR5fq	by default not set.
			setted: report barcode R5 in fastq format.
		-gz	bydeault not set. input is gzip or not.
		-N	bydeault not set. N is allowed or not.
		\n";
	exit;
}

my ($L5Lmin,$L5Lmax,$L5suffix)=(0,0,"");
($L5Lmin,$L5Lmax,$L5suffix)=split/,/,$opts{"L5"} if($opts{'L5'} ne "");
my ($R5Lmin,$R5Lmax,$R5suffix)=(0,0,"");
($R5Lmin,$R5Lmax,$R5suffix)=split/,/,$opts{"R5"} if($opts{'R5'} ne "");
my ($L3Lmin,$L3Lmax,$L3prefix)=(0,0,"");
($L3Lmin,$L3Lmax,$L3prefix)=split/,/,$opts{"L3"} if($opts{'L3'} ne "");
my ($R3Lmin,$R3Lmax,$R3prefix)=(0,0,"");
($R3Lmin,$R3Lmax,$R3prefix)=split/,/,$opts{"R3"} if($opts{'R3'} ne "");
$L5suffix = uc($L5suffix);
$R5suffix = uc($R5suffix);
$L3prefix = uc($L3prefix);
$R3prefix = uc($R3prefix);
print "L5:$L5Lmin,$L5Lmax,$L5suffix\n";
print "R5:$R5Lmin,$R5Lmax,$R5suffix\n";
print "L3:$L3Lmin,$L3Lmax,$L3prefix\n";
print "R3:$R3Lmin,$R3Lmax,$R3prefix\n";

my $FlagR5fq=0;
$FlagR5fq=1 if(exists $opts{'BR5fq'});
my $Flagv=0;
$Flagv=1 if(exists $opts{'v1'});
my $FlagVall=0;
if(exists $opts{'V2'}){
	$FlagVall=1;
	$Flagv=0;
}
my $CharN="ATGC";
$CharN="ATGCN" if(defined $opts{'N'});
print "v:$Flagv; V:$FlagVall; N:$CharN\n";
my @infile1 = split/,/,shift;
my @infile2 = split/,/,shift;

# my $key = shift;
my $out_fwd = shift;
my $out_rev = shift;
my $key = $out_fwd;
$key =~ s/.1.fastq.gz$//g;
open(OUT, "|gzip >$key.barcode.gz") or die $!;
open(OUT1, "|gzip >$key.barcode.log.gz") or die $!;
#############output in all, L5,R5,L3,R3 order
open(OUT2, ">$key.barcodeUMI.saturation") or die $!;
open(OUT5, ">$key.Fragment.saturation") or die $!;

open(OUT3, "|gzip>$out_fwd") or die $!;
open(OUT4, "|gzip>$out_rev") or die $!;
open(OUT6, "|gzip>$key.nobarcode.1.fastq.gz") or die $!;
open(OUT7, "|gzip>$key.nobarcode.2.fastq.gz") or die $!;

open(OUT8, "|gzip>$key.BarcodeR5.fastq.gz") if($FlagR5fq==1);

my %Barcode;
my %UniqFrag;
my %UniqUMI;
my $Line=0;
my $BarcodeNum=0;
for(my$i=0;$i<@infile1;$i++){
	if(defined $opts{'gz'}){
		$infile1[$i] = "zcat $infile1[$i]|";
		$infile2[$i] = "zcat $infile2[$i]|";
		print "input file @infile1 @infile2 if gziped\n";
	}
	open(IN1, $infile1[$i]) or die $!;
	open(IN2, $infile2[$i]) or die $!;
	while (<IN1>) {
		chomp;
		my $L11 = $_;
		my $L12 = <IN1>;
		my $L13 = <IN1>;
		my $L14 = <IN1>;
		my $L21 = <IN2>;
		my $L22 = <IN2>;
		my $L23 = <IN2>;
		my $L24 = <IN2>;
		chomp $L11;
		chomp $L12;
		chomp $L13;
		chomp $L14;
		chomp $L21;
		chomp $L22;
		chomp $L23;
		chomp $L24;
		my ($L5Barcode,$R5Barcode,$L3Barcode,$R3Barcode,$LRead,$RRead,$Barcode);
		if($L5Lmax>0 && $L12=~/^([$CharN]{$L5Lmin,$L5Lmax})$L5suffix/){
			$L5Barcode = $1;
			$Barcode .="$L5Barcode:";
		}elsif($L5Lmax>0){
			print OUT1 "Error, wrong L-5Barcode $L11\t$L12\n";
			if($Flagv!=0){
				print OUT6 "$L11\n$L12\n$L13\n$L14\n";
				print OUT7 "$L21\n$L22\n$L23\n$L24\n";
			}
			next;
		}
		if($R5Lmax>0 && $L22=~/^([$CharN]{$R5Lmin,$R5Lmax})$R5suffix/){
			$R5Barcode = $1;
			$Barcode .="$R5Barcode:";
		}elsif($R5Lmax>0){
			print OUT1 "Error, wrong R-5Barcode $L21\t$L22\n";
			if($Flagv!=0){
				print OUT6 "$L11\n$L12\n$L13\n$L14\n";
				print OUT7 "$L21\n$L22\n$L23\n$L24\n";
			}
			next;
		}
		if($L3Lmax>0 && $L12=~/$L3prefix([$CharN]{$L3Lmin,$L3Lmax})$/){
			$L3Barcode = $1;
			$Barcode .="$L3Barcode:";
#		print "L3:$1\n";
		}elsif($L3Lmax>0){
			print OUT1 "Error, wrong L-3Barcode $L11\t$L12\n";
			if($Flagv!=0){
				print OUT6 "$L11\n$L12\n$L13\n$L14\n";
				print OUT7 "$L21\n$L22\n$L23\n$L24\n";
			}
			next;
		}
		if($R3Lmax>0 && $L22=~/$R3prefix([$CharN]{$R3Lmin,$R3Lmax})$/){
			$R3Barcode = $1;
			$Barcode .="$R3Barcode:";
		}elsif($R3Lmax>0){
			print OUT1 "Error, wrong R-3Barcode $L21\t$L22\n";
			if($Flagv!=0){
				print OUT6 "$L11\n$L12\n$L13\n$L14\n";
				print OUT7 "$L21\n$L22\n$L23\n$L24\n";
			}
			next;
		}
		my $L5Len = length($L5Barcode);
		my $L3Len = length($L3Barcode);
		my @name = split//,$L11;
		shift @name;
		my $name = join("","@",$Barcode,@name);
		my $LSeq = substr $L12,$L5Len;
		my $LQ = substr $L14,$L5Len;
		$LSeq = substr $LSeq,0,(length($LSeq)-$L3Len);
		$LQ = substr $LQ,0,(length($LSeq)-$L3Len);
		$LRead = substr $LSeq,0,10;
		print OUT3 "$name\n$LSeq\n$L13\n$LQ\n" if($FlagVall==0);

		my $R5Len = length($R5Barcode);
		my $R3Len = length($R3Barcode);
		@name = split//,$L21;
		shift @name;
		my $name = join("","@",$Barcode,@name);
		my $RSeq = substr $L22,$R5Len;
		my $RQ = substr $L24,$R5Len;
		$RSeq = substr $RSeq,0,(length($RSeq)-$R3Len);
		$RQ = substr $RQ,0,(length($RSeq)-$R3Len);
		$RRead = substr $RSeq,0,10;
		print OUT4 "$name\n$RSeq\n$L23\n$RQ\n" if($FlagVall==0);

		if($FlagR5fq==1){
			my $R5BarcodeQ=substr $L24,0,$R5Len;
			print OUT8 "$name\n$R5Barcode\n$L23\n$R5BarcodeQ\n";
		}

		$UniqFrag{"$LRead\t$RRead"}++;
		$UniqUMI{'0'}{"$Barcode"}++;
		my @TempB = split/:/,$Barcode;
		$BarcodeNum=scalar(@TempB);
		for(my$i=1;$i<=$BarcodeNum;$i++){
			$UniqUMI{$i}{$TempB[$i-1]}++;
		}
#	print OUT1 "$L5Barcode\t$Read3\t$len\t$R5Barcode\n";
		$Line++;
		if($Line%1000000 == 0){
			print "read $Line\n";
			my $uniqFrag = scalar(keys %UniqFrag);
			print OUT5 "$Line\t$uniqFrag\n";
			print "$Line\t$uniqFrag\n";

			my $out2;
			for(my$i=0;$i<=$BarcodeNum;$i++){
				my $uniqUMI = scalar(keys %{$UniqUMI{$i}}) + 0;
				$out2 .= "$uniqUMI\t";
			}
			chop $out2;
			print OUT2 "$Line\t$out2\n";
			print "$Line\t$out2\n";
		}
		$Barcode{"$LRead\t$RRead"}{$Barcode}++;
	}
	close IN1;
	close IN2;
}



my $uniqFrag = scalar(keys %UniqFrag);
print OUT5 "$Line\t$uniqFrag\n";
print "$Line\t$uniqFrag\n";

my $out2;
for(my$i=0;$i<=$BarcodeNum;$i++){
	my $uniqUMI = scalar(keys %{$UniqUMI{$i}}) + 0;
	$out2 .= "$uniqUMI\t";
}
chop $out2;
print OUT2 "$Line\t$out2\n";
print "$Line\t$out2\n";

foreach my $r(keys %Barcode){
	foreach my $b(keys %{$Barcode{$r}}){
		print OUT "$r\t$b\t$Barcode{$r}{$b}\n";
	}
}


close OUT;
close OUT1;
close OUT2;
close OUT3;
close OUT4;
close OUT5;
close OUT6;
close OUT7;
close OUT8;

