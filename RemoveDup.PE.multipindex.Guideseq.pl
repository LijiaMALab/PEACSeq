#!/usr/bin/perl
BEGIN {
	push @INC, "/home/luzhk/Documents/bin/";
}
use strict;
use Getopt::Long;
my %opts;
my $program=`basename $0`;
chomp $program;
my $usage=<<USAGE; #******* Instruction of this program *********# 

Usage: $program  H1WG.SCP1.H1.R1.fastq.sam[sort by name] key
		remove PE reads with same barcode,same start/end/length
	-L<int>		 length of the index
	-F<int>		 index file
	-b		input bam format
	-p		1/2(F/R. primer direction for guideseq, bydefault 1)
	-s		primer extension sequence
	-help            output help information
USAGE

GetOptions(\%opts, "L:i","b!","F:s","p:i","s:s","help!");
die $usage if ( $opts{'L'} eq "" ||@ARGV==0 || defined($opts{"help"}));

#****************************************************************#
#--------------------Main-----Function-----Start-----------------#
#****************************************************************#

my $infile = shift;
$infile = "samtools view -h $infile|" if(exists $opts{'b'});
open(IN3, $infile) or die $!;
my $outfile = shift;
open(OUT, "|samtools view -Sb - -o $outfile.bam") or die $!;
open(OUT1, ">$outfile.dup") or die $!;
open(OUT2, ">$outfile.nodup.log.summary") or die $!;
open(OUT3, ">$outfile.nodup.wrongIndex") or die $!;
open(OUT4, ">$outfile.unmap") or die $!;
open(OUT5, ">$outfile.guideseq.sam") or die $!;
open(OUT6, ">$outfile.guideseqE.sam") or die $!;

my $len=5;
$len = $opts{'L'} if($opts{'L'} ne "");
print "index length:$len\n";
my $Primer = 1;
$Primer = 2 if($opts{'p'} eq "2");
my $PrimerSeqFlag = 0;
my $PrimerSeq = "";
if($opts{'s'} ne ""){
	$PrimerSeqFlag = 1;
	$PrimerSeq = uc($opts{'s'}) if($opts{'s'} ne "");
}
my @R5S = split//,$PrimerSeq;

my %MyIndex;
if($opts{'F'} ne ""){
	open(IN1, $opts{'F'}) or die $!;
	while (<IN1>) {
		chomp;
		my @info = split/\t/;
		$MyIndex{$info[1]} = 1;
		if($len != length($info[1])){
			print "Error, index has different length\n";
		}
	}
	close IN1;
}

my $head;
my %Count;
my %Reads;
while (<IN3>) {
	if (/^@/){
		$head .="$_";
		print OUT "$_";
		print OUT5 "$_";
		print OUT6 "$_";
	}else{
		my $L1 = $_;
		chomp $L1;
		my @info1 = split/\t/,$L1;
		next if($info1[8]==0);
		my @temp11 = split/:/,$info1[0];
		my $L2 = <IN3>;
		chomp $L2;
		my @info2 = split/\t/,$L2;
		my @temp21 = split/:/,$info2[0];
		if ($info1[0] ne $info2[0]){
			print "Error, $L1 $L2 doen have same name\n";
			next;
		}
		my @temp111 = split//,$temp11[0];
		my @temp211 = split//,$temp21[0];
		$temp11[0]="";
		$temp21[0]="";
		for(my $i=0;$i<$len;$i++){
			$temp11[0] .= "$temp111[$i]";
			$temp21[0] .= "$temp211[$i]";
		}
		$MyIndex{$temp11[0]} = 1 if($opts{'F'} eq "");
		$MyIndex{$temp21[0]} = 1 if($opts{'F'} eq "");
		if ($temp11[0] eq ""){
			print "Error, cant find Index for $info1[0]\n";
			next;
		}
		if ((($info1[1]&2) == 0) ||(($info2[1]&2) == 0)){
			next;
		}
		if ((($info1[1]&4) == 4) ||(($info1[1]&8) == 8)){
			$Count{$temp11[0]}{"unmap\t\t"}++;
#			print OUT4 "$L1\n$L2\n";
			next;
		}
		my $strand="+";
		if(($info1[1]&16) == 0){
			if(($info1[1]&64) == 0){
				$strand = "-";
			}
			my $Len = length($info1[9]);
			$Count{$temp11[0]}{"$info1[2]\t$info1[3]\t$info1[7]\t$strand\t$Len"}++;
			$Reads{"$info1[2]\t$info1[3]\t$info1[7]\t$strand"}{$temp11[0]}{'c'}++;
			if($Reads{"$info1[2]\t$info1[3]\t$info1[7]\t$strand"}{$temp11[0]}{'c'}==1){
				$Reads{"$info1[2]\t$info1[3]\t$info1[7]\t$strand"}{$temp11[0]}{'d'} = "$L1\n$L2";
			}
			if($Count{$temp11[0]}{"$info1[2]\t$info1[3]\t$info1[7]\t$strand\t$Len"} == 1){
				if($MyIndex{$temp11[0]} ==1){
					print OUT "$L1\n$L2\n";
				}else{
					print OUT3 "$L1\n$L2\n";
				}
			}else{
				print OUT1 "$L1\n$L2\n";
			}
		}elsif(($info2[1]&16) == 0){
			if(($info2[1]&64) == 0){
				$strand = "-";
			}
			my $Len = length($info2[9]);
			$Count{$temp21[0]}{"$info2[2]\t$info2[3]\t$info2[7]\t$strand\t$Len"}++;
			$Reads{"$info2[2]\t$info2[3]\t$info2[7]\t$strand"}{$temp21[0]}{'c'}++;
			if($Reads{"$info2[2]\t$info2[3]\t$info2[7]\t$strand"}{$temp21[0]}{'c'}==1){
				$Reads{"$info2[2]\t$info2[3]\t$info2[7]\t$strand"}{$temp21[0]}{'d'} = "$L1\n$L2";
			}
			if($Count{$temp21[0]}{"$info2[2]\t$info2[3]\t$info2[7]\t$strand\t$Len"} == 1){
				if($MyIndex{$temp21[0]} ==1){
					print OUT "$L1\n$L2\n";
				}else{
					print OUT3 "$L1\n$L2\n";
				}
			}else{
				print OUT1 "$L1\n$L2\n";
			}
		}else{
			print "Error, $info1[0] dont have first strand map\n";
		}
	}
}
close IN3;
close OUT;
close OUT1;
close OUT3;
close OUT4;

$"="\t";
my $n=0;
foreach my $p(keys %Reads){
	foreach my $index(keys %{$Reads{$p}}){
		my ($L1,$L2) = split/\n/,$Reads{$p}{$index}{'d'};
		my @info1 = split/\t/,$L1;
		my @info2 = split/\t/,$L2;
		$n++;
		my $PrimerE = $Primer;
		if($PrimerSeqFlag==1){
			my $Seq;
			if(($info1[1]&128)==128){
				$Seq = $info1[9];
				if(($info1[1]&16)==16){
					$Seq =~ tr/ATGCN/TACGN/;
					$Seq = reverse $Seq;
				}
			}elsif(($info2[1]&128)==128){
				$Seq = $info2[9];
				if(($info2[1]&16)==16){
					$Seq =~ tr/ATGCN/TACGN/;
					$Seq = reverse $Seq;
				}
			}else{
				print "Error, cant find sencond strand $L1\n$L2\n";
			}
			my @Temp = split//,$Seq;
			my $R5SFlag=0;
			for(my $i=0;$i<@R5S;$i++){
				if($Temp[$i] ne $R5S[$i]){
					$R5SFlag=1;
				}
			}
			$PrimerE +=2 if($R5SFlag==0);
		}
		my $name = "${index}_${n}_$Reads{$p}{$index}{'c'}_$PrimerE";
		$info1[0] = $name;
		$info2[0] = $name;
		print OUT5 "@info1\n@info2\n";
		print OUT6 "@info1\n@info2\n" if($PrimerE>2);
	}
}
close OUT5;
close OUT6;

my %Index;
foreach my $index(keys %Count){
	foreach my $p(keys %{$Count{$index}}){
		my @info = split/\t/,$p;
		if($info[0] eq "unmap"){
			$Index{$index}{'unmap'} = $Count{$index}{$p};
		}else{
			$Index{$index}{'uniq'}++;
			$Index{$index}{'total'}+=$Count{$index}{$p};
		}
	}
}

foreach my $i (keys %Index){
	print OUT2 "$i\t$Index{$i}{'uniq'}\t$Index{$i}{'total'}\t$Index{$i}{'unmap'}\n";
}
close OUT2;
