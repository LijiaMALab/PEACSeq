#########################Alignment,Identification,Visualization###################
genomedir="/Data/Database/hg38/"
genome="GRCh38"
T1_Seq="GGGTGGGGGGAGTTTGCTCCNGG"
T1_Name="VEGFAT1"
T2_Seq="GACCCCCTCCACCCCGCCTCNGG"
T2_Name="VEGFAT2"
T3_Seq="GGTGAGTGAGTGTGTGCGTGNGG"
T3_Name="VEGFAT3"
T4_Seq="GAGTCCGAGCAGAAGAAGAANGG"
T4_Name="EMX1T4"
T5_Seq="GTCATCTTAGTCATTACCTGNGG"
T5_Name="RNF2T5"
T6_Seq="GGAATCCCTTCTGCAGCACCNGG"
T6_Name="FANCFT6"
sample="Cg1 Cg2 Lg10 Lg11 Lg12 Lg13"
primerF="F"
primerR="R"
for s in $sample
do
for F in $primerF
do
for R in $primerR
do
Seq=${s}_Seq
Name=${s}_Name
echo cat ../$s.$F.PEbwa.nodup.guideseq.sam ../$s.$R.PEbwa.nodup.guideseq.sam \> $s-$F-$R.bwa.sam
cat ../$s.$F.PEbwa.nodup.guideseq.sam ../$s.$R.PEbwa.nodup.guideseq.sam > $s-$F-$R.bwa.sam
python /Data2/Data/PEGUIDESeq/guideseqMalab/guideseq/guideseq.py identify --aligned $s-$F-$R.bwa.sam --genome $genomedir/$genome.fa --outfolder outputbwa --target_sequence ${!Seq} --description ${!Name}
python /Data2/Data/PEGUIDESeq/guideseqMalab/guideseq/guideseq.py visualize --infile outputbwa/identified/$s-$F-$R\_identifiedOfftargets.txt --outfolder outputbwa/ --title $s-$F-$R
cat outputbwa/identified/$s-$F-$R\_identifiedOfftargets.txt| awk -F "\t" '$28!=""||$1~/^#/ {print $0}' > $s.$F.$R.PEACSeq
cat outputbwa/identified/$s-$F-$R\_identifiedOfftargets.txt| awk -F "\t" '$28!="" {print $0}'|cut -f28-30|sort -u > $s.$F.$R.PEACSeq.bed
perl /Data2/Data/PEGUIDESeq/DNA.Translocation.pl -b -q 50 -F ../$s.$F.nodup.sort.bam -R ../$s.$R.nodup.sort.bam $s.$F.$R.PEACSeq $s.$F.$R.PEACSeq.DSB3 &
perl /Data2/Data/PEGUIDESeq/DNA.Translocation.pl -q 50 -F ../$s.$F.nodup.sort.bam -R ../$s.$R.nodup.sort.bam $s.$F.$R.PEACSeq $s.$F.$R.PEACSeq.DSB2 &
wait
perl /Data2/Data/PEGUIDESeq/PEACSeq.Nor.pl ../$s.$F.bwa.bam.log ../$s.$R.bwa.bam.log $s.$F.$R.PEACSeq.DSB3 $s.$F.$R.PEACSeq.DSB.nor
done
done
wait
done

sample="Ag1 Ag2 Ag3 Lg1 Lg2 Lg3 Lg4"
sample="Ag4 Ag5 Lg5 Lg6 Lg7 Lg8 Lg9"
sample="Cg1 Cg2 Lg10 Lg11 Lg12 Lg13"
primerF="F"
primerR="R"
for F in $primerF
do
for R in $primerR
do
Contrl="WT ODN"
for c in $Contrl
do
echo cat ../$c.$F.PEbwa.nodup.guideseq.sam ../$c.$R.PEbwa.nodup.guideseq.sam \> $c-$F-$R.bwa.sam
cat ../$c.$F.PEbwa.nodup.guideseq.sam ../$c.$R.PEbwa.nodup.guideseq.sam > $c-$F-$R.bwa.sam
for s in $sample
do
Seq=${s}_Seq
Name=${s}_Name
python /Data2/Data/PEGUIDESeq/guideseqMalab/guideseq/guideseq.py identify --aligned $c-$F-$R.bwa.sam --genome $genomedir/$genome.fa --outfolder outputbwa-$c-$s --target_sequence ${!Seq} --description ${!Name}
python /Data2/Data/PEGUIDESeq/guideseqMalab/guideseq/guideseq.py visualize --infile outputbwa-$c-$s/identified/$c-$F-$R\_identifiedOfftargets.txt --outfolder outputbwa-$c-$s/ --title $c-$F-$R
cat outputbwa-$c-$s/identified/$c-$F-$R\_identifiedOfftargets.txt| awk -F "\t" '$28!=""||$1~/^#/ {print $0}' > $c-$s.$F.$R.PEACSeq
cat outputbwa-$c-$s/identified/$c-$F-$R\_identifiedOfftargets.txt| awk -F "\t" '$28!="" {print $0}'|cut -f28-30|sort -u > $c-$s.$F.$R.PEACSeq.bed

perl /Data2/Data/PEGUIDESeq/DNA.Translocation.pl -b -q 50 -F ../$c.$F.nodup.sort.bam -R ../$c.$R.nodup.sort.bam $c-$s.$F.$R.PEACSeq $c-$s.$F.$R.PEACSeq.DSB3 &
perl /Data2/Data/PEGUIDESeq/DNA.Translocation.pl -q 50 -F ../$s.$F.nodup.sort.bam -R ../$s.$R.nodup.sort.bam $s.$F.$R.PEACSeq $s.$F.$R.PEACSeq.DSB2 &
wait
perl /Data2/Data/PEGUIDESeq/PEACSeq.Nor.pl ../$c.$F.bwa.bam.log ../$c.$R.bwa.bam.log $c-$s.$F.$R.PEACSeq.DSB3 $c-$s.$F.$R.PEACSeq.DSB.nor
done
done
wait
done
done
#exit



sample="Ag1 Ag2 Ag3 Lg1 Lg2 Lg3 Lg4"
sample="Ag4 Ag5 Lg5 Lg6 Lg7 Lg8 Lg9"
sample="Cg1 Cg2 Lg10 Lg11 Lg12 Lg13"
Contrl="WT ODN"
for c in $Contrl
do
primerF="F"
primerR="R"
for s in $sample
do
for F in $primerF
do
for R in $primerR
do
perl /Ubin/bin/fish.pl $c-$s.$F.$R.PEACSeq.DSB.nor $s.$F.$R.PEACSeq.DSB.nor -bait 28,29,30 -fish 28,29,30 -all -cb 1:30,37:48 > $c.$s.$F.$R.PEACSeq.DSB.nor
done
done
done
done
