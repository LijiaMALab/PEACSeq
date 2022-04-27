bin=/Ubin/bin
Genome=GRCh38
Genomedir=/Data/Database/hg38

ln -s R22008657-293-F-F293-Lg10_combined_R1.fastq.gz Lg10.F.1.fastq.gz
ln -s R22008657-293-F-F293-Lg10_combined_R2.fastq.gz Lg10.F.2.fastq.gz
ln -s R22008658-293-R-R293-ODN_combined_R1.fastq.gz ODN.R.1.fastq.gz
ln -s R22008658-293-R-R293-ODN_combined_R2.fastq.gz ODN.R.2.fastq.gz
ln -s R22008658-293-R-R293-WT_combined_R1.fastq.gz WT.R.1.fastq.gz
ln -s R22008658-293-R-R293-WT_combined_R2.fastq.gz WT.R.2.fastq.gz
#exit

sample="Lg10.F ODN.F WT.F"
for s in $sample
do
perl /Ubin/bin/ParseSampleBarcode.53End.PE.pl $s.1.fastq.gz $s.2.fastq.gz  $s.Barcode -L5 6,8,AGAT -gz -N
###################iguide
cutadapt -n 1 -e 0.2 -O 12 -m 16 -g ^AGATGTGTATAAGAGACAG -G ^GCTCGCGTTTAATTGAGTTGTCATATGTTAATAA -o $s.1.fastq.clipper1.gz -p $s.2.fastq.clipper1.gz --untrimmed-output=$s.1.fastq.clipper2.gz --untrimmed-paired-output=$s.2.fastq.clipper2.gz $s.Barcode.1.fastq.gz $s.Barcode.2.fastq.gz > $s.log
cutadapt -n 1 -e 0.1 -O 2 -m 16 -a TTATTAACATATGACAACTCAATTAAACGCGAGC -A CTGTCTCTTATACACATCT -o $s.1.fastq.clipper.gz -p $s.2.fastq.clipper.gz  $s.1.fastq.clipper1.gz $s.2.fastq.clipper1.gz > $s.log
###perl /Data/PEGUIDESeq/AdapterFilter.pl -gz -R5 3,0.8,TCATAC $s.1.fastq.clipper.gz $s.2.fastq.clipper.gz $s.filter
bwa mem -t 100 $Genomedir/GRCh38 $s.1.fastq.clipper.gz $s.2.fastq.clipper.gz | samtools view -Sb -F 3332 -o $s.PEbwa.bam
perl ../RemoveDup.PE.multipindex.Guideseq.pl -L 7 -b $s.PEbwa.bam $s.PEbwa.nodup -p 1 -s CG
wait
done
sample="Lg10.R ODN.R WT.R"
for s in $sample
do
perl /Ubin/bin/ParseSampleBarcode.53End.PE.pl $s.1.fastq.gz $s.2.fastq.gz  $s.Barcode -L5 6,8,AGAT -gz -N
###################iguide
cutadapt -n 1 -e 0.2 -O 12 -m 16 -g ^AGATGTGTATAAGAGACAG -G ^TCGCGTATACCGTTATTAACATATGACAACTCAAT -o $s.1.fastq.clipper1.gz -p $s.2.fastq.clipper1.gz --untrimmed-output=$s.1.fastq.clipper2.gz --untrimmed-paired-output=$s.2.fastq.clipper2.gz $s.Barcode.1.fastq.gz $s.Barcode.2.fastq.gz > $s.log
cutadapt -n 1 -e 0.1 -O 2 -m 16 -a ATTGAGTTGTCATATGTTAATAACGGTATACGCGA -A CTGTCTCTTATACACATCT -o $s.1.fastq.clipper.gz -p $s.2.fastq.clipper.gz  $s.1.fastq.clipper1.gz $s.2.fastq.clipper1.gz > $s.log
###perl /Data/PEGUIDESeq/AdapterFilter.pl -gz -R5 3,0.8,GGT $s.1.fastq.clipper.gz $s.2.fastq.clipper.gz $s.filter
bwa mem -t 100 $Genomedir/GRCh38 $s.1.fastq.clipper.gz $s.2.fastq.clipper.gz | samtools view -Sb -F 3332 -o $s.PEbwa.bam
perl ../RemoveDup.PE.multipindex.Guideseq.pl -L 7 -b $s.PEbwa.bam $s.PEbwa.nodup -p 2 -s TA
wait
done

sample="Lg10 ODN WT"
for s in $sample
do
FR="F R"
for p in $FR
do
samtools sort $s.$p.PEbwa.bam -o $s.$p.bwa.sort.bam  -@ 10
samtools flagstat $s.$p.bwa.sort.bam > $s.$p.bwa.bam.log
samtools index $s.$p.bwa.sort.bam
igvtools count -z 10 -w 2 $s.$p.bwa.sort.bam $s.$p.bwa.tdf /Data/Database/hg38/GRCh38.chrom.sizes

samtools view -Sb $s.$p.PEbwa.nodup.guideseq.sam|samtools sort - -@ 10 -o $s.$p.nodup.sort.bam
samtools flagstat $s.$p.nodup.sort.bam > $s.$p.nodup.bam.log
samtools index $s.$p.nodup.sort.bam
igvtools count -z 10 -w 2 $s.$p.nodup.sort.bam $s.$p.nodup.tdf /Data/Database/hg38/GRCh38.chrom.sizes

samtools view -q 50 -Sb $s.$p.PEbwa.nodup.guideseq.sam|samtools sort - -@ 20 -o $s.$p.nodupfilter.sort.bam
samtools flagstat $s.$p.nodupfilter.sort.bam > $s.$p.nodupfilter.bam.log
samtools index $s.$p.nodupfilter.sort.bam
igvtools count -z 10 -w 2 $s.$p.nodupfilter.sort.bam $s.$p.nodupfilter.tdf /Data/Database/hg38/GRCh38.chrom.sizes

samtools view -Sb $s.$p.PEbwa.nodup.guideseqE.sam|samtools sort - -@ 10 -o $s.$p.nodupE.sort.bam
samtools flagstat $s.$p.nodupE.sort.bam > $s.$p.nodupE.bam.log
samtools index $s.$p.nodupE.sort.bam
igvtools count -z 10 -w 2 $s.$p.nodupE.sort.bam $s.$p.nodupE.tdf /Data/Database/hg38/GRCh38.chrom.sizes

bwa mem -t 46 ../Genome/T1 $s.$p.1.fastq.clipper.gz $s.$p.2.fastq.clipper.gz | samtools view -Sb -F 3332 -o $s.$p.vector.bam
samtools sort $s.$p.vector.bam -o $s.$p.vector.sort.bam  -@ 10
samtools flagstat $s.$p.vector.sort.bam > $s.$p.vector.bam.log
samtools index $s.$p.vector.sort.bam
igvtools count -z 10 -w 2 $s.$p.vector.sort.bam $s.$p.vector.tdf ../Genome/T1.chrom.sizes

C0=$(zcat $s.$p.1.fastq.gz|wc -l)
C1=$(zcat $s.$p.1.fastq.clipper1.gz|wc -l)
C2=$(zcat $s.$p.1.fastq.clipper2.gz|wc -l)
C3=$(zcat $s.$p.1.fastq.clipper.gz|wc -l)
C4=$(samtools view $s.$p.bwa.sort.bam|wc -l)
C5=$(samtools view $s.$p.nodup.sort.bam|wc -l)
C6=$(samtools view $s.$p.nodupfilter.sort.bam|wc -l)
C7=$(samtools view $s.$p.nodupE.sort.bam|wc -l)
C8=$(samtools view $s.$p.vector.sort.bam|wc -l)
CR0=$((C0 / 4))
CR1=$((C1 / 4))
CR2=$((C2 / 4))
CR3=$((C3 / 4))
CR4=$((C4 / 2))
CR5=$((C5 / 2))
CR6=$((C6 / 2))
CR7=$((C7 / 2))
CR8=$((C8 / 2))
echo -e "start" > $s.$p.summary
echo -e "$s\tTotal\t$CR0" >> $s.$p.summary
echo -e "$s\tclipper1\t$CR1" >> $s.$p.summary
echo -e "$s\tclipper2(removed)\t$CR2" >> $s.$p.summary
echo -e "$s\tclipper(final)\t$CR3" >> $s.$p.summary
echo -e "$s\tAllMapped\t$CR4" >> $s.$p.summary
echo -e "$s\tnodup\t$CR5" >> $s.$p.summary
echo -e "$s\tnodupFilter\t$CR6" >> $s.$p.summary
echo -e "$s\tnodupE\t$CR7" >> $s.$p.summary
echo -e "$s\tVector\t$CR8" >> $s.$p.summary
wait
done
done
more *F.summary *R.summary > 2nd.summary
