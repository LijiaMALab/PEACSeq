# coding:utf-8

# 配置文件路径
configfile: "config.yaml"

# print(config)
# 读取配置文件信息
Bin = config['Bin']
Genome = config['Genome']
Genomedir = config['Genomedir']
ChromSize = config['ChromSize']
targets = config['Targets']
primer_F = config['primer_F']
primer_R = config['primer_R']
# adapter_F = config['adapter_F']
# adapter_R = config['adapter_R']

adapters = config['adapters']
thread = config['thread']

# 软件IGV安装路径
igvtools = config['igvtools']

# 所有正向反向测序样本
SAMPLES = [f'{i}.{p}' for i in config['Samples'] for p in [primer_F, primer_R]]

# snakemake rules
rule all:
    input:
        expand("{sample_combine}.{primer_F}.{primer_R}.PEACSeq.DSB2", sample_combine=config['Samples'], primer_F=primer_F, primer_R=primer_R)

rule  barcode_identify:
    input:
        fwd=expand("data/{sample_combine}.1.fastq.gz", sample_combine=SAMPLES),
        rev=expand("data/{sample_combine}.2.fastq.gz", sample_combine=SAMPLES),
    output:
        fwd="barcode_identify/{sample_combine}.Barcode.1.fastq.gz",
        rev="barcode_identify/{sample_combine}.Barcode.2.fastq.gz"
    shell:
        """
        mkdir -p barcode_identify
        perl {Bin}/ParseSampleBarcode.53End.PE.pl data/{wildcards.sample_combine}.1.fastq.gz data/{wildcards.sample_combine}.2.fastq.gz {output.fwd} {output.rev} -L5 6,8,AGAT -gz -N
        """

rule cutadapt:
    input:
        fwd="barcode_identify/{sample_combine}.Barcode.1.fastq.gz",
        rev="barcode_identify/{sample_combine}.Barcode.2.fastq.gz"
    output:
        fwd_trim=temp("cutadapt/{sample_combine}.1.fastq.clipper1.gz"),
        rev_trim=temp("cutadapt/{sample_combine}.2.fastq.clipper1.gz"),
        fwd_untrim=temp("cutadapt/{sample_combine}.1.fastq.clipper2.gz"),
        rev_untrim=temp("cutadapt/{sample_combine}.2.fastq.clipper2.gz")
    log: "cutadapt/{sample_combine}.cutadapt"
    threads: thread
    params:
        five_prime_end_read1_adapter = lambda wildcards : adapters[wildcards.sample_combine]['five_prime_end_read1_adapter'],
        five_prime_end_read2_adapter = lambda wildcards : adapters[wildcards.sample_combine]['five_prime_end_read2_adapter']
    shell:
        """
        mkdir -p cutadapt
        cutadapt -j {threads} -n 1 -e 0.2 -O 12 -m 16 -g ^{params.five_prime_end_read1_adapter} -G ^{params.five_prime_end_read2_adapter}  -o {output.fwd_trim} -p {output.rev_trim} --untrimmed-output={output.fwd_untrim} --untrimmed-paired-output={output.rev_untrim} {input.fwd} {input.rev} > {log}.log
        """


rule cutadapt2:
    input:
        fwd_trim="cutadapt/{sample_combine}.1.fastq.clipper1.gz",
        rev_trim="cutadapt/{sample_combine}.2.fastq.clipper1.gz"
    output:
        fwd_trim="cutadapt/{sample_combine}.1.fastq.clipper.gz",
        rev_trim="cutadapt/{sample_combine}.2.fastq.clipper.gz",
    params:
        three_prime_end_read1_adapter = lambda wildcards : adapters[wildcards.sample_combine]['three_prime_end_read1_adapter'],
        three_prime_end_read2_adapter = lambda wildcards : adapters[wildcards.sample_combine]['three_prime_end_read2_adapter'],
    log:
        "cutadapt/{sample_combine}.cutadapt2"
    threads: thread
    shell:
        """
        cutadapt -j {threads} -n 1 -e 0.1 -O 2 -m 16 -a {params.three_prime_end_read1_adapter} -A {params.three_prime_end_read2_adapter} -o {output.fwd_trim} -p {output.rev_trim}  {input.fwd_trim} {input.rev_trim} > {log}.log
        """

rule bwamem:
    input:
        fwd="cutadapt/{sample_combine}.1.fastq.clipper.gz",
        rev="cutadapt/{sample_combine}.2.fastq.clipper.gz"
    output:
        bam="bwamem/{sample_combine}.PEbwa.bam",
    threads: thread
    shell:
        """
        mkdir -p bwamem
        bwa mem -t {threads} {Genomedir}/{Genome} {input.fwd} {input.rev} | samtools view -Sb -F 3332 -o {output.bam}
        """

rule guidseqBamDedup:
    input:
        bam="bwamem/{sample_combine}.PEbwa.bam",
    output:
        sam=temp("guidseqBamDedup/{sample_combine}.PEbwa.bam.nodup.guideseq.sam"),
        bam="guidseqBamDedup/{sample_combine}.nodup.sort.bam",
        tdf="guidseqBamDedup/{sample_combine}.nodup.tdf",
    log:
        "guidseqBamDedup/{sample_combine}.nodup.bam.log"
    params:
        nodup="guidseqBamDedup/{sample_combine}.PEbwa.bam.nodup",
        chain=lambda wildcards : f"-p 1 -s {config['extend_seq'][wildcards.sample_combine.split('.')[0]][0]}" if primer_F in wildcards.sample_combine else f"-p 2 -s {config['extend_seq'][wildcards.sample_combine.split('.')[0]][1]}"
    threads: thread
    shell:
        """
        mkdir -p guidseqBamDedup
        perl {Bin}/RemoveDup.PE.multipindex.Guideseq.pl -L 7 -b {input.bam} {params.nodup} {params.chain}
        samtools view -Sb {output.sam}|samtools sort - -@ {threads} -o {output.bam}
        samtools flagstat {output.bam} > {log}
        samtools index {output.bam}
        {igvtools} count -z 10 -w 2 {output.bam} {output.tdf} {ChromSize}
        """


rule guidseqAnalysis:
    input:
        sam=expand("guidseqBamDedup/{sample}.{primer}.PEbwa.bam.nodup.guideseq.sam", sample=config['Samples'], primer=[primer_F, primer_R])
    output:
        dbs2="{sample}.{primer_F}.{primer_R}.PEACSeq.DSB2"
    params:
        name = lambda wildcards : targets[wildcards.sample]['name'],
        seq = lambda wildcards : targets[wildcards.sample]['seq'],
    shell:
        """
        mkdir -p outputbwa
        cp guidseqBamDedup/{wildcards.sample}.{primer_F}.PEbwa.bam.nodup.guideseq.sam  {wildcards.sample}-{primer_F}-{primer_R}.bwa.sam
        awk '!/^@/' guidseqBamDedup/{wildcards.sample}.{primer_R}.PEbwa.bam.nodup.guideseq.sam >> {wildcards.sample}-{primer_F}-{primer_R}.bwa.sam
        rm -rf guidseqBamDedup/{wildcards.sample}.{primer_F}.PEbwa.bam.nodup.guideseq.sam guidseqBamDedup/{wildcards.sample}.{primer_R}.PEbwa.bam.nodup.guideseq.sam
        python {Bin}/guideseq.py identify --aligned {wildcards.sample}-{primer_F}-{primer_R}.bwa.sam --genome {Genomedir}/{Genome}.fa --outfolder outputbwa --target_sequence {params.seq} --description {params.name}
        python {Bin}/guideseq.py visualize --infile outputbwa/identified/{wildcards.sample}-{primer_F}-{primer_R}\_identifiedOfftargets.txt --outfolder outputbwa/ --title {wildcards.sample}-{primer_F}-{primer_R}
        cat outputbwa/identified/{wildcards.sample}-{primer_F}-{primer_R}\_identifiedOfftargets.txt| awk -F \"\\t\" '$28!=\"\"||$1~/^#/ {{print $0}}' > {wildcards.sample}.{primer_F}.{primer_R}.PEACSeq
        cat outputbwa/identified/{wildcards.sample}-{primer_F}-{primer_R}\_identifiedOfftargets.txt| awk -F \"\\t\" '$28!=\"\" {{print $0}}'|cut -f 28-30|sort -u > {wildcards.sample}.{primer_F}.{primer_R}.PEACSeq.bed
        perl {Bin}/DNA.Translocation.pl -q 50 -F guidseqBamDedup/{wildcards.sample}.{primer_F}.nodup.sort.bam -R guidseqBamDedup/{wildcards.sample}.{primer_R}.nodup.sort.bam {wildcards.sample}.{primer_F}.{primer_R}.PEACSeq {output.dbs2}
        rm -rf {wildcards.sample}-{primer_F}-{primer_R}.bwa.sam
        """
