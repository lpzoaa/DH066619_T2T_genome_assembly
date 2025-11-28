#!/bin/bash
#SBATCH -p xhacnormalb
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 96

mkdir -p $1
cd ${1}

mkdir -p 00.ref
mkdir -p 01.fastq
mkdir -p 02.bam

ln -sf ~/project/onion_pangenome/rawdata/${1}/hic/${1}_1.filtered.fq.gz ./01.fastq/
ln -sf ~/project/onion_pangenome/rawdata/${1}/hic/${1}_2.filtered.fq.gz ./01.fastq/

seqkit sliding -W 1000000000 -s 1000000000 -g -S '' ~/project/onion_pangenome/result/02.assembly/01.contig/${1}/${1}.ctg.pt_mt.filtered.fa >00.ref/${1}.1000.fa

bwa index 00.ref/${1}.1000.fa

bwa mem -5SP -t 128 00.ref/${1}.1000.fa 01.fastq/${1}_1.filtered.fq.gz 01.fastq/${1}_2.filtered.fq.gz|samblaster|samtools view - -@ 128 -S -h -b -F 3340 -o 02.bam/${1}.bam

source activate
conda activate haphic

~/mambaforge/envs/HapHiC/utils/filter_bam 02.bam/${1}.bam 1 --nm 3 --threads 96|samtools view - -b -@ 96 -o 02.bam/${1}.filtered.bam

mkdir -p ${1}_haphic
cd ${1}_haphic

~/mambaforge/envs/HapHiC/haphic pipeline --threads 96 ../00.ref/${1}.1000.fa  ../02.bam/${1}.filtered.bam 8
