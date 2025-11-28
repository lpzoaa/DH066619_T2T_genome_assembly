#!/bin/bash
#SBATCH -p xhacnormalb
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32

ln -s ../../00.ref/DH066619.1000.fa .
samtools faidx DH066619.1000.fa
juicer pre -a -q 1 -o out_JBAT ../../02.bam/DH066619.filtered.bam scaffolds.raw.agp DH066619.1000.fa.fai >out_JBAT.log 2>&1
(java -jar -Xmx128G juicer_tools.1.9.9_jcuda.0.8.jar pre out_JBAT.txt out_JBAT.hic.part <(cat out_JBAT.log | grep PRE_C_SIZE | awk '{print $2" "$3}')) && (mv out_JBAT.hic.part out_JBAT.hic)
