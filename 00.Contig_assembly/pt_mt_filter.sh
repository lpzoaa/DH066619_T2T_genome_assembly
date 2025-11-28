#!/bin/bash
#SBATCH -p xhhcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4

seqkit seq -M 1000000 ${1}.ctg.filtered.fa >${1}.1M.contig.fa

minimap2 -cx asm5 -t 4 ${1}.1M.contig.fa mt.fa >mt.paf
minimap2 -cx asm5 -t 4 ${1}.1M.contig.fa pt.fa >pt.paf

awk '$10/$7>0.8' mt.paf |cut -f 6 >mt.contig.list
awk '$10/$7>0.8' pt.paf |cut -f 6 >pt.contig.list

seqkit grep -f <(cat mt.contig.list pt.contig.list) -v ${1}.ctg.filtered.fa >${1}.ctg.pt_mt.filtered.fa
