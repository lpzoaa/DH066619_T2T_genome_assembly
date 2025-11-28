#!/bin/bash
#SBATCH -p xhacnormalb
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8

cd ${1}

fastp.0.23.1 -i ~/project/onion_pangenome/rawdata/${1}/hic/${2}_1.fq.gz -I ~/project/onion_pangenome/rawdata/${1}/hic/${2}_2.fq.gz -o ${2}_1.filtered.fq.gz -O ${2}_2.filtered.fq.gz -j ${2}.json -h ${2}.html
