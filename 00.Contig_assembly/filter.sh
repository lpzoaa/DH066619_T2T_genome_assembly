#!/bin/bash
#SBATCH -p xhhcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
awk  '$1=="S"' ${1}.bp.p_ctg.noseq.gfa |perl -p -e 's/LN\:i\://;s/rd\:i\://'|awk '$4>50000&&$5>24'|cut -f 2 >${1}.filtered.contig.list
