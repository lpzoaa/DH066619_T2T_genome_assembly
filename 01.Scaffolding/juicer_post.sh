#!/bin/bash
#SBATCH -p xhacnormalb
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
juicer post -o DH066619 DH066619.review.assembly out_JBAT.liftover.agp ../../00.ref/DH066619.1000.fa
