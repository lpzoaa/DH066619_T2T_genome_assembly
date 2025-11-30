export PATH=./utils:${PATH}

mkdir -p logs

~/anaconda3/bin/snakemake --snakefile gap_patching.smk  --configfile config.yaml --jobs 40  --cluster-config cluster.yaml --keep-going --cluster "sbatch -p {cluster.queue} -c {cluster.nCPUs} -n 1 -N 1 -o {cluster.output} -e {cluster.error}"
