# DH066619-T2T genome assembly
[![Snakemake](https://img.shields.io/badge/Snakemake-Workflow-blue.svg)](https://snakemake.github.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

## LG-Patch（Large Genome Patch）
A Snakemake-based pipeline for efficient **gap patching** in large genomes.
The workflow was developed to address the prohibitive computational cost of performing whole-genome alignments on exceptionally large genomes while reducing misalignments caused by abundant repetitive sequences.
Instead of aligning all contigs to entire chromosomes, the pipeline adopts a **contig-end anchoring + candidate enrichment** strategy to identify high-confidence gap-bridging contigs.

## Overview
Gap patching for large  genomes is computationally demanding due to:
- high memory and runtime requirements for chromosome-level alignment
- large numbers of repetitive alignments that obscure true patching candidates
- 

**LG-Patch** 
- Works efficiently on plant genomes reaching tens of gigabases  
- Localized alignment strategy massively reduces computation  
- Modular, reproducible Snakemake workflow  

## Installation
Clone repository:
```
git clone https://github.com/<your_repo>/OnionGapPatch.git
cd OnionGapPatch
```

## Configuration
Edit `config.yaml`:
```
raw_genome: "DH066619_2.FINAL.split.fa"
chr_list: "chr.list"
threads: 32
de_utg_threshold: 0.001
de_ont_threshold: 0.01
```

## Run
```
snakemake -s gap_patching.smk --cores 32
```

## Output
- gap-patched chromosomes  
- updated AGP  
- final merged genome  
