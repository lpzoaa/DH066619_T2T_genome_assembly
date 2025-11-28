# Onion Gap Patch
[![Snakemake](https://img.shields.io/badge/Snakemake-Workflow-blue.svg)](https://snakemake.github.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/Platform-Linux%20%7C%20HPC-orange.svg)]()

A Snakemake-based pipeline for efficient **gap patching** in very large plant genomes.
The workflow is designed for assemblies where whole-genome alignment is computationally prohibitive.
Instead of aligning all contigs to entire chromosomes, the pipeline adopts a **contig-end anchoring + candidate enrichment** strategy to identify high-confidence gap-bridging contigs.

## Overview
Gap patching for large eukaryotic genomes is computationally demanding due to:
- high memory and runtime requirements for chromosome-level alignment
- large numbers of repetitive alignments that obscure true patching candidates

**OnionGapPatch** solves this by introducing:
1. **Contig-end extraction**  
2. **Candidate enrichment** using positional logic  
3. **Per-chromosome patching**  
4. **Two-round scaffolding**  
   - Round 1: HiFi/ONT hybrid unitigs  
   - Round 2: ONT-only contigs  

## Key Features
- Works efficiently on plant genomes reaching tens of gigabases  
- Localized alignment strategy massively reduces computation  
- Perl scripts for contig-end extraction, enrichment, PAF filtering, and AGP manipulation  
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

## Citation
If you use this workflow, please cite this repository.

## Contact
pp l
