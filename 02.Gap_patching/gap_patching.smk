###############################
# Snakemake workflow for two-round gap patching
# Round 1: patch with HiFi/UTG contigs
# Round 2: patch with ONT contigs (only if gaps remain)
###############################

import os

configfile: "config.yaml"

RAW_GENOME   = config["raw_genome"]
RAW_AGP      = config["raw_agp"]
UTG_CTG      = config["utg_ctg"]
ONT_CTG      = config["ont_ctg"]
PTG_ROUND1   = config["ptg_round1"]
CHR_LIST     = config["chr_list"]

DE_utg_THR   = config.get("de_utg_threshold", 0.001)
DE_ont_THR   = config.get("de_ont_threshold", 0.01)
MS_THR       = config.get("ms_threshold", 15000)
MAX_GAP_MAG  = config.get("max_gap_margin", 5000000)
THREADS      = config.get("threads", 8)

# Load chromosome IDs (AGP column 1) from chr_list
CHROMS = [
    l.strip() for l in open(CHR_LIST)
    if l.strip() and not l.startswith("#")
]


###############################
# Final target
###############################

rule all:
    """
    Final target: patched genome after two rounds of gap patching.
    """
    input:
        "round2/final.genome.fa"


###############################
# Indexing rules
###############################

rule index_raw_genome:
    """
    Build a faidx index for the raw genome assembly.
    """
    input:
        RAW_GENOME
    output:
        RAW_GENOME + ".fai"
    shell:
        "samtools faidx {input}"

rule index_utg:
    """
    Build a faidx index for the HiFi/UTG contigs.
    """
    input:
        UTG_CTG
    output:
        UTG_CTG + ".fai"
    shell:
        "samtools faidx {input}"

rule index_ont:
    """
    Build a faidx index for the ONT contigs.
    """
    input:
        ONT_CTG
    output:
        ONT_CTG + ".fai"
    shell:
        "samtools faidx {input}"

rule index_ptg_round1:
    """
    Build a faidx index for the primary contig reference used in round 1.
    """
    input:
        PTG_ROUND1
    output:
        PTG_ROUND1 + ".fai"
    shell:
        "samtools faidx {input}"


###############################
# Round 1: candidate enrichment (UTG)
###############################

rule gap_candidates_round1:
    """
    Enrich candidate CTG (UTG) contigs for gap patching using flanking mappings.
    Restricted to chromosomes listed in the user-provided chromosome list.
    """
    input:
        agp      = RAW_AGP,
        ctg      = UTG_CTG,
        ctg_fai  = UTG_CTG + ".fai",
        ptg      = PTG_ROUND1,
        chr_list = CHR_LIST
    output:
        "round1/utg.Flanking_region.list"
    params:
        prefix = "round1/utg",
        de     = DE_utg_THR,
        ms     = MS_THR,
        maxgap = MAX_GAP_MAG,
    threads: 32
    shell:
        r"""
        mkdir -p round1

        gap_patch_candidate_enrichment.pl \
          -agp {input.agp} \
          -ctg {input.ctg} \
          -ptg {input.ptg} \
          -chr-list {input.chr_list} \
          -prefix {params.prefix} \
          -threads {threads} \
          --de {params.de} \
          --ms {params.ms} \
          --max-gap {params.maxgap}
        """


###############################
# Round 1: per-chromosome reference extraction
###############################

rule extract_ref_chr_round1:
    """
    Construct a chromosome-specific reference by concatenating all
    FASTA records that belong to the same logical chromosome.
    Records whose IDs equal the chromosome ID, or start with 'chr:'
    or 'chr_' are merged using 'samtools faidx -r'.
    """
    input:
        fa  = RAW_GENOME,
        fai = RAW_GENOME + ".fai"
    output:
        "round1/ref/{chr}.fa"
    params:
        chr = lambda wildcards: wildcards.chr
    shell:
        r"""
        mkdir -p round1/ref

        # Build a list of FASTA record IDs belonging to this chromosome
        awk -v C={params.chr} '($1==C) || index($1, C":")==1 || index($1, C"_")==1 {{print $1}}' {input.fai} > round1/ref/{wildcards.chr}.ids

        # Extract and concatenate all relevant records into a single FASTA
        samtools faidx {input.fa} -r round1/ref/{wildcards.chr}.ids > {output}
        """


###############################
# Round 1: per-chromosome UTG candidate extraction
###############################

rule extract_utg_for_chr:
    """
    For each chromosome, extract gap-patching UTG candidates whose
    flanking mappings lie within the AGP-defined components for that
    chromosome. Component IDs are taken from W-type entries in the AGP.
    """
    input:
        agp     = RAW_AGP,
        utg     = UTG_CTG,
        utg_fai = UTG_CTG + ".fai",
        flank   = "round1/utg.Flanking_region.list"
    output:
        "round1/utg/{chr}.utg.fa"
    params:
        chr = lambda wildcards: wildcards.chr
    shell:
        r"""
        mkdir -p round1/utg

        # 1) Get component contig IDs (AGP column 6) associated with this chromosome
        #    and with component type W (sequence).
        awk -v C={params.chr} '$1==C && $5=="W" {{print $6}}' {input.agp} \
          | sort -u > round1/utg/{wildcards.chr}.agp.ctg_ids

        # 2) From the UTG flanking region list, select candidate UTG contigs
        #    whose PTG anchors involve these AGP contigs.
        grep -F -f round1/utg/{wildcards.chr}.agp.ctg_ids {input.flank} \
          | cut -f 1 \
          | sort -u > round1/utg/{wildcards.chr}.utg_ids

        # 3) Extract the corresponding UTG contig sequences
        samtools faidx {input.utg} -r round1/utg/{wildcards.chr}.utg_ids > {output}
        """


###############################
# Round 1: ragtag patch per chromosome
###############################

rule ragtag_round1:
    """
    Perform the first round of homology-based patching using UTG contigs
    for each chromosome.
    """
    input:
        ref = "round1/ref/{chr}.fa",
        utg = "round1/utg/{chr}.utg.fa"
    output:
        agp   = "round1/{chr}_round1/ragtag.patch.agp",
        fasta = "round1/{chr}_round1/ragtag.patch.fasta"
    params:
        threads = 16
    shell:
        r"""
        mkdir -p round1/{wildcards.chr}_round1
        
        source activate && conda activate ragtag
        ragtag.py patch \
          --debug \
          -i 0.99 \
          --remove-small \
          -q 10 \
          -d 5000000 \
          -u \
          --aligner minimap2 \
          --mm2-params "-cx asm5 -t {params.threads}" \
          -o round1/{wildcards.chr}_round1 \
          {input.ref} \
          {input.utg}

        filter_paf_by_de.pl round1/{wildcards.chr}_round1/ragtag.patch.asm.paf 

        ragtag.py patch \
          --debug \
          -i 0.99 \
          --remove-small \
          -q 10 \
          -d 5000000 \
          -u \
          --aligner minimap2 \
          --mm2-params "-cx asm5 -t {params.threads}" \
          -o round1/{wildcards.chr}_round1 \
          {input.ref} \
          {input.utg}
        """


###############################
# Round 1: split patched scaffolds & merge
###############################

rule split_scaffold_round1:
    """
    Split patched scaffolds for each chromosome at gap positions
    and rename resulting contigs in a chromosome-specific manner.
    """
    input:
        agp   = "round1/{chr}_round1/ragtag.patch.agp",
        fasta = "round1/{chr}_round1/ragtag.patch.fasta"
    output:
        agp   = "round1/{chr}.round1.agp",
        fasta = "round1/{chr}.round1.contig.fa"
    params:
        chr = lambda wildcards: wildcards.chr
    shell:
        r"""
        split_scaffold_by_agp.pl \
          {input.agp} \
          {input.fasta} \
          {params.chr} \
          {output.agp} \
          {output.fasta}
        """

rule merge_round1:
    """
    Merge per-chromosome round 1 AGP and contig FASTA into genome-wide
    round 1 assemblies.
    """
    input:
        agps   = expand("round1/{chr}.round1.agp",       chr=CHROMS),
        fastas = expand("round1/{chr}.round1.contig.fa", chr=CHROMS)
    output:
        agp   = "round1/all.round1.agp",
        fasta = "round1/all.round1.contig.fa"
    shell:
        r"""
        cat {input.agps}   > {output.agp}
        cat {input.fastas} > {output.fasta}
        """


###############################
# Round 2: candidate enrichment (ONT)
###############################

rule gap_candidates_round2:
    """
    Enrich candidate ONT contigs for a second round of gap patching,
    using the round 1 patched assembly as the PTG reference.
    """
    input:
        agp      = "round1/all.round1.agp",
        ctg      = ONT_CTG,
        ctg_fai  = ONT_CTG + ".fai",
        ptg      = "round1/all.round1.contig.fa",
        chr_list = CHR_LIST
    output:
        "round2/ont.Flanking_region.list"
    params:
        prefix = "round2/ont",
        de     = DE_ont_THR,
        ms     = MS_THR,
        maxgap = MAX_GAP_MAG,
    threads: 32
    shell:
        r"""
        mkdir -p round2

        gap_patch_candidate_enrichment.pl \
          -agp {input.agp} \
          -ctg {input.ctg} \
          -ptg {input.ptg} \
          -chr-list {input.chr_list} \
          -prefix {params.prefix} \
          -threads {threads} \
          --de {params.de} \
          --ms {params.ms} \
          --max-gap {params.maxgap}
        """


###############################
# Round 2: per-chromosome ONT candidate extraction
###############################

rule extract_ont_for_chr:
    """
    For each chromosome, extract gap-patching ONT contigs based on
    the round 1 AGP structure and flanking-region mappings.

    If the round 1 AGP for a chromosome contains no N/U gaps,
    an empty ONT FASTA is created and no second-round patching
    is performed for that chromosome.
    """
    input:
        agp     = "round1/all.round1.agp",
        ont     = ONT_CTG,
        ont_fai = ONT_CTG + ".fai",
        flank   = "round2/ont.Flanking_region.list"
    output:
        "round2/ont/{chr}.ont.fa"
    params:
        chr = lambda wildcards: wildcards.chr
    shell:
        r"""
        mkdir -p round2/ont

        # Check whether this chromosome still has any N/U gaps after round 1
        if awk -v C={params.chr} '$1==C && ($5=="N" || $5=="U")' {input.agp} | grep -q .; then
            # --- There are still gaps: extract ONT candidates ---

            # 1) Get component contig IDs (W-type) for this chromosome in the round 1 AGP
            awk -v C={params.chr} '$1==C && $5=="W" {{print $6}}' {input.agp} \
              | sort -u > round2/ont/{wildcards.chr}.agp.ctg_ids

            # 2) Identify ONT contigs whose flanks anchor into these components
            grep -F -f round2/ont/{wildcards.chr}.agp.ctg_ids {input.flank} \
              | cut -f 1 \
              | sort -u > round2/ont/{wildcards.chr}.ont_ids

            # 3) Extract the ONT contig sequences if any candidates exist
            if [ -s round2/ont/{wildcards.chr}.ont_ids ]; then
                samtools faidx {input.ont} -r round2/ont/{wildcards.chr}.ont_ids > {output}
            else
                # No ONT candidates found: create an empty FASTA as a placeholder
                : > {output}
            fi
        else
            # --- No gaps after round 1: skip ONT extraction for this chromosome ---
            : > {output}
        fi
        """


###############################
# Round 2: ragtag patch per chromosome
###############################

rule ragtag_round2:
    """
    Perform the second round of homology-based patching using ONT contigs
    for each chromosome, starting from the round 1 patched reference.

    If no gaps remain after round 1 (or no ONT candidates are available),
    round 1 patched sequences are propagated directly as round 2 outputs
    for that chromosome.
    """
    input:
        ref = "round1/{chr}_round1/ragtag.patch.fasta",
        ont = "round2/ont/{chr}.ont.fa"
    output:
        fasta = "round2/{chr}_round2/ragtag.patch.fasta",
        agp   = "round2/{chr}_round2/ragtag.patch.agp"
    params:
        threads = 16
    shell:
        r"""
        mkdir -p round2/{wildcards.chr}_round2

        if [ ! -s {input.ont} ]; then
            # No second-round patching: propagate round 1 outputs
            cp {input.ref}  {output.fasta}
            cp round1/{wildcards.chr}_round1/ragtag.patch.agp {output.agp}
        else
            # Perform second-round patching with ONT contigs
            source activate && conda activate ragtag
            ragtag.py patch \
              --debug \
              -i 0.99 \
              --remove-small \
              -q 10 \
              -d 5000000 \
              -u \
              --aligner minimap2 \
              --mm2-params "-cx asm5 -t {params.threads}" \
              -o round2/{wildcards.chr}_round2 \
              {input.ref} \
              {input.ont}

            filter_paf_by_de.pl round1/{wildcards.chr}_round1/ragtag.patch.asm.paf

            ragtag.py patch \
              --debug \
              -i 0.99 \
              --remove-small \
              -q 10 \
              -d 5000000 \
              -u \
              --aligner minimap2 \
              --mm2-params "-cx asm5 -t {params.threads}" \
              -o round2/{wildcards.chr}_round2 \
              {input.ref} \
              {input.ont}

            # Ensure outputs are in the expected locations/names
              if [ ! "{output.fasta}" -ef "round2/{wildcards.chr}_round2/ragtag.patch.fasta" ]; then
                cp round2/{wildcards.chr}_round2/ragtag.patch.fasta {output.fasta}
                cp round2/{wildcards.chr}_round2/ragtag.patch.agp   {output.agp}
              fi
        fi
        """

rule collapse_chr_round2:
    """
    For each chromosome, collapse all RagTag objects in the round 2
    patching result into a single chromosome sequence and AGP.
    """
    input:
        agp   = "round2/{chr}_round2/ragtag.patch.agp",
        fasta = "round2/{chr}_round2/ragtag.patch.fasta"
    output:
        agp   = "round2/{chr}.final.agp",
        fasta = "round2/{chr}.final.fa"
    params:
        chr_id = lambda wildcards: wildcards.chr
    shell:
        r"""
        collapse_chr_from_agp.pl \
          {params.chr_id} \
          {input.agp} \
          {input.fasta} \
          {output.agp} \
          {output.fasta}
        """



###############################
# Final merged genome after round 2
###############################
rule merge_round2:
    """
    Concatenate per-chromosome collapsed round 2 sequences and AGPs
    into the final genome assembly.
    """
    input:
        fastas = expand("round2/{chr}.final.fa",  chr=CHROMS),
        agps   = expand("round2/{chr}.final.agp", chr=CHROMS)
    output:
        fasta = "round2/final.genome.fa",
        agp   = "round2/final.genome.agp"
    shell:
        r"""
        # Merge FASTA
        cat {input.fastas} > {output.fasta}
        #
        # Merge AGP
        cat {input.agps} > {output.agp}
        """
        #
