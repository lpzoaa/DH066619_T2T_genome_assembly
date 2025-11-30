#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

# ----------------------------------------------------------------------
# Enrich candidate contigs for gap patching using flanking mappings.
#
# This script:
#   1) Parses an AGP file and collects contig placement information
#      for a user-specified set of chromosomes/scaffolds (column 1).
#   2) Extracts 20 kb flanking regions from long contigs (ctg FASTA).
#   3) Maps these flanks to a primary contig (PTG) assembly using minimap2.
#   4) Identifies contigs whose left and right flanks both map to the same
#      scaffold region, and estimates the inferred gap span.
#   5) Outputs a list of candidate contigs suitable for gap patching.
#
# Usage:
#   perl gap_patch_candidate_enrichment.pl \
#     -agp      raw.agp \
#     -ctg      ctg.fa \
#     -ptg      ptg.fa \
#     -chr-list chr_ids.txt \
#     -threads  32 \
#     -prefix   utg \
#     --de      0.001 \
#     --ms      15000 \
#     --max-gap 5000000
#
# Required:
#   -agp       AGP file describing scaffold structure
#   -ctg       contig FASTA used as gap-filling candidates
#   -ptg       primary contig (PTG) reference FASTA
#   -chr-list  text file with one chromosome/scaffold ID per line.
#              Only AGP entries whose column 1 is in this list are used.
#   -prefix    prefix for intermediate and output files
#
# Optional:
#   -threads   number of minimap2 threads          [1]
#   --de       maximum allowed per-base error rate in alignments [0.001]
#   --ms       minimum ms:i (score) threshold      [15000]
#   --max-gap  maximum allowed (gap_span - ctg_len) margin [5000000]
#
# Output:
#   <prefix>.Flanking_region.list
#     Columns:
#       1: ctg_id
#       2: left_PTG_id
#       3: left_PTG_pos
#       4: right_PTG_id
#       5: right_PTG_pos
#       6: inferred_gap_length (left + mid + right)
#       7: ctg_length (from ctg.fa.fai)
# ----------------------------------------------------------------------

my ($agp_file, $ctg_fa, $threads, $ptg_fa, $prefix,
    $chr_list_file, $de_thr, $ms_thr, $max_gap_margin);

$threads       = 1;
$de_thr        = 0.001;
$ms_thr        = 15000;
$max_gap_margin = 5_000_000;

GetOptions(
    'agp=s'      => \$agp_file,
    'ctg=s'      => \$ctg_fa,
    'ptg=s'      => \$ptg_fa,
    'chr-list=s' => \$chr_list_file,
    'threads=i'  => \$threads,
    'prefix=s'   => \$prefix,
    'de=f'       => \$de_thr,
    'ms=i'       => \$ms_thr,
    'max-gap=i'  => \$max_gap_margin,
) or die "Error in command line arguments\n";

die <<"USAGE" unless $agp_file && $ctg_fa && $ptg_fa && $chr_list_file && $prefix;
Usage: $0 \\
  -agp      <AGP_file> \\
  -ctg      <contig_fasta> \\
  -ptg      <PTG_fasta> \\
  -chr-list <chr_id_list.txt> \\
  -prefix   <output_prefix> \\
  [-threads INT] [--de FLOAT] [--ms INT] [--max-gap INT]

USAGE

# ----------------------------------------------------------------------
# 1. Read chromosome/scaffold IDs to be retained (AGP column 1 filter)
# ----------------------------------------------------------------------
my %keep_chr;
open(my $CH, "<", $chr_list_file) or die "Cannot open chr-list '$chr_list_file': $!\n";
while (<$CH>) {
    chomp;
    next if /^\s*$/;
    next if /^#/;
    my ($id) = split;
    $keep_chr{$id} = 1;
}
close $CH;

# ----------------------------------------------------------------------
# 2. Parse AGP file and collect contig placement information
#    Only keep W-type components whose object (column 1) is in chr list.
# ----------------------------------------------------------------------
open(my $AGP, "<", $agp_file) or die "Cannot open AGP '$agp_file': $!\n";

my %scaffolds;  # contig_id -> { scaffold, start, end, strand, length }

while (<$AGP>) {
    chomp;
    next if /^#/;
    my @fields = split /\t/;

    # AGP columns:
    # 0: object (scaffold/chromosome ID)
    # 1: object_beg
    # 2: object_end
    # 3: part_number
    # 4: component_type (W/N/U/...)
    # 5: component_id (contig ID)
    # 6: component_beg
    # 7: component_end
    # 8: orientation
    my ($obj_id, $start, $end, $part_number, $component_type,
        $contig_id, $contig_start, $contig_end, $strand) = @fields;

    # Restrict to selected chromosomes/scaffolds
    next unless exists $keep_chr{$obj_id};

    # Only retain W-type components (sequence segments)
    next if $component_type ne 'W';

    $scaffolds{$contig_id} = {
        scaffold => $obj_id,
        start    => $start,
        end      => $end,
        strand   => $strand,
        length   => $contig_end,  # scaffold-level span of this component
    };
}
close $AGP;

# ----------------------------------------------------------------------
# 3. Retrieve contig lengths from the FASTA index (ctg.fa.fai)
# ----------------------------------------------------------------------
my %ctg_length;
my $fai_file = "$ctg_fa.fai";
open(my $FAI, "<", $fai_file) or die "Cannot open FAI '$fai_file': $!\n";
while (<$FAI>) {
    chomp;
    my ($name, $length, @rest) = split /\t/;
    $ctg_length{$name} = $length;
}
close $FAI;

# ----------------------------------------------------------------------
# 4. Extract 20 kb flanking segments from contigs:
#    - 3' terminal 20 kb:   suffix "_R"
#    - 5' initial  20 kb:   suffix "_F"
#    Only contigs with length >= 40 kb are retained.
# ----------------------------------------------------------------------
system("seqkit seq -m 40000 $ctg_fa | seqkit subseq -r -20000:-1 | perl -p -e 'if(/\\>/){s/\\n/_R\\n/}' > ${ctg_fa}.20k_R")
    and die "Error running seqkit for 3' flanks\n";

system("seqkit seq -m 40000 $ctg_fa | seqkit subseq -r 1:20000     | perl -p -e 'if(/\\>/){s/\\n/_F\\n/}' > ${ctg_fa}.20k_F")
    and die "Error running seqkit for 5' flanks\n";

# ----------------------------------------------------------------------
# 5. Map flanking contig segments to PTG reference using minimap2
# ----------------------------------------------------------------------
minimap2($threads, $ptg_fa, "${ctg_fa}.20k_R", $de_thr, $ms_thr);
minimap2($threads, $ptg_fa, "${ctg_fa}.20k_F", $de_thr, $ms_thr);

# ----------------------------------------------------------------------
# 6. Collect primary flank alignments:
#    For each contig:
#      - key   = contig ID
#      - value = PTG target ID
#      - pos   = alignment coordinate on PTG
# ----------------------------------------------------------------------
my $R_list = `perl -p -e 's/\\_R//' ${ctg_fa}.20k_R.filtered.sam | cut -f 1,3,4`;
my $F_list = `perl -p -e 's/\\_F//' ${ctg_fa}.20k_F.filtered.sam | cut -f 1,3,4`;

my %F_hash;
for my $line (split /\n/, $F_list) {
    next if $line =~ /^\s*$/;
    my ($key, $rname, $pos) = split /\t/, $line;

    # Store all PTG hits and their positions for each contig
    push @{$F_hash{$key}}, [$rname, $pos];
}

# ----------------------------------------------------------------------
# 7. Identify flanking pairs mapping to the same scaffold and estimate gaps
# ----------------------------------------------------------------------
open(my $FR, '>', "$prefix.Flanking_region.list")
    or die "Cannot open output '$prefix.Flanking_region.list': $!\n";

for my $line (split /\n/, $R_list) {
    next if $line =~ /^\s*$/;
    my ($key, $rname, $pos) = split /\t/, $line;

    next unless exists $F_hash{$key};

    for my $f_entry (@{$F_hash{$key}}) {
        my ($f_rname, $f_pos) = @$f_entry;

        # Both flanks must map to PTG contigs whose placements are known in AGP,
        # and those PTG contigs must belong to the same scaffold/chromosome.
        next unless exists $scaffolds{$f_rname};
        next unless exists $scaffolds{$rname};
        next unless $scaffolds{$f_rname}{scaffold} eq $scaffolds{$rname}{scaffold};
        next if $f_rname eq $rname;   # ignore self-pairs

        my ($gap_len, $left, $mid, $right);

        # Case 1: left flank precedes right flank on the scaffold
        if ($scaffolds{$rname}{start} > $scaffolds{$f_rname}{end}) {

            $mid = $scaffolds{$rname}{start} - $scaffolds{$f_rname}{end} + 1;

            if ($scaffolds{$f_rname}{strand} eq '+') {
                $left = $scaffolds{$f_rname}{length} - $f_pos + 1;
            } else {
                $left = $f_pos;
            }

            if ($scaffolds{$rname}{strand} eq '+') {
                $right = $pos;
            } else {
                $right = $scaffolds{$rname}{length} - $pos + 1;
            }

            $gap_len = $left + $mid + $right;

        } else {
            # Case 2: reverse ordering of flanks on the scaffold
            $mid = $scaffolds{$f_rname}{start} - $scaffolds{$rname}{end} + 1;

            if ($scaffolds{$f_rname}{strand} eq '+') {
                $left = $f_pos;
            } else {
                $left = $scaffolds{$f_rname}{end} - $scaffolds{$f_rname}{start} - $f_pos + 1;
            }

            if ($scaffolds{$rname}{strand} eq '+') {
                $right = $scaffolds{$rname}{end} - $scaffolds{$rname}{start} - $pos + 1;
            } else {
                $right = $pos;
            }

            $gap_len = $left + $mid + $right;
        }

        # Require that the inferred gap length is not excessively larger than
        # the contig length (allowing a user-defined margin).
        next unless exists $ctg_length{$key};
        if ( ($gap_len - $max_gap_margin) < $ctg_length{$key} ) {
            print $FR join("\t",
                $key,
                $f_rname,
                $f_pos,
                $rname,
                $pos,
                $gap_len,
                $ctg_length{$key}
            ), "\n";
        }
    }
}

close $FR;

# ----------------------------------------------------------------------
# 8. Run minimap2 and retain primary, high-identity, long alignments
#    Output format is SAM (not PAF).
# ----------------------------------------------------------------------
sub minimap2 {
    my ($threads, $ptg, $ctg, $de_thr, $ms_thr) = @_;

    # Use minimap2 in HiFi mapping mode and retain only primary alignments (tp:A:P)
    open my $pipe, '-|', "minimap2 -ax map-hifi -t $threads $ptg $ctg | grep 'tp:A:P'"
        or die "Could not open minimap2 pipe: $!";

    my $out_sam = "${ctg}.filtered.sam";
    open my $out, '>', $out_sam
        or die "Could not open output SAM '$out_sam': $!";

    while (<$pipe>) {
        # SAM optional tags: look for de:f: and ms:i:
        my ($de) = /de:f:([0-9.eE+-]+)/;
        my ($ms) = /ms:i:(\d+)/;

        # Skip records without the required tags
        next unless defined $de && defined $ms;

        # Retain only highly accurate and sufficiently long/high-scoring alignments
        next unless $de < $de_thr;
        next unless $ms > $ms_thr;

        print $out $_;
    }

    close $pipe;
    close $out;
}

