#!/usr/bin/env perl
use strict;
use warnings;

# ----------------------------------------------------------------------
# Collapse all objects from a per-chromosome RagTag patch into a
# single chromosome sequence and AGP entry set.
#
# Inputs:
#   1) chr_id       - final chromosome ID to use in AGP/FASTA (e.g. "scaffold_1")
#   2) in.agp       - per-chromosome AGP from ragtag.patch (e.g. scaffold_1_round2/ragtag.patch.agp)
#   3) in.fa        - per-chromosome FASTA from ragtag.patch (e.g. scaffold_1_round2/ragtag.patch.fasta)
#
# Outputs:
#   4) out.agp      - collapsed AGP where:
#                       * column 1 is always chr_id
#                       * object coordinates (columns 2–3) are shifted
#                         by the cumulative length of preceding objects
#   5) out.fa       - collapsed FASTA with a single sequence:
#                       * header: >chr_id
#                       * sequence: concatenation of all input sequences
#                         in the order they appear in in.fa
#
# Usage:
#   collapse_chr_from_agp.pl \
#     scaffold_1 \
#     round2/scaffold_1_round2/ragtag.patch.agp \
#     round2/scaffold_1_round2/ragtag.patch.fasta \
#     round2/scaffold_1.final.agp \
#     round2/scaffold_1.final.fa
# ----------------------------------------------------------------------

if (@ARGV != 5) {
    die "Usage: $0 <chr_id> <in.agp> <in.fa> <out.agp> <out.fa>\n";
}

my ($chr_id, $agp_file, $fa_file, $out_agp, $out_fa) = @ARGV;

# ----------------------------------------------------------------------
# 1. Load all sequences from the per-chromosome FASTA.
#    We also record their order to preserve the concatenation order.
# ----------------------------------------------------------------------
my %seq;        # object_id -> sequence
my @seq_order;  # order of object IDs as they appear in FASTA

{
    open(my $FA, "<", $fa_file) or die "Cannot open FASTA '$fa_file': $!\n";
    my $id = "";
    my $s  = "";

    while (<$FA>) {
        chomp;
        if (/^>(\S+)/) {
            if ($id ne "") {
                $seq{$id} = $s;
                push @seq_order, $id;
            }
            $id = $1;
            $s  = "";
        } else {
            $s .= $_;
        }
    }
    if ($id ne "") {
        $seq{$id} = $s;
        push @seq_order, $id;
    }
    close $FA;
}

# ----------------------------------------------------------------------
# 2. Read AGP and group lines by object (column 1).
# ----------------------------------------------------------------------
my %agp_by_obj;

open(my $AGP, "<", $agp_file) or die "Cannot open AGP '$agp_file': $!\n";
while (<$AGP>) {
    chomp;
    next if /^#/;  # skip AGP comments; we will write a fresh header

    my @f = split /\t/;
    my $obj = $f[0];

    push @{$agp_by_obj{$obj}}, \@f;
}
close $AGP;

# ----------------------------------------------------------------------
# 3. Open output AGP and FASTA for the collapsed chromosome.
# ----------------------------------------------------------------------
open(my $OUT_AGP, ">", $out_agp) or die "Cannot open '$out_agp': $!\n";
open(my $OUT_FA,  ">", $out_fa)  or die "Cannot open '$out_fa':  $!\n";

# Write a simple AGP header
print $OUT_AGP "## agp-version 2.1\n";
print $OUT_AGP "# Collapsed chromosome $chr_id from per-chromosome RagTag patch\n";

# ----------------------------------------------------------------------
# 4. Concatenate sequences and rewrite AGP coordinates.
# ----------------------------------------------------------------------
my $final_seq = "";
my $offset    = 0;   # cumulative length of previous objects
my $part_no   = 1;   # global AGP part number for this chromosome

for my $obj (@seq_order) {

    unless (exists $seq{$obj}) {
        die "Sequence '$obj' missing in FASTA '$fa_file'\n";
    }
    unless (exists $agp_by_obj{$obj}) {
        die "Object '$obj' missing in AGP '$agp_file'\n";
    }

    my $obj_seq = $seq{$obj};
    my $obj_len = length($obj_seq);

    # Append this object sequence to the final chromosome
    $final_seq .= $obj_seq;

    # Rewrite its AGP lines with shifted coordinates and new chromosome ID
    for my $f (@{$agp_by_obj{$obj}}) {
        my @g = @$f;  # copy

        # Columns (0-based):
        # 0: object
        # 1: object_start
        # 2: object_end
        # 3: part_number
        # 4: component_type (W, N, U, etc.)
        my ($old_obj, $s, $e, $old_part, $type) = @g[0..4];

        my $new_start = $s + $offset;
        my $new_end   = $e + $offset;

        $g[0] = $chr_id;     # new chromosome ID
        $g[1] = $new_start;  # shifted start
        $g[2] = $new_end;    # shifted end
        $g[3] = $part_no++;  # renumber parts sequentially

        print $OUT_AGP join("\t", @g), "\n";
    }

    # Advance offset by the full length of this object sequence
    $offset += $obj_len;
}

# ----------------------------------------------------------------------
# 5. Write the final collapsed chromosome sequence to FASTA.
# ----------------------------------------------------------------------
print $OUT_FA ">$chr_id\n";
(my $tmp = $final_seq) =~ s/(.{1,60})/$1\n/g;
print $OUT_FA $tmp;

close $OUT_AGP;
close $OUT_FA;

exit 0;

