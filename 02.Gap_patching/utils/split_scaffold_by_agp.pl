#!/usr/bin/perl
use strict;
use warnings;

# ----------------------------------------------------------------------
# Split all AGP objects (scaffolds) into contigs at gap (N/U) entries.
#
# For each AGP object (e.g. scf00000000, seq00000005_RagTag):
#   - Consecutive W components with no intervening gaps are merged into a
#     single contig.
#   - Gaps (N/U lines) define contig boundaries.
#   - New contig sequences are extracted from the corresponding object
#     sequence in the FASTA file using the object coordinates (columns 2–3).
#   - A new AGP is written with:
#       * column 1  = chromosome ID derived from a user-specified prefix
#                     (e.g. "scaffold_2", then "scaffold_2_2", "scaffold_2_3"...)
#       * column 6  = new contig IDs (e.g. scaffold_2_1, scaffold_2_2, ...).
#
# Usage:
#   perl split_scaffold_by_agp.pl \
#     <input.agp> <input.fa> <chr_prefix> <output.agp> <output.fa>
#
# Example:
#   perl split_scaffold_by_agp.pl \
#     ragtag.patch.agp ragtag.patch.fasta scaffold_2 \
#     scaffold_2.round1.agp scaffold_2.round1.contig.fa
#
#   In this example:
#     - The first AGP object is mapped to chromosome "scaffold_2".
#     - The second AGP object is mapped to chromosome "scaffold_2_2".
#     - The third AGP object is mapped to chromosome "scaffold_2_3", etc.
# ----------------------------------------------------------------------

if (@ARGV != 5) {
    die "Usage: $0 <input.agp> <input.fa> <chr_prefix> <output.agp> <output.fa>\n";
}

my ($agp_file, $fa_file, $chr_prefix, $out_agp, $out_fa) = @ARGV;

# -----------------------------
# Load all sequences from FASTA
# -----------------------------
my %seq;
{
    open(my $FA, "<", $fa_file) or die "Cannot open FASTA '$fa_file': $!\n";
    my $id = "";
    my $s  = "";

    while (<$FA>) {
        chomp;
        if (/^>(\S+)/) {
            if ($id ne "") {
                $seq{$id} = $s;
            }
            $id = $1;
            $s  = "";
        } else {
            $s .= $_;
        }
    }
    $seq{$id} = $s if $id ne "";
    close $FA;
}

# -----------------------------
# Read AGP and group lines by object (scaffold ID in column 1)
# -----------------------------
my %agp_by_obj;
my @obj_order;

open(my $AGP, "<", $agp_file) or die "Cannot open AGP '$agp_file': $!\n";
while (<$AGP>) {
    chomp;
    next if /^#/;
    next if /^\s*$/;

    my @f = split /\t/;
    my $obj = $f[0];

    if (!exists $agp_by_obj{$obj}) {
        push @obj_order, $obj;
        $agp_by_obj{$obj} = [];
    }
    push @{$agp_by_obj{$obj}}, \@f;
}
close $AGP;

# -----------------------------
# Open output files
# -----------------------------
open(my $OUT_AGP, ">", $out_agp) or die "Cannot open '$out_agp': $!\n";
open(my $OUT_FA,  ">", $out_fa)  or die "Cannot open '$out_fa': $!\n";

# -----------------------------
# Process each AGP object in file order
# -----------------------------
my $obj_index = 0;

for my $obj_id (@obj_order) {
    $obj_index++;

    # Determine the new chromosome (object) ID.
    # First scaffold:  chr_prefix (e.g. "scaffold_2")
    # Second scaffold: chr_prefix_2 (e.g. "scaffold_2_2")
    # Third scaffold:  chr_prefix_3, etc.
    my $new_chr_id;
    if ($obj_index == 1) {
        $new_chr_id = $chr_prefix;
    } else {
        $new_chr_id = $chr_prefix . "_" . $obj_index;
    }

    die "Scaffold '$obj_id' not found in FASTA\n"
        unless exists $seq{$obj_id};

    my $scaf_seq = $seq{$obj_id};

    # ---- Parse this object's AGP lines into contigs and gaps ----
    my @lines = @{$agp_by_obj{$obj_id}};

    # Each contig is defined by object coordinates [start,end] with no gaps inside.
    my @contigs;
    my @gaps;      # gap records in order of appearance
    my $current    = undef;

    for my $f (@lines) {
        my ($obj, $s, $e, $part, $type) = @{$f}[0..4];

        if ($type eq 'W') {
            # Sequence component: extend current contig or start a new one
            if (!defined $current) {
                $current = {
                    start => $s,
                    end   => $e,
                };
            } else {
                # Merge consecutive W blocks into a single contig
                $current->{end} = $e;
            }
        }
        elsif ($type eq 'N' || $type eq 'U') {
            # Gap line: close current contig (if any) and record the gap
            if (defined $current) {
                push @contigs, $current;
                $current = undef;
            }
            my $gap_len  = $f->[5];
            my $gap_type = $type;
            my $info7    = $f->[6] // '';
            my $info8    = $f->[7] // '';
            my $info9    = $f->[8] // '';
            push @gaps, {
                len   => $gap_len,
                type  => $gap_type,
                info7 => $info7,
                info8 => $info8,
                info9 => $info9,
            };
        }
        else {
            # Other component types (if present) are ignored in this splitting logic.
            next;
        }
    }
    # If the object ends with W, close the last contig
    push @contigs, $current if defined $current;

    # ---- Build new contig sequences and write to FASTA & AGP ----
    my $chr_pos   = 1;    # running coordinate along the new chromosome
    my $part_num  = 1;    # AGP part number
    my $ctg_index = 0;    # contig index within this new chromosome
    my $gap_index = 0;    # index into @gaps

    for my $c (@contigs) {
        $ctg_index++;

        my $ctg_start = $c->{start};
        my $ctg_end   = $c->{end};
        my $ctg_len   = $ctg_end - $ctg_start + 1;

        # Sanity check
        if ($ctg_len <= 0) {
            warn "Skipping contig with non-positive length: $obj_id:$ctg_start-$ctg_end\n";
            next;
        }

        # Extract contig sequence from the object sequence using object coordinates
        my $ctg_seq = substr($scaf_seq, $ctg_start - 1, $ctg_len);

        # Define a new contig ID: <new_chr_id>_<ctg_index>, e.g. scaffold_2_1, scaffold_2_2, ...
        my $new_ctg_id = $new_chr_id . "_" . $ctg_index;

        # Write FASTA (wrap at 60 nt per line)
        print $OUT_FA ">$new_ctg_id\n";
        my $len = length($ctg_seq);
        for (my $i = 0; $i < $len; $i += 60) {
            print $OUT_FA substr($ctg_seq, $i, 60), "\n";
        }

        # Write AGP W-line for this contig
        my $chr_start = $chr_pos;
        my $chr_end   = $chr_pos + $ctg_len - 1;

        print $OUT_AGP join("\t",
            $new_chr_id,   # column 1: new chromosome ID
            $chr_start,    # column 2: object start
            $chr_end,      # column 3: object end
            $part_num,     # column 4: part number
            'W',           # column 5: component type
            $new_ctg_id,   # column 6: component ID (new contig)
            1,             # column 7: component start
            $ctg_len,      # column 8: component end
            '+',           # column 9: orientation
        ), "\n";

        $chr_pos  = $chr_end + 1;
        $part_num++;

        # If there is a recorded gap after this contig, reflect it in the new AGP
        if ($gap_index < @gaps) {
            my $g     = $gaps[$gap_index++];
            my $g_len = $g->{len};

            my $g_start = $chr_pos;
            my $g_end   = $chr_pos + $g_len - 1;

            print $OUT_AGP join("\t",
                $new_chr_id,   # object
                $g_start,
                $g_end,
                $part_num,
                $g->{type},    # N or U
                $g->{len},     # gap length
                $g->{info7},   # linkage or gap_type
                $g->{info8},   # yes/no
                $g->{info9},   # additional info
            ), "\n";

            $chr_pos  = $g_end + 1;
            $part_num++;
        }
    }
}

close $OUT_AGP;
close $OUT_FA;

