#!/usr/bin/env perl
use strict;
use warnings;

# This script filters a PAF file based on the "de:f" tag.
# Only alignments with de:f < 0.001 are retained.
# The original PAF file is overwritten in-place with the filtered content.
#
# Usage:
#   perl filter_paf_by_de.pl input.paf
#
# Notes:
#   - PAF lines without a "de:f" tag are discarded.
#   - Header or comment lines (starting with '#') are preserved as-is.
#   - The script writes to a temporary file in the same directory and then renames it.

my $in_paf = shift or die "Usage: $0 <input.paf>\n";

# Temporary output file (same directory as input)
my $tmp_paf = "$in_paf.tmp.de_lt_0.001";

open my $IN,  '<', $in_paf  or die "Cannot open input PAF file '$in_paf': $!\n";
open my $OUT, '>', $tmp_paf or die "Cannot open temporary file '$tmp_paf': $!\n";

while (my $line = <$IN>) {
    chomp $line;

    # Preserve empty lines and comment/header lines if present
    if ($line =~ /^\s*$/) {
        print $OUT "$line\n";
        next;
    }
    if ($line =~ /^#/) {
        print $OUT "$line\n";
        continue;
    }

    my @f = split /\t/, $line;

    # PAF main fields are 12 columns; tags start from column 13 (index 12)
    my $keep = 0;
    for my $i (12 .. $#f) {
        # Example tag format: de:f:0.0005
        if ($f[$i] =~ /^de:f:(\S+)/) {
            my $de = $1;
            # Ensure it looks numeric
            next unless $de =~ /^-?\d*\.?\d+(e-?\d+)?$/i;
            if ($de < 0.001) {
                $keep = 1;
            }
            last;
        }
    }

    # Only write alignments that pass the de:f threshold
    if ($keep) {
        print $OUT join("\t", @f), "\n";
    }
}

close $IN;
close $OUT;

# Replace the original PAF file with the filtered one
rename $tmp_paf, $in_paf
  or die "Failed to replace '$in_paf' with filtered file '$tmp_paf': $!\n";

