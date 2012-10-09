#! /usr/loca/bin/perl

use strict;
use List::Util qw(max min sum);

my $INPUT_tumor = $ARGV[0];
my $INPUT_normal = $ARGV[1];
my $MIN_TUMOR_DEPTH = $ARGV[2];
my $MIN_NORMAL_DEPTH = $ARGV[3];
my $MIN_TUMOR_VARIANT_READ = $ARGV[4];
my $MIN_TUMOR_ALLELE_FREQ = $ARGV[5];
my $MAX_NORMAL_ALLELE_FREQ = $ARGV[6];



my %base2ind = ();
$base2ind{"A"} = 0;
$base2ind{"C"} = 1;
$base2ind{"G"} = 2;
$base2ind{"T"} = 3;

open(IN, $INPUT_tumor) || die "cannot open $!";
my %ID2bases_tum = ();
while(<IN>) {

	s/[\r\n]//g;
	my @F = split("\t", $_);

	my $ref = $F[2];
    my $idx_ref = $base2ind{$ref};
    my $key = join("\t", @F[0 .. 2]);

    ##########
    my @baseNums = ($F[4] + $F[5], $F[6] + $F[7], $F[8] + $F[9], $F[10] + $F[11]);
    my $depth = sum(@baseNums);
    my @baseNums_p = ($F[4], $F[6], $F[8], $F[10]);
    my @baseNums_n = ($F[5], $F[7], $F[9], $F[11]);

	if ($depth >= $MIN_TUMOR_DEPTH) {

		# determine which is the most frequent mismatch
		##########
		$baseNums[$base2ind{$ref}] = 0;
            
        my $misNum = max(@baseNums);
        my $idx_tum = min(grep{$baseNums[$_] == $misNum} 0 .. $#baseNums);
		my $mis = substr("ACGT", $idx_tum, 1); 
		##########

        if (($misNum / $depth) >= $MIN_TUMOR_ALLELE_FREQ and $misNum >= $MIN_TUMOR_VARIANT_READ) {

           	$ID2bases_tum{$key} = $mis . "\t" . sum(@baseNums_p) . "," . $F[2 * $idx_tum + 4] . "," . sum(@baseNums_n) . "," . $F[2 * $idx_tum + 5];

        }

   	}

}
close(IN);




open(IN, $INPUT_normal) || die "cannot open $!";
my %ID2bases = ();
while(<IN>) {

    s/[\r\n]//g;
    my @F = split("\t", $_);

    my $ref = $F[2];
    my $idx_ref = $base2ind{$ref};
    my $key = join("\t", @F[0 .. 2]);

    ##########
    my @baseNums = ($F[4] + $F[5], $F[6] + $F[7], $F[8] + $F[9], $F[10] + $F[11]);
    my $depth = sum(@baseNums);
    my @baseNums_p = ($F[4], $F[6], $F[8], $F[10]);
    my @baseNums_n = ($F[5], $F[7], $F[9], $F[11]);


	if (exists $ID2bases_tum{$key} and $depth >= $MIN_NORMAL_DEPTH) {

		my @bases_tum = split("\t", $ID2bases_tum{$key});
		my $mis = $bases_tum[0];

        my $idx_tum = $base2ind{$mis};
        my $misNum = $F[2 * $idx_tum + 4] + $F[2 * $idx_tum + 5];
    
        if (($misNum / $depth) < $MAX_NORMAL_ALLELE_FREQ) {
            print $key . "\t" .  $mis . "\t" . $bases_tum[1] . "\t" . sum(@baseNums_p) . "," . $F[2 * $idx_tum + 4] . "," . sum(@baseNums_n) . "," . $F[2 * $idx_tum + 5] . "\n";
            # $ID2bases{$key} = $mis . "\t" . $bases_tum[1] . "\t" . sum(@baseNums_p) . "," . $F[2 * $idx_tum + 4] . "," . sum(@baseNums_n) . "," . $F[2 * $idx_tum + 5];
        }
   }
}
close(IN);

