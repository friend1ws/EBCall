#! /usr/local/bin/perl

use strict;

my $input = $ARGV[0];
my $input_ins = $ARGV[1];
my $input_del = $ARGV[2];

my $thres_af = 0.1;
my $thres_depth = 10;
my $neighbor = 5;

my %pos2filt = ();
open(IN, $input_ins) || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    my @F = split("\t", $_);
   
    if ($F[5] > $thres_depth and $F[7] > $thres_af) { 

        for (my $i = $F[1] - $neighbor; $i <= $F[2] + $neighbor; $i++) {
            $pos2filt{$F[0] . "\t" . $i} = 1;
        }
    }

}
close(IN);


open(IN, $input_del) || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    my @F = split("\t", $_);

    if ($F[5] > $thres_depth and $F[7] > $thres_af) {
        for (my $i = $F[1] - $neighbor; $i <= $F[2] + $neighbor; $i++) {
            $pos2filt{$F[0] . "\t" . $i} = 1;
        }
    }
}   
close(IN);
    
    
open(IN, $input) || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    my @F = split("\t", $_);

    my $flag = 0;
    for (my $i = $F[1]; $i <= $F[2]; $i++) {
        if (exists $pos2filt{$F[0] . "\t" . $i}) {
            $flag = 1;
        }
    }

    if ($flag == 0) {
        print join("\t", @F) . "\n";
    }
}
close(IN);


