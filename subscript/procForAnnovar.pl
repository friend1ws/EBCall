# conver three candidate variation files (mutation, insertion and deletion) for process via Annovar 
use strict;

my $mutationFile = $ARGV[0];
my $insersionFile = $ARGV[1];
my $deletionFile = $ARGV[2];

print "Chr" . "\t" . "Start" . "\t" . "End" . "\t" . "Ref" . "\t" . "Obs" . "\t" . "bases_tumor" . "\t" . "bases_normal" . "\t" . "misRate_tumor" . "\t" . "strandRatio_tumor" . "\t" . "depth_tumor" . "\t" . "variantNum_tumor" . "\t" . "misRate_normal" . "\t" .  "strandRatio_normal" . "\t" . "depth_normal" . "\t" . "variantNum_normal" . "\t" . "p-value" . "\t" . "p-value (+strand)" . "\t" . "p-value (-strand)" . "\t" . "p-value(Fisher)" . "\t" .  "alpha (+starnd)" . "\t" . "beta (+strand)" . "\t" . "alpha (-strand)" . "\t" . "beta(-strand)" . "\n";

# mutation
if (-e $mutationFile) { 

open(IN, $mutationFile) || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    my @F = split("\t", $_);
    my $chr = $F[0];
    my $pos = $F[1];
    print $chr . "\t" . $pos . "\t" . $pos . "\t" . join("\t", @F[2 .. $#F]) . "\n";
}
close(IN);

}


# insertion
if (-e $insersionFile) { 

open(IN, $insersionFile) || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    my @F = split("\t", $_);

    my $chr = $F[0];
    my $pos = $F[1];
    my $ref = $F[2];
    my $mis = $F[3];

    print $chr . "\t" . $pos . "\t" . $pos . "\t" . "-" . "\t" . $mis . "\t" . join("\t", @F[4 .. $#F]) . "\n";

}
close(IN);

}

# deletion 
if (-e $deletionFile) { 

open(IN, $deletionFile) || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    my @F= split("\t", $_);

    my $chr = $F[0];
    my $pos = $F[1];
    my $ref = $F[2];
    my $mis = $F[3];

    print $chr . "\t" . ($pos + 1) . "\t" . ($pos + length($mis)) . "\t" . $mis . "\t" . "-" . "\t" . join("\t", @F[4 .. $#F]) . "\n";
}
close(IN);

}    
