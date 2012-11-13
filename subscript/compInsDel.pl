use strict;
use List::Util qw(max);

my $INPUT_var_tum = $ARGV[0];
my $INPUT_var_nor = $ARGV[1];
my $INPUT_depth_tum = $ARGV[2];
my $INPUT_depth_nor = $ARGV[3];
my $MIN_TUMOR_DEPTH = $ARGV[4];
my $MIN_NORMAL_DEPTH = $ARGV[5];
my $MIN_TUMOR_VARIANT_READ = $ARGV[6];
my $MIN_TUMOR_ALLELE_FREQ = $ARGV[7];
my $MAX_NORMAL_ALLELE_FREQ = $ARGV[8];


open(IN, $INPUT_var_tum) || die "cannot open $!";
my %ID2vars = ();
while(<IN>) {

	s/[\r\n]//g;
	my @F = split("\t", $_);
    my $key = join("\t", @F[0 .. 2]);
    my $var = $F[4];
    my $depth = $F[3]; 
     
	if ($depth >= $MIN_TUMOR_DEPTH and ($F[5] + $F[6]) / $depth >= $MIN_TUMOR_ALLELE_FREQ and $F[5] + $F[6] >= $MIN_TUMOR_VARIANT_READ) {
        $ID2vars{$key}->{$var} =  [(0) x 8];
        $ID2vars{$key}->{$var}->[1] = $F[5];
        $ID2vars{$key}->{$var}->[3] = $F[6];
   	}

}
close(IN);


open(IN, $INPUT_depth_tum) || die "cannot open $!";
while(<IN>) {
    s/[\r\n]//g;
    my @F = split("\t", $_);
    my $key = join("\t", @F[0 .. 2]);
    my $depth = $F[3];

    if (exists $ID2vars{$key}) {
        foreach my $var (keys %{$ID2vars{$key}}) {
            $ID2vars{$key}->{$var}->[0] = $F[4];
            $ID2vars{$key}->{$var}->[2] = $F[5];
        }
    }

}
close(IN);


open(IN, $INPUT_depth_nor) || die "cannot open $!";
while(<IN>) {
    s/[\r\n]//g;
    my @F = split("\t", $_);
    my $key = join("\t", @F[0 .. 2]);
    my $depth = $F[4] + $F[5];

    if (exists $ID2vars{$key} and $depth >= $MIN_NORMAL_DEPTH) {
  
       foreach my $var (keys %{$ID2vars{$key}}) {
           $ID2vars{$key}->{$var}->[4] = $F[4];
           $ID2vars{$key}->{$var}->[6] = $F[5];
       }  
    } else {
        delete $ID2vars{$key};
    } 
  
}  
close(IN);


open(IN, $INPUT_var_nor) || die "cannot open $!";
while(<IN>) {

    s/[\r\n]//g;
    my @F = split("\t", $_);
    my $key = join("\t", @F[0 .. 2]);
    my $var = $F[4];
    my $depth = $F[3];
   
    if (exists $ID2vars{$key}->{$var} and ($F[5] + $F[6]) / $depth < $MAX_NORMAL_ALLELE_FREQ) {
 
        $ID2vars{$key}->{$var}->[5] = $F[5];
        $ID2vars{$key}->{$var}->[7] = $F[6];
    } else {
        delete $ID2vars{$key}->{$var};
    }


}
close(IN);

 

foreach my $key (sort chrpos keys %ID2vars) {

    foreach my $var (sort keys %{$ID2vars{$key}}) {
        my @tmpBases = @{$ID2vars{$key}->{$var}};
        if ($tmpBases[4] > 0 || $tmpBases[6] > 0) {
            print $key . "\t" . $var . "\t" . join(",", @tmpBases[0 .. 3]) . "\t" . join(",", @tmpBases[4 .. 7]) . "\n";
        }
    }
}


sub chrpos {

    my @posa = split("\t", $a);
    my @posb = split("\t", $b);

    
    if ($posa[0] > $posb[0]) {
        return 1;
    } elsif ($posa[0] < $posb[0]) {
        return -1;
    } else {
        if ($posa[1] > $posb[1]) {
            return 1;
        } else {
            return -1;
        }
    }
}

