
%base2ind = ();
$base2ind{"A"} = 0;
$base2ind{"C"} = 1;
$base2ind{"G"} = 2;
$base2ind{"T"} = 3;

my $INPUT = $ARGV[0];
my $MIN_NORMAL_DEPTH = $ARGV[1];

open(IN, $INPUT) || die "cannot open $!";
$ID2bases_tum = ();
<IN>; # skip header line
while(<IN>) {

  s/[\r\n]//g;
  @curRow = split("\t", $_);

  $ref = $curRow[2];

  # the numbers of each bases
  ##########
  $A = $curRow[4] + $curRow[5];
  $C = $curRow[6] + $curRow[7];
  $G = $curRow[8] + $curRow[9];
  $T = $curRow[10] + $curRow[11];
  $depth = $A + $C + $G + $T;

  if ($depth >= $MIN_NORMAL_DEPTH) {

    # determine which is the most frequent mismatch
    ##########
    @Nums = ($A, $C, $G, $T);
    $Nums[$base2ind{$ref}] = 0;

    $mis = "N";
    $misNum = 0;
    foreach $base (keys %base2ind) {
	
      if ($Nums[$base2ind{$base}] > $misNum) {
        $misNum = $Nums[$base2ind{$base}];
        $mis = $base;
      }
    }

    $misRate = ($misNum / $depth);
    if ($mis ne "N" and $misRate > 0.07) {

      # ratio of mismatch bases in each strand
      ##########
      $sratio = 0;
      if ($mis eq "A") {
        $sratio = $curRow[4] / ($curRow[4] + $curRow[5]);
      } elsif ($mis eq "C") {
        $sratio = $curRow[6] / ($curRow[6] + $curRow[7]);
      } elsif ($mis eq "G") {
        $sratio = $curRow[8] / ($curRow[8] + $curRow[9]);
      } elsif ($mis eq "T") {
        $sratio = $curRow[10] / ($curRow[10] + $curRow[11]);
      } else {
        print $mis . $mis. "\t" . "something is wrong!\n";
      }
            
      print $curRow[0] . "\t" . $curRow[1] . "\t" . $curRow[2] . "\t" . $mis . "\t" . $curRow[3] . "\t" . $A . "\t" . $C . "\t" . $G . "\t" . $T . "\t"  . $misRate . "\t" . $sratio . "\n";

    }
  }
}
close(IN);

