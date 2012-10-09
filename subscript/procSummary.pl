open(IN, $ARGV[0]) || die "cannot open $!";
$_ = <IN>;
s/[\r\n]//g;
@curRow = split("\t", $_);

print join("\t", @curRow[0 .. ($#curRow - 1)]). "\t" . "bases_tumor" . "\t" . "bases_normal" . "\t" . "misRate_tumor" . "\t" . "strandRatio_tumor" . "\t" . "depth_tumor" . "\t" . "variantNum_tumor" . "\t" . "misRate_normal" . "\t" .  "strandRatio_normal" . "\t" . "depth_normal" . "\t" . "variantNum_normal" . "\t" . "p-value" . "\t" . "p-value (+strand)" . "\t" . "p-value (-strand)" . "\t" . "p-value(Fisher)" . "\t" .  "alpha (+starnd)" . "\t" . "beta (+strand)" . "\t" . "alpha (-strand)" . "\t" . "beta(-strand)" . "\n"; 
while(<IN>) {
    print $_;
}
close(IN);

