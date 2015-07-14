#! /usr/local/bin/perl

use strict;

my $input_mut = $ARGV[0];
my $ref_list = $ARGV[1];
my $TH_BASE = $ARGV[2];
my $TH_MAP =$ARGV[3];
my $OUTPUT_DIR =$ARGV[4];
my $PATH_TO_SAMTOOLS =$ARGV[5];
my $TEMP_BED = $OUTPUT_DIR . "/temp.base.bed";
my $TEMP_FILE = $OUTPUT_DIR . "/temp.base.txt";

my $filterQuals = "";
for (my $i = 33; $i < 33 + $TH_BASE; $i++) {
    $filterQuals .= chr($i);
}

my %pos2var = ();
open(TEMP_BED, ">" . $TEMP_BED) || die "cannot open $!";
open(IN, $input_mut) || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    my @F = split("\t", $_);
    $pos2var{$F[0] . "\t" . $F[1]} = $F[3];
    print TEMP_BED $F[0] . "\t" . ($F[1] - 1) . "\t" . $F[1] . "\n";
}
close(TEMP_BED);

my @refList = ();
open(IN, $ref_list) || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    push @refList, $_;
}
close(IN);

my %key2refInfo = ();
for (my $i = 0; $i <= $#refList; $i++) {
    # my $ret = system($PATH_TO_SAMTOOLS . "/samtools mpileup -A -B -d10000000 -q " . $TH_MAP . " -l " . $TEMP_BED . " -Q 0 " . $refList[$i] . " >> " . $TEMP_FILE . "." . $i);
    my $ret = system($PATH_TO_SAMTOOLS . "/samtools mpileup -B -d10000000 -q " . $TH_MAP . " -l " . $TEMP_BED . " -Q 0 " . $refList[$i] . " >> " . $TEMP_FILE . "." . $i);
    if ($ret != 0) {
      print STDERR "ERROR CODE: ".$ret."\n";
      exit 1
    }
    open(IN, $TEMP_FILE . "." . $i) || die "cannot open $!";
    while(<IN>) {
        s/[\r\n]//g;
        my @F = split("\t", $_);
        my $var = $pos2var{$F[0] . "\t" . $F[1]};
        my $key = $F[0] . "\t" . $F[1] . "\t" . $var;
        if (not exists $key2refInfo{$key}) {
            $key2refInfo{$key} = [ ("0,0,0,0") x ($#refList + 1) ];
        }
        $key2refInfo{$key}->[$i] = &getBaseInfo($_, $var);
    }
}


open(IN, $input_mut) || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    my @F = split("\t", $_);

    $DB::single = 1;
    my $key = $F[0] . "\t" . $F[1] . "\t" . $F[3];
    my $refInfo = exists $key2refInfo{$key} ? join("\t", @{$key2refInfo{$key}}) : join("\t", ("0,0,0,0") x ($#refList + 1)); 
    print join("\t", @F) . "\t" . $refInfo . "\n";
}
close(IN);


sub getBaseInfo {

    my @curRow = split("\t", $_[0]);

    my %insertion_p = ();
    my %insertion_n = ();
    my %deletion_p = ();
    my %deletion_n = ();
 
    my $depth = $curRow[3];
    
    if ($depth > 0) {

        # remove insertion and deletion
        while($curRow[4] =~ m/\+([0-9]+)/g) {

            my $num = $1;
            my $site = pos $curRow[4];
 
            my $var = substr($curRow[4], $site, $num);

            my $strand = "1";
            if ($var =~ tr/[acgtn]/[ACGTN]/) {
                 $strand = "-1";
            }

            my $key = join("\t", @curRow[0 .. 3]) . "\t" . $var;

            if (exists $insertion_p{$key} or exists $insertion_n{$key}) {

                if ($strand == "1") {
                    $insertion_p{$key} = $insertion_p{$key} + 1;
                } else {
                    $insertion_n{$key} = $insertion_n{$key} + 1;
                }

            } else {

                if ($strand == "1") {
                    $insertion_p{$key} = 1;
                    $insertion_n{$key} = 0;
                } else {
                    $insertion_p{$key} = 0;
                    $insertion_n{$key} = 1;
                }
            }


            substr($curRow[4], $site - length($num) - 1, $num + length($num) + 1, "");
        }


        while($curRow[4] =~ m/\-([0-9]+)/g) {

            my $num = $1;
            my $site = pos($curRow[4]);
            my $var = substr($curRow[4], $site, $num);

            my $strand = "1";
            if ($var =~ tr/[acgtn]/[ACGTN]/) {
                $strand = "-1";
            }

            my $key = join("\t", @curRow[0 .. 3]) . "\t" . $var;
            if (exists $deletion_p{$key} or exists $deletion_n{$key}) {

                if ($strand == "1") {
                    $deletion_p{$key} = $deletion_p{$key} + 1;
                } else {
                    $deletion_n{$key} = $deletion_n{$key} + 1;
                }

            } else {

                if ($strand == "1") {
                    $deletion_p{$key} = 1;
                    $deletion_n{$key} = 0;
                } else {
                    $deletion_p{$key} = 0;
                    $deletion_n{$key} = 1;
                }
            }

            substr($curRow[4], $site - length($num) - 1, $num + length($num) + 1, "");
        }

        # remove start marks!
        $curRow[4] =~ s/\^.//g;
 
        # remove end marks!
        $curRow[4] =~ s/\$//g;
 

        # for debugging !!
        if (length($curRow[4]) != length($curRow[5])) {
            print "something is wrong!!!\n";
            print join("\t", @curRow) . "\n";
            # print length($curRow[4]) . "\t" . length($curRow[5]) . "\n";
        }

        my %base2qual = ();
        $base2qual{"A"} = "";
        $base2qual{"C"} = "";
        $base2qual{"G"} = "";
        $base2qual{"T"} = "";
        $base2qual{"N"} = "";
        $base2qual{"a"} = "";
        $base2qual{"c"} = "";
        $base2qual{"g"} = "";
        $base2qual{"t"} = "";
        $base2qual{"n"} = "";

        my @bases = split("", $curRow[4]);
        my @quals = split("", $curRow[5]);
        for (my $ii = 0; $ii <= $#bases; $ii++) {
            $base2qual{$bases[$ii]} .= $quals[$ii];
        }

        # my $depth_p = length($base2qual{"A"}) + length($base2qual{"C"}) + length($base2qual{"G"}) + length($base2qual{"T"}) + length($base2qual{"N"});
        # my $depth_n = length($base2qual{"a"}) + length($base2qual{"c"}) + length($base2qual{"g"}) + length($base2qual{"t"}) + length($base2qual{"n"});

        # print join("\t", @curRow[0 .. 3]) . "\t";
        # print $depth_p . "\t" . $depth_n . "\n";

        if ($filterQuals ne "") {    
            foreach my $base (keys %base2qual) {
                $base2qual{$base} =~ s/[$filterQuals]//g;
            }
        }

        my $depth_p = length($base2qual{"A"}) + length($base2qual{"C"}) + length($base2qual{"G"}) + length($base2qual{"T"}) + length($base2qual{"N"});
        my $depth_n = length($base2qual{"a"}) + length($base2qual{"c"}) + length($base2qual{"g"}) + length($base2qual{"t"}) + length($base2qual{"n"});

        # print join("\t", @curRow[0 .. 3]) . "\t";

        my $var_p = uc($_[1]);
        my $var_n = lc($_[1]);
        return $depth_p . "," . length($base2qual{$var_p}) . "," . $depth_n . "," . length($base2qual{$var_n});

        # print length($base2qual{"A"}) . "\t" . length($base2qual{"a"}) . "\t" . length($base2qual{"C"}) . "\t" . length($base2qual{"c"}) . "\t";
        # print length($base2qual{"G"}) . "\t" . length($base2qual{"g"}) . "\t" . length($base2qual{"T"}) . "\t" . length($base2qual{"t"}) . "\n";

        # foreach my $ins (keys %insertion_p) {
        #     print $ins . "\t" . $insertion_p{$ins} . "\t" . $insertion_n{$ins} . "\n";
        # }

        # foreach my $del (keys %deletion_p) {
        #     print $del . "\t" . $deletion_p{$del} . "\t" . $deletion_n{$del} . "\n";
        # }

    }
    else {
        return "0,0,0,0";
    }

}


