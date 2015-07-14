#! /usr/local/bin/perl

use strict;

my $input_mut = $ARGV[0];
my $ref_list = $ARGV[1];
my $insdel = $ARGV[2]; # 1: insertion, 2: deletion
my $TH_MAP = $ARGV[3];
my $OUTPUT_DIR =$ARGV[4];
my $PATH_TO_SAMTOOLS =$ARGV[5];

my $TEMP_BED = $OUTPUT_DIR;
my $TEMP_FILE = $OUTPUT_DIR;

if ($insdel == 1) {
    $TEMP_BED = $TEMP_BED . "/temp.ins.bed";
    $TEMP_FILE = $TEMP_FILE . "/temp.ins.txt";
} elsif ($insdel == 2) {
    $TEMP_BED = $TEMP_BED . "/temp.del.bed";
    $TEMP_FILE = $TEMP_FILE . "/temp.del.txt";
} else {
    die "the index for insertion or deletion is wrongly specified.\n";
}

my %pos2var = ();
open(TEMP_BED, ">" . $TEMP_BED) || die "cannot open $!";
open(IN, $input_mut) || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    my @F = split("\t", $_);

    my $var = $F[3];
    if ($insdel == 1) {
        $var = "+" . $var;
    } elsif ($insdel == 2) {
        $var = "-" . $var;
    } else {
        die "the index for insertion or deletion is wrongly specified.\n";
    }
    if (exists $pos2var{$F[0] . "\t" . $F[1]}) {
        $pos2var{$F[0] . "\t" . $F[1]} = $pos2var{$F[0] . "\t" . $F[1]} . "\t" . $var;
    } else {
        $pos2var{$F[0] . "\t" . $F[1]} = $var;
    }

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
        my $vars = $pos2var{$F[0] . "\t" . $F[1]};
        if ($vars eq "") {
            print "At the key:\n";
            print join("\t", @F) . "\n";
            die "No variant extracted. Something is wrong.\n";
        }

        foreach my $var (split("\t", $vars)) {

            my $key = $F[0] . "\t" . $F[1] . "\t" . $var;
            if (not exists $key2refInfo{$key}) {
                $key2refInfo{$key} = [ ("0,0,0,0") x ($#refList + 1) ];
            }
            $key2refInfo{$key}->[$i] = &getInDelInfo($_, $var);
            # print $key . "\t" . $key2refInfo{$key}->[$i] . "\n";
        }
    }
}

$DB::single = 1;

open(IN, $input_mut) || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    my @F = split("\t", $_);

    my $var = $F[3];
    if ($insdel == 1) {
        $var = "+" . $var;
    } elsif ($insdel == 2) {
        $var = "-" . $var;
    } else {
        die "the index for insertion or deletion is wrongly specified.\n";
    }

    # print join("\t", @F) . "\n";
    my $key = $F[0] . "\t" . $F[1] . "\t" . $var;
    my $refInfo = exists $key2refInfo{$key} ? join("\t", @{$key2refInfo{$key}}) : join("\t", ("0,0,0,0") x ($#refList + 1));
    print join("\t", @F) . "\t" . $refInfo . "\n";
}
close(IN);




sub getInDelInfo {

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

            # my $key = join("\t", @curRow[0 .. 3]) . "\t" . $var;
            my $key = "+" . $var;

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

            # my $key = join("\t", @curRow[0 .. 3]) . "\t" . $var;
            my $key = "-" . $var;
            
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
            # print length($curRow[4]) . "\t" . length($curRow[5]) . "\n";
            print $curRow[0] ."\t". $curRow[1] ."\t". $curRow[2] ."\t". $curRow[3] ."\t". length($curRow[4]) . "\t" . length($curRow[5]) . "\n";
            exit 1;
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

        my $depth_p = length($base2qual{"A"}) + length($base2qual{"C"}) + length($base2qual{"G"}) + length($base2qual{"T"}) + length($base2qual{"N"});
        my $depth_n = length($base2qual{"a"}) + length($base2qual{"c"}) + length($base2qual{"g"}) + length($base2qual{"t"}) + length($base2qual{"n"});

        my $misNum_p = 0;  
        my $misNum_n = 0;      
        if ($insdel == 1) {
            # $misNum_p = exists $insertion_p{$_[1]} ? $insertion_p{$_[1]} : 0;
            # $misNum_n = exists $insertion_n{$_[1]} ? $insertion_n{$_[1]} : 0;
        
           # conservative evaluation
           foreach my $key (keys %insertion_p) {
               $misNum_p = $misNum_p + $insertion_p{$key};
           }
           foreach my $key (keys %insertion_n) {
               $misNum_n = $misNum_n + $insertion_n{$key};
           }
        }
        elsif ($insdel == 2) {
            # $misNum_p = exists $deletion_p{$_[1]} ? $deletion_p{$_[1]} : 0;
            # $misNum_n = exists $deletion_n{$_[1]} ? $deletion_n{$_[1]} : 0;
           
            # conservative evaluation
            foreach my $key (keys %deletion_p) {
                $misNum_p = $misNum_p + $deletion_p{$key};
            }
            foreach my $key (keys %deletion_n) {
                $misNum_n = $misNum_n + $deletion_n{$key};
            }
        }

        return $depth_p . "," . $misNum_p . "," . $depth_n . "," . $misNum_n;
    }
    else {
        return "0,0,0,0";
    }

}


