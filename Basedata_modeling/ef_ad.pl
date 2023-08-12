
use warnings;
use autodie;
use strict;

open my $filename1, "<", "efp_ct.txt";
my %hash;
while ( my $pos = <$filename1> ) {
        $hash{$pos}=' ';
        }
close $filename1;
open my $filename2, '<', '../ambRm_chi.assoc' ;
open my $out, '>', '../ef_ambRm_chi.assoc';
my $line;
while ( $line = <$filename2> ) {
        if ( $line !~ /CHR/ ) {
                my @e = split(" ", $line);
                my $bpos = $e[2]."\n";
                if (exists ($hash{$bpos})) {
                        print "$e[0]\t$e[1] corrected \n ";
                        print $out "$e[0]\t$e[1]\t$e[2]\t$e[6]\t$e[4]\t$e[5]\t$e[3]\t$e[7]\t$e[8]\t$e[9]\n";
                }else{ print $out "$e[0]\t$e[1]\t$e[2]\t$e[3]\t$e[4]\t$e[5]\t$e[6]\t$e[7]\t$e[8]\t$e[9]\n";
        }
        }else{print $out "CHR\tSNP\tBP\tA1\tF_A\tF_U\tA2\tCHISQ\tP\tOR\n";}
}
close $filename2;
close $out;
