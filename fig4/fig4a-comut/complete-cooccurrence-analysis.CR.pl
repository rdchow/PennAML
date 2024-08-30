#!/usr/bin/perl
# essentially we reorder the genes so that the alphabetically first gene is "A"
use strict;
use warnings;
my $fh;
my $ofh;

open ($fh,"<mutual-exclusivity-table.diagnosis.tsv") || die "c";
open ($ofh,">mutual-exclusivity-table.diagnosis.2.tsv") || die "c";
my $head = <$fh>;
print $ofh $head;
my %hash;
while (<$fh>){
    s/[\r\n]//g;
    my $line = $_;
    my @lines = split ("\t",$line);
    my $gene1 = $lines[0];
    my $gene2 = $lines[1];

    if ((uc($gene1) cmp uc($gene2)) == -1){ 
        print $ofh $line,"\n";
        my $combo = $gene1.":".$gene2;
        $hash{$combo} = 1;
    }
    elsif ((uc($gene1) cmp uc($gene2)) == 1){
        print $ofh "$gene2\t$gene1\t$lines[2]\t$lines[4]\t$lines[3]";
        for (my $i = 5; $i < scalar @lines; $i++){
            print $ofh "\t$lines[$i]";
        }
        print $ofh "\n";
        my $combo = $gene2.":".$gene1;
        $hash{$combo} = 1;
    }
}
close $fh;

open ($fh,"<mutual-exclusivity-table.CR1.tsv") || die "c";
open ($ofh,">mutual-exclusivity-table.CR1.2.tsv") || die "c";
my $head2 = <$fh>;
print $ofh $head2;
my %hash2;
while (<$fh>){
    s/[\r\n]//g;
    my $line = $_;
    my @lines = split ("\t",$line);
    my $gene1 = $lines[0];
    my $gene2 = $lines[1];

    if ((uc($gene1) cmp uc($gene2)) == -1){ 
        print $ofh $line,"\n";
        my $combo = $gene1.":".$gene2;
        $hash2{$combo} = 1;
    }
    elsif ((uc($gene1) cmp uc($gene2)) == 1){
        print $ofh "$gene2\t$gene1\t$lines[2]\t$lines[4]\t$lines[3]";
        for (my $i = 5; $i < scalar @lines; $i++){
            print $ofh "\t$lines[$i]";
        }
        print $ofh "\n";
        my $combo = $gene2.":".$gene1;
        $hash2{$combo} = 1;
    }
}
close $fh;

# add any gene pairs that were covered in the "diagnosis" section but not in "CR1"
foreach my $k (keys %hash){
    if (!defined $hash2{$k}){
        my @genes = split(":",$k);
print $genes[0],":",$genes[1],"\n";
        print $ofh "$genes[0]\t$genes[1]\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
    }
}
close $ofh;


