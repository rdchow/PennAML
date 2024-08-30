#!/usr/bin/perl
use strict;
use warnings;
my %hash;
my $fh;
open ($fh,"<../fig2/uniqueStage-sample-names.txt") || die "c";
while (<$fh>){
    s/[\r\n]//g;
    my $line = $_;
    my @lines = split ("\t",$line);
    $hash{$lines[0]} = 1;
}
close $fh;

my $ofh;
open ($fh,"<sample-info.p6.deid.ELNlab.txt") || die "c";
open ($ofh,">sample-info.p6.deid.ELNlab.uniqFilt.txt") || die "a";
my $head = <$fh>;
print $ofh $head;
while (<$fh>){
    s/[\r\n]//g;
    my $line = $_;
    my @lines = split ("\t",$line);
    if (defined $hash{$lines[0]}){
        print $ofh $line,"\n";
    }
}
close $fh;
close $ofh;

open ($fh,"<AML_data.muts.chroms.txt") || die "c";
open ($ofh,">AML_data.muts.chroms.uniqFilt.txt") || die "a";
while (<$fh>){
    s/[\r\n]//g;
    my $line = $_;
    my @lines = split ("\t",$line);
    if (defined $hash{$lines[0]}){
        print $ofh $line,"\n";
    }
}
close $fh;
close $ofh;