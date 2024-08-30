#!/usr/bin/perl
use strict;
use warnings;
my $fh;

# reformat mut freq table for CALDER analysis
# specificially for pt 132: diagnosis, REL1, REF1

open ($fh,"<../../fig2/pt-level-mut-freqs.txt") || die "c";
my %patients;
<$fh>;
while (<$fh>){
    s/[\r\n]//g;
    my $line = $_;
    my @lines = split("\t",$line);
    # specifically take pt 132
    if ($lines[0] eq "132"){
        $patients{$lines[0]} = 1;
    }
}
close $fh;

foreach my $k (keys %patients){
    my $ofh;
    open ($fh,"<../../fig2/pt-level-mut-freqs.txt") || die "c";
    open ($ofh,">./input-matrices/pt.$k.txt") || die "a";
    <$fh>;
    my %hash; # patient-specific hash of the variants found at diagnosis, CR1 and/or REL1
    my %cr1check; # hash to determine whether to print the CR1 row
    while (<$fh>){
        s/[\r\n]//g;
        my $line = $_;
        my @lines = split("\t",$line);
        my $variant = $lines[1].":".$lines[2];

        # if this variant is in the current patient under examination at diag or REL1/REF1
        if ($lines[0] eq $k && ($lines[3] !=0 || $lines[6] != 0 || $lines[5] != 0)){
            # if this patient was seuqenced at CR1, save to the cr1check hash
            if ($lines[4] ne "NA"){
                $cr1check{$k} = 1;
            }

            $hash{$variant} = [$lines[3],$lines[4],$lines[5],$lines[6]]; 
            # Diagnosis, CR1, REL1, REF1  VAFs

            
        }
    }
    close $fh;

    # print the calder-formated matrices
    
    foreach my $j (sort keys %hash){
        print $ofh "\t$j\t$j";
    }
    print $ofh "\n";

    # print the diagnosis VAFs
    print $ofh "diagnosis";
    foreach my $j (sort keys %hash){
        #if ($hash{$j}[0] > 0){
            print $ofh "\t",10000-($hash{$j}[0]*100),"\t",$hash{$j}[0]*100;
        #}
        #else {
             #print $ofh "\t",0,"\t",0;
       #}
    }
    print $ofh "\n";

    # print the CR1 VAFs, if this patient was sequenced at CR1
    if (defined $cr1check{$k}){
        print $ofh "CR1";
        foreach my $j (sort keys %hash){
            # if ($hash{$j}[1] > 0){
                print $ofh "\t",10000-($hash{$j}[1]*100),"\t",$hash{$j}[1]*100;
            # }
            # else {
            #     print $ofh "\t",0,"\t",0;
            # }
        }
        print $ofh "\n";
    }


    # print the REL1 VAFs
    print $ofh "REL1";
    foreach my $j (sort keys %hash){
        # if ($hash{$j}[2] > 0){
            print $ofh "\t",10000-($hash{$j}[2]*100),"\t",$hash{$j}[2]*100;
        # }
        # else {
        #     print $ofh "\t",0,"\t",0;
        # }
    }
    print $ofh "\n";

     # print the REF1 VAFs
    print $ofh "REF1";
    foreach my $j (sort keys %hash){
        # if ($hash{$j}[2] > 0){
            print $ofh "\t",10000-($hash{$j}[3]*100),"\t",$hash{$j}[3]*100;
        # }
        # else {
        #     print $ofh "\t",0,"\t",0;
        # }
    }
    print $ofh "\n";
    close $ofh;
}



