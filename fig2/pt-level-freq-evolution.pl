#!/usr/bin/perl
# compare variants on patient-specific level, across timepoints

use strict;
use warnings;

my %keepsamples; 
my $fh;
my $ofh;
my %ptvarhash; # patient-specific hash of variants and VAFs
my %check;
my %ptsamples; # hash with array detailing if a patient has a given sample type or not

# first pass we are keeping a catalog of all variants found in each patient
open ($ofh,">pt-level-mut-freqs.txt") || die "a";
open ($fh,"<../preprocess/AML_data.muts.chroms.txt") || die "c";
while (<$fh>){
    s/[\r\n]//g;
    my $line = $_;
    my @lines = split ("\t",$line);
    my $sample = $lines[0];
    my $gene = $lines[1];
    my $var = $lines[2];
    my $patient = $lines[6];
    my $status = "other";

    # only looking at point mutations/indels (or empty, eg no muts)
    if ($lines[4] eq "MISSENSE" || $lines[4] eq "TRUNC" || $lines[4] eq "INFRAME"|| $lines[4] eq "OTHER" || $lines[4] eq ""){
        if ($lines[5] =~ /^dn/){
            $status = "dn";
        }
        elsif ($lines[5] =~ /^CR1/){
            $status = "CR1";
        }
        elsif ($lines[5] =~ /^CR2/){
            $status = "CR2";
        }
        elsif ($lines[5] =~ /^rel1/){
            $status = "REL1";
        }
        elsif ($lines[5] =~ /^ref1/ || $lines[5] =~ /primref/ || $lines[5] eq "ref.AML" || $lines[5] eq "refc.AML"){
            $status = "REF1";
        }
        elsif ($lines[5] =~ /^rel2/){
            $status = "REL2";
        }
        elsif ($lines[5] =~ /^ref2/){
            $status = "REF2";
        }
        
        my $test = $sample.":".$gene;
        my $ptstatus = $patient.":".$status;
        if (!defined $ptsamples{$patient}){
            $ptsamples{$patient} = [0,0,0,0,0,0,0];
        }

        if ($status eq "dn") {$ptsamples{$patient}[0]++;}
        if ($status eq "CR1") {$ptsamples{$patient}[1]++;}
        if ($status eq "REL1") {$ptsamples{$patient}[2]++;}
        if ($status eq "REF1") {$ptsamples{$patient}[3]++;}
        if ($status eq "CR2") {$ptsamples{$patient}[4]++;}
        if ($status eq "REL2") {$ptsamples{$patient}[5]++;}
        if ($status eq "REF2") {$ptsamples{$patient}[6]++;}

        # only count the earliest sample for a given disease stage for each patient
        # if this is the first time seeing this combination of patient:stage, then save that sample to %keepsamples
        if (!defined $check{$ptstatus}){
            $check{$ptstatus} = 1;
            $keepsamples{$sample} = "$patient\t$status";
        }

        # if this sample is in the "pass list" from above
        # keep a hash of all variants identified in a given patient
        if (defined $keepsamples{$sample}){
            # if this is the first time seeing a given gene variant in this patient, then initialize to all NA
            if (!defined $ptvarhash{"$patient.$gene.$var"} && $gene ne ""){
                $ptvarhash{"$patient.$gene.$var"} = ["NA","NA","NA","NA","NA","NA","NA"]; # diagnosis, CR1, rel1, ref1, CR2, rel2, ref2
            }
        }
    }       
}
close $fh;

#foreach my $k (sort keys %ptsamples){
#    print "$k\t$ptsamples{$k}[0]\t$ptsamples{$k}[1]\n";
#}

open ($fh,"<../preprocess/AML_data.muts.chroms.txt") || die "c";
while (<$fh>){
    s/[\r\n]//g;
    my $line = $_;
    my @lines = split ("\t",$line);
    my $sample = $lines[0];
    my $gene = $lines[1];
    my $var = $lines[2];
    my $faf = $lines[3];
    my $patient = $lines[6];
    my $status = "other";

    # only looking at point mutations/indels
    # only look at samples "passing" from above
    if (defined $keepsamples{$sample} && ($lines[4] eq "MISSENSE" || $lines[4] eq "TRUNC" || $lines[4] eq "INFRAME"|| $lines[4] eq "OTHER")){ 
        if ($lines[5] =~ /^dn/){
            $status = "dn";
        }
        elsif ($lines[5] =~ /^CR1/){
            $status = "CR1";
        }
        elsif ($lines[5] =~ /^CR2/){
            $status = "CR2";
        }
        elsif ($lines[5] =~ /^rel1/){
            $status = "REL1";
        }
        elsif ($lines[5] =~ /^ref1/ || $lines[5] =~ /primref/ || $lines[5] eq "ref.AML" || $lines[5] eq "refc.AML"){
            $status = "REF1";
        }
        elsif ($lines[5] =~ /^rel2/){
            $status = "REL2";
        }
        elsif ($lines[5] =~ /^ref2/){
            $status = "REF2";
        }

        if ($lines[4] ne ""){
            if ($status eq "dn") {$ptvarhash{"$patient.$gene.$var"}[0] = $faf;}
            if ($status eq "CR1") {$ptvarhash{"$patient.$gene.$var"}[1] = $faf;}
            if ($status eq "REL1") {$ptvarhash{"$patient.$gene.$var"}[2] = $faf;}
            if ($status eq "REF1") {$ptvarhash{"$patient.$gene.$var"}[3] = $faf;}
            if ($status eq "CR2") {$ptvarhash{"$patient.$gene.$var"}[4] = $faf;}
            if ($status eq "REL2") {$ptvarhash{"$patient.$gene.$var"}[5] = $faf;}
            if ($status eq "REF2") {$ptvarhash{"$patient.$gene.$var"}[6] = $faf;}
        }
        else {
            if ($status eq "dn") {$ptvarhash{"$patient.$gene.$var"}[0] = "NA";}
            if ($status eq "CR1") {$ptvarhash{"$patient.$gene.$var"}[1] = "NA";}
            if ($status eq "REL1") {$ptvarhash{"$patient.$gene.$var"}[2] = "NA";}
            if ($status eq "REF1") {$ptvarhash{"$patient.$gene.$var"}[3] = "NA";}
            if ($status eq "CR2") {$ptvarhash{"$patient.$gene.$var"}[4] = "NA";}
            if ($status eq "REL2") {$ptvarhash{"$patient.$gene.$var"}[5] = "NA";}
            if ($status eq "REF2") {$ptvarhash{"$patient.$gene.$var"}[6] = "NA";}
        }
    }
}
close $fh;

# change NA to 0 if the patient-status combo was sequenced
foreach my $k (sort {$a cmp $b} keys %ptvarhash){
    my @kk = split("\Q.",$k);
    my $pt= $kk[0];
    my $gene = $kk[1];
    my $var = $kk[2];
    my @dat = @{$ptvarhash{$k}};

    # if the sample exists for that timepoint in a given patient and otherwise did not have a variant VAF found, initialize the variant VAF to 0
    for (my $i = 0; $i < scalar@dat; $i++){
        if ($ptsamples{$pt}[$i] > 0 && $dat[$i] eq "NA"){ 
            $dat[$i] = 0;
        }
    }
    @{$ptvarhash{$k}} = @dat;
}     


print $ofh "patientID\tgene\tVariant\tdiagnosis\tCR1\tREL1\tREF1\tCR2\tREL2\tREF2\n";
# print the output
foreach my $k (sort {$a cmp $b} keys %ptvarhash){
    my @kk = split("\Q.",$k);
    print $ofh "$kk[0]\t$kk[1]\t$kk[2]\t";
    my @dat = @{$ptvarhash{$k}};

    for (my $i =0; $i < scalar @dat -1;$i++){
        print $ofh "$dat[$i]\t";
    }
    print $ofh $dat[scalar@dat-1],"\n";
}
close $ofh;

open ($ofh,">uniqueStage-sample-names.txt") || die "c";
foreach my $k (sort keys %keepsamples){
    print $ofh "$k\t$keepsamples{$k}\n";
}
close $ofh;