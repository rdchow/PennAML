#!/usr/bin/perl
# calculate mutation frequency among samples at diagnosis
use strict;
use warnings;
my $fh;
# already filtered to only include one sample per disease stage
open ($fh,"<../preprocess/AML_data.muts.chroms.uniqFilt.dropDiscordantPanelGenes.txt") || die "c";
my $ofh;
open ($ofh,">mut-frequencies.txt") || die "c";
my %samples; # tracks total number of samples in each cateogry
my $diagct = 0;
my $cr1ct = 0;
my $cr2ct = 0;
my $rel1ct = 0;
my $ref1ct = 0;
my $rel2ct = 0;
my $ref2ct = 0;
my %vars;
my %check; # tracks for unique pt-stage and sample-gene combos
my %samplecheck2;

my %keepsamples;
while (<$fh>){
    s/[\r\n]//g;
    my @lines = split("\t",$_);
    my $sample = $lines[0];
    my $gene = $lines[1];

    my $patient = $lines[6];
    my $status = "other";
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

    # initialize hash that counts # of mut samples in each category
    if (!defined $vars{$gene}){ 
        $vars{$gene} = [0,0,0,0,0,0,0]; # diagnosis, CR1, rel1, ref1, CR2, rel2, ref2
    }

    # only count the earliest sample for a given disease stage for each patient
    # if this is the first time seeing this combination of patient:stage, then save that sample to %keepsamples
    if (!defined $check{$ptstatus}){
        $check{$ptstatus} = 1;
        $keepsamples{$sample} = "$patient\t$status";
    }      

    # if this sample is in the "pass list" from above
    if (defined $keepsamples{$sample}){
        if (!defined $check{$test}){ # if this gene has not yet been counted for this sample
            $check{$test} = 1;
            
            # diagnosis, CR1, rel1, ref1, CR2, rel2, ref2
            if ($status eq "dn"){ # diagnosis samples
                $vars{$gene}[0]++;
            }
            elsif ($status eq "CR1"){ # CR1 samples 
                $vars{$gene}[1]++;
            }
            elsif ($status eq "REL1") { # rel1 samples 
                $vars{$gene}[2]++;
            }
            elsif ($status eq "REF1") { # ref1 samples 
                $vars{$gene}[3]++;
            }
            elsif ($status eq "CR2"){ # CR2+ samples 
                $vars{$gene}[4]++;
            }
            elsif ($status eq "REL2"){ # rel2 samples 
                $vars{$gene}[5]++;
            }
            elsif ($status eq "REF2"){ # ref2 samples 
                $vars{$gene}[6]++;
            }
        }

        if (!defined $samplecheck2{$sample}){
            $samplecheck2{$sample} = 1;
            # count samples
            if ($status eq "dn"){ # diagnosis samples
                $diagct++;
            }
            elsif ($status eq "CR1"){ # CR1 samples 
                $cr1ct++;
            }
            elsif ($status eq "REL1") { # rel1 samples 
                $rel1ct++;
            }
            elsif ($status eq "REF1") { # rel1 samples 
                $ref1ct++;
            }
            elsif ($status eq "CR2"){ # CR2 samples 
                $cr2ct++;
            }
            elsif ($status eq "REL2"){ # rel2 samples 
                $rel2ct++;
            }
            elsif ($status eq "REF2"){ # ref2 samples 
                $ref2ct++;
            }
        }
    }
}
close $fh;  
print $ofh "$diagct diagnosis samples, $cr1ct CR1 samples, $rel1ct Rel1 samples, $ref1ct Ref1 samples, ";
print $ofh "$cr2ct CR2 samples, $rel2ct Rel2 samples, $ref2ct Ref2 samples\n";
print $ofh "gene\tDiagnosisMutCt\tDiagnosisMutFreq\t";
print $ofh "CR1mutCt\tCR1mutFreq\tRel1MutCt\tRel1MutFreq\tRef1MutCt\tRef1MutFreq\t";
print $ofh "CR2mutCt\tCR2mutFreq\tRel2MutCt\tRel2MutFreq\tRef2MutCt\tRef2MutFreq\t";
print $ofh "CR1_2mutCt\tCR1_2mutFreq\tRel1_2MutCt\tRel1_2MutFreq\tRef1_2MutCt\tRef1_2MutFreq\n";

foreach my $k (sort {$vars{$b}[0] <=> $vars{$a}[0]} keys %vars){
    if ($k ne ""){
        print $ofh "$k\t$vars{$k}[0]\t",$vars{$k}[0]/$diagct,"\t";
        print $ofh "$vars{$k}[1]\t",$vars{$k}[1]/$cr1ct,"\t";
        print $ofh "$vars{$k}[2]\t",$vars{$k}[2]/$rel1ct,"\t";
        print $ofh "$vars{$k}[3]\t",$vars{$k}[3]/$ref1ct,"\t";
        print $ofh "$vars{$k}[4]\t",$vars{$k}[4]/$cr2ct,"\t";
        print $ofh "$vars{$k}[5]\t",$vars{$k}[5]/$rel2ct,"\t";
        print $ofh "$vars{$k}[6]\t",$vars{$k}[6]/$ref2ct,"\t";

        print $ofh ($vars{$k}[1] + $vars{$k}[4]),"\t",($vars{$k}[1] + $vars{$k}[4])/($cr1ct+$cr2ct),"\t";
        print $ofh ($vars{$k}[2] + $vars{$k}[5]),"\t",($vars{$k}[2] + $vars{$k}[5])/($rel1ct+$rel2ct),"\t";
        print $ofh ($vars{$k}[3] + $vars{$k}[6]),"\t",($vars{$k}[3] + $vars{$k}[6])/($ref1ct+$ref2ct),"\n";
    }
}
close $ofh;

