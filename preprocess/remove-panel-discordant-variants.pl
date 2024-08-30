#!/usr/bin/perl
# flag diagnosis:CR1:REL1 combos with different sequencing panels that were used

# essentially, will store the NGS panel used for each sample
# Then go back through the mutation calls and delete any call that is only there because a newer panel was used

use strict;
use warnings;
my $fh;
my %samples;

my %panels; #[diag, cr1, rel1]

open ($fh,"<sample-info.p6.ELNlab.uniqFilt.txt") || die "c";
<$fh>;
while (<$fh>){
    s/[\r\n]//g;
    my $line = $_;
    my @lines = split ("\t",$line);
    $samples{$lines[0]} = [$lines[1],$lines[5],$lines[31]]; # patientID, NGS panel, disease stage

    if (!defined $panels{$lines[1]}){
        $panels{$lines[1]} = ["NA","NA","NA"];
    }
    if ($lines[31] eq "dn"){$panels{$lines[1]}[0] = $lines[5];}
    if ($lines[31] eq "CR1"){$panels{$lines[1]}[1] = $lines[5];}
    if ($lines[31] eq "REL1"){$panels{$lines[1]}[2] = $lines[5];}
}
close $fh;   

# flag patients with discordant NGS panels
my $ofh;
my %discordants;
open ($ofh,">discordant-samples.txt") || die "c";
print $ofh "patientID\tdn_panel\tCR1_panel\tREL1_panel\n";
foreach my $k (sort keys %panels){
    my @dat = @{$panels{$k}};

    if ( (($dat[0] ne $dat[1]) && ($dat[0] ne "NA" && $dat[1] ne "NA")) || (($dat[0] ne $dat[2]) && ($dat[0] ne "NA" && $dat[2] ne "NA")) || (($dat[1] ne $dat[2]) && ($dat[1] ne "NA" && $dat[2] ne "NA")) ){
        print $ofh "$k\t$dat[0]\t$dat[1]\t$dat[2]\n";
        $discordants{$k} = 1;
    }

}
close $ofh;

my %newgenes; # genes that are only in the newer NGS panel
open ($fh,"<v2.3-unique-panel-genes.txt") || die "c";
while (<$fh>){
    s/[\r\n]//g;
    s/"//g;
    my $line = $_;
    my @lines = split ("\t",$line);
    my $gene = $lines[0];
    $newgenes{$gene} = 1;
    #print "$gene\n";
}
close $fh;

#foreach my $g (sort keys %newgenes){
#    print "$g\n";
#}

# drop mutations from CR1/ REL1 samples that are only in the new panel, and when the diag samples used the old panel
# print the final mutation analysis set
# Removes TPMT mutations as well
my $ofh2;
open ($ofh,">AML_data.muts.chroms.uniqFilt.dropDiscordantPanelGenes.txt") || die "c";
open ($ofh2,">dropped-variants.txt") || die "c";
open ($fh,"<AML_data.muts.chroms.uniqFilt.txt") || die "c";
my %checkids; # checks if a certain sample has been seen yet
while (<$fh>){
    s/[\r\n]//g;
    s/"//g;
    my $line = $_;
    my @lines = split("\t",$line);
    my $gene = $lines[1];
    
    if ($gene ne "TPMT"){
        if ($lines[4] ne "FUSION" && $lines[4] ne "CNA"){
            if ($lines[5] =~ /^dn/){ # for diagnosis samples, can print all variants
                print $ofh $line,"\n";
                $checkids{$lines[0]} = 1;
            }
            else { # nondiagnosis samples
                # if this pt had discordant NGS panels used
                if (defined $discordants{$lines[6]} ) { 
                    
                    # if the variant was from a new Panel-only gene, and this sample has not yet been printed:
                    if (defined $newgenes{$gene}) {
                        print $ofh2 $line,"\n"; # this goes to a separate file, just to keep track of removed variants
                        print $ofh "$lines[0]\t\t\tNA\t\t$lines[5]\t$lines[6]\n";
                        $checkids{$lines[0]} = 1;
                    }

                    if (!defined $newgenes{$gene}){ # if this gene variant is present in the old panel
                        $checkids{$lines[0]} = 1;
                        print $ofh $line,"\n";
                    }
                    
                }
                else { # samples were longitudinally sequenced with the same NGS panel
                    print $ofh $line,"\n";
                    $checkids{$lines[0]} = 1;
                }

            }
        }
        else { # if fusion or CNA, can print direct
            print $ofh $line,"\n";
            $checkids{$lines[0]} = 1;
        }
    }
    else {
        if (!defined $checkids{$lines[0]}){
            print $ofh "$lines[0]\t\t\t\t\t$lines[5]\t$lines[6]\n";
            $checkids{$lines[0]} = 1;
        }
    }

}

close $fh;
close $ofh;
close $ofh2;


