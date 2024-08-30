#!/usr/bin/perl
use strict;
use warnings;
$|=1;
# convert timescape files to clevRvis files
# first get the hierarchical relationships
my $fh;
my @ptinds;
open ($fh,"<dotlist.txt") || die "c";
while (<$fh>){
    s/[\r\n]//g;
    my $line = $_;
    my @info = split ("\Q.",$line);
    my $pt = $info[1];
    push(@ptinds,$pt);
}
close $fh;

foreach my $ptind (@ptinds){

    #my $ptind = 98;
    my $fh2;
    my $ofh;
    my %labels; # captures the gene variants associated with a clone
    my %parents; # captures the parental clone

    open ($fh2,"<./output-res/pt.$ptind.tree1.dot") || die "c";
    while (<$fh2>){
        s/[\r\n]//g;
        my $line = $_;

        # if this line describes the clone labels:
        if ($line =~ /\Q[label/){
            my @dat = split (" ",$line);
            my $clone = $dat[0];
            my $label = $dat[1];
            $label =~ s/\Q[label="//g;
            $label =~ s/"\Q];//g;

            #print "$clone\t$label\n";
            $labels{$clone} = $label;
        }

        # if this line describes the clone hierarchy:
        if ($line =~ /->/){
            my @dat = split (" -> ",$line);
            my $root = $dat[0];
            my $rootn = $labels{$root};
            my $destination = $dat[1];
            $destination =~ s/;//g;
            my $destn = $labels{$destination};
            $parents{$destn} = $rootn;
        }
    }
    close $fh2;
  
    ####################
    # sum up the CCFs for the daughter clones
    open ($fh2,"<./output-res/pt.$ptind.soln1.csv") || die "c";
    my $crcheck = 0;
    while (<$fh2>){
        s/[\r\n]//g;

        if ($_ =~ ","){
            my @lines = split (",",$_);
            if ($lines[0] eq "CR1"){
                $crcheck = 1;
            }
        }
    }
    close $fh2;


    open ($fh2,"<./output-res/pt.$ptind.soln1.csv") || die "c";
    open ($ofh,">./output-res/pt.$ptind.soln1.clevRvis.txt") || die "a";
   <$fh2>;

    my %hh; # stores merged clone Names
    my %bb; # stores column indices. CloneName => #index
    my $head = <$fh2>;
    $head =~ s/[\r\n]//g;
    print "================\n",$ptind,"\n";
 
    my @headers2 = split (",",$head);
    print $ofh "time";
    for (my $i =1; $i < scalar @headers2; $i++){
        $bb{$headers2[$i]} = $i;
        print $ofh "\t$headers2[$i]";
    }
    print $ofh "\n";
    my $linemax = 0;
    if ($crcheck == 1){
        $linemax = 4;
    }
    elsif ($crcheck == 0){
        $linemax = 3;
    }

    my $linecount = 0;
    while ($linecount < $linemax && defined(my $line = <$fh2>)){
        $line =~ s/,/\t/g;
        print $ofh $line;
        $linecount++;
    }
    close $fh2;
    close $ofh;

    ######
    # Now print the parental clone indices in order
    open ($ofh,">./output-res/pt.$ptind.soln1.clevRvis.indices.txt") || die "a";
    for (my $i = 1; $i < scalar @headers2; $i++){
        if (defined $parents{$headers2[$i]}){
            print $ofh "$headers2[$i]\t$bb{$parents{$headers2[$i]}}\n";
        }
        else {
            print $ofh "$headers2[$i]\t0\n";
        }
    }
    close $ofh;

}
