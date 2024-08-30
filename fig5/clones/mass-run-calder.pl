#!/usr/bin/perl
use strict;
use warnings;
my $fh;
open ($fh,"<file-list.txt") || die "c";
while (<$fh>){
    s/[\r\n]//g;
    my $file = $_;
    my $pt = $file;
    $pt =~ s/.txt//g;

    chdir("./calder/");
    
    print "==================\n$file\n";
    system("java -jar calder.jar -i ../input-matrices/$file -o ../output-res -v glpk");
    system("mv ../output-res/pt_soln1.csv ../output-res/$pt.soln1.csv");
    system("mv ../output-res/pt_tree1.dot ../output-res/$pt.tree1.dot");

    #print "==================\n$file\n";
    #system("python soln_to_timescape.py ../output-res/$pt.soln1.csv ../output-res/$pt.tree1.dot ../output-res/$pt");
}
close $fh;
