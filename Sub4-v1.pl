#!/usr/bin/perl
# File: Sub4-v1.pl
$|=1; #buffer off
use strict;
use warnings;

# PseudoCode/strategy
# Convert string to array, look for START (ATG) and STOP codons (TAA, TAG, TGA)
# NOTES: need to test and address exceptions if any, 

# Main to test Sub4
# Here we can do this x6 for each frame or call a new Sub
my $testFrame1 = "ATGTATTAAATGTGAGGCCCATGATCATAACATAACTGTGTATGTCTTAGAGGACCAAACCCCCCTCCTTCC"; #this comes as input from elsewher
my @tF1 = split ('', $testFrame1);

my @resultFrame1 = sub4(@tF1); # call Sub4
my $numberOfORFsF1 = (scalar(@resultFrame1)/2);

my $minORFLenght = 50; # min ORF lenght, Note: this is not yet in use, can be calculate while (Stop - Start position > $minORFLenght)
# Reformat results as needed

print "Number of ORFs: ".$numberOfORFsF1."\n";
print "Result array contents (Start/Stop position pairs): @resultFrame1\n";
#Regenerate the sequences based on Start/Stop, e.g while i < $numberOfORFsF1, then use start stop pair in a for to construct sequence.

exit; # end of main program

# Sub4
# Primary author: Roy Carambula
sub sub4 { #change to better name
    my @seq = @_;
    my @results = ();
    my $len = scalar(@seq);
    my $hits = 0; # tracks # of start stop position pairs
    my $i = 0;
    while ($i < $len) { # Find start ATG 
        if (($seq[$i] eq 'A') && ($seq[$i+1] eq 'T') && ($seq[$i+2] eq 'G')) {
            my $start = $i + 1;
            print "Start position is $start\n"; #will be removed
            my $j = $i + 3;
            while ($j < $len) { # Find stop  (TAA, TAG, TGA)
                if ((($seq[$j] eq 'T') && ($seq[$j+1] eq 'A') && ($seq[$j+2] eq 'A'))
                 || (($seq[$j] eq 'T') && ($seq[$j+1] eq 'A') && ($seq[$j+2] eq 'G'))
                 || (($seq[$j] eq 'T') && ($seq[$j+1] eq 'G') && ($seq[$j+2] eq 'A'))) {
                    my $stop = $j + 1;
                    $results[$hits] = $start;
                    $results[$hits+1] = $stop;
                    print "Stop position is $stop\n"; #will be removed
                    $hits = $hits + 2;
                    $j = $len;
                }
                $j = $j + 3;
                }
        }
        $i = $i + 3;
    }
    return(@results); # return array pairs of start, stop position values
} # End Sub4
