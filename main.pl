use strict;
use warnings;

#Sub 1: Read and extract FASTA from file - Author: Marcus Agostinelli
## Comment from Erick: read file and return an array with the contents of the file.
## Then pass the output from read_fasta into extract_fasta.
## extract_fasta returns a string with the sequence in FASTA format.
sub read_fasta {
  my($file)=@_;
  my @fasta_data=();
  unless(open(READ_FASTA,$file)){
    print STDERR "Cannot open file $file \n\n";
    }
    @fasta_data=<READ_FASTA>;
    close READ_FASTA;
    return @fasta_data;
}

sub extract_fasta {
  my(@fasta)= @_;
  my $seq='';
  foreach my $line (@fasta){
    if($line =~ /^\s*$/){
    next;
    } elsif($line =~ /^\s*#/){
      next;
      } elsif($line =~ /^>/){
      next;
      } else{
      $seq .=$line;
      }
    }
    $seq =~ s/\s//g;
    return $seq;
}

#Sub 2: Create a reverse complementary strand of DNA - Author: Courtney Kelley
## Comment from Erick: Take output from extract_fasta and pass it to revcomp_strand.
## This will provide us the reverse complement of the sequence in addition to the sequence.
sub revcomp_strand{

    #Create variable for DNA string input
    ## Comment from Erick: Modifying to read string.
    my $forward = $_[0];

    #Reverse the string
    my $reverse = reverse($forward);

    #Copy the reverse strand into the revcomp variable
    my $revcomp = $reverse;

    #Create the complementary strand
    $revcomp =~ tr/ATGC/TACG/;

    #Return the new string
    return $revcomp;

}
# Sub4
# Primary author: Roy Carambula
## Comment from Erick: Once we have the sequence and the reverse complement of the sequence,
## we pass the sequence to the findORF subroutine. This will give us the start and stop positions of the ORF.

# PseudoCode/strategy
# Convert string to array, look for START (ATG) and STOP codons (TAA, TAG, TGA)
# NOTES: need to test and address exceptions if any,
sub findORF { #change to better name
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

## Main program execution - Author: Erick Galinkin
my $infile;
if (@ARGV) {
  $infile = $ARGV[0];
  print "Reading data from $infile...\n";
} else {
  print "Usage: orfs.pl <input_file>\n";
  exit;
}

my @data = read_fasta($infile);
my $sequenceData = extract_fasta(@data);
my $revcompSeq = revcomp_strand($sequenceData);
#Program works up to this point.
exit;

# Main to test Sub4
# Here we can do this x6 for each frame or call a new Sub
my $testFrame1 = "ATGTATTAAATGTGAGGCCCATGATCATAACATAACTGTGTATGTCTTAGAGGACCAAACCCCCCTCCTTCC"; #this comes as input from elsewher
my @tF1 = split ('', $testFrame1);

my @resultFrame1 = findORF(@tF1); # call Sub4
my $numberOfORFsF1 = (scalar(@resultFrame1)/2);

my $minORFLength = 50; # min ORF lenght, Note: this is not yet in use, can be calculate while (Stop - Start position > $minORFLenght)
# Reformat results as needed

print "Number of ORFs: ".$numberOfORFsF1."\n";
print "Result array contents (Start/Stop position pairs): @resultFrame1\n";
#Regenerate the sequences based on Start/Stop, e.g while i < $numberOfORFsF1, then use start stop pair in a for to construct sequence.

exit; # end of main program
