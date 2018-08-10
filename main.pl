use strict;
use warnings;

#Sub 1: Read and extract FASTA from file - Author: Marcus Agostinelli
##Used Beginning Perl for Bioinformatics by James D.Tisdall for reference
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

# Sub 3 - Author: Erick Galinkin
# Take in a string and return an array with the appropriate 3-gram.
sub convert3Gram {
  my $input = $_[0];
  my $frame = $_[1];
  my @seq = ();
  for (my $i=($frame-1); $i<length($input)-3; $i+=3){
    my $threeGram = substr($input, $i, 3);
    push @seq, $threeGram;
  }
  return @seq
}

# Printresults sub - Author: Erick Galinkin
# Cleanly print results from the frames.
sub printresults {
  my @frame = @{$_[0]};
  my $frameID = $_[1];
  my $numberofORFs = (scalar(@frame)/2);
  print "$numberofORFs ORFs found in $frameID\n";
  for(my $i=0; $i < ($numberofORFs*2); $i++){
    my $start = $frame[$i];
    $i++;
    my $stop = $frame[$i];
    print "Start: $start     Stop: $stop\n";
  }
  print "\n";
}

# Sub4
# Primary author: Roy Carambula
## Comment from Erick: Once we have the sequence and the reverse complement of the sequence,
## we pass the sequence to the findORF subroutine. This will give us the start and stop positions of the ORF.

# PseudoCode/strategy
# Convert string to array, look for START (ATG) and STOP codons (TAA, TAG, TGA)
# NOTES: need to test and address exceptions if any,
sub findORF { #change to better name
    ## Comment from Erick: Integrating my sub into this so we can avoid complex array references.
    ## Takes input as a string per output from Courtney and Marcus's subs and converts it to an array.
    ## Adding ability to put in Frame number.
    my $input = $_[0];
    my @seq = convert3Gram($input, $_[1]);
    my @results = ();
    my $len = scalar(@seq);
    my $hits = 0; # tracks # of start stop position pairs
    my $i = 0;
    while ($i < $len) { # Find start ATG
        if ($seq[$i] eq 'ATG') {
            my $start = $i + 1;
            my $j = $i + 1;
            while ($j < $len) { # Find stop  (TAA, TAG, TGA)
                if (($seq[$j] eq 'TAA') || ($seq[$j] eq 'TAG') || ($seq[$j] eq 'TGA')) {
                    my $stop = ($j * 3) + 1;
                    $results[$hits] = $start;
                    $results[$hits+1] = $stop;
                    $hits = $hits + 2;
                    $j = $len;
                }
                $j = $j + 1;
                }
        }
        $i = $i + 1;
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

#6 reading frames.
my @frame1 = findORF($sequenceData, 1);
my @frame2 = findORF($sequenceData, 2);
my @frame3 = findORF($sequenceData, 3);
my @frameneg1 = findORF($revcompSeq, 1);
my @frameneg2 = findORF($revcompSeq, 2);
my @frameneg3 = findORF($revcompSeq, 3);

printresults(\@frame1, "Frame +1");
printresults(\@frame2, "Frame +2");
printresults(\@frame3, "Frame +3");
printresults(\@frameneg1, "Frame -1");
printresults(\@frameneg2, "Frame -2");
printresults(\@frameneg3, "Frame -3");

exit;
