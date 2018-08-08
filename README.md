# bifs617
Group project
#Reading Fasta!

use strict;
use warnings;

my @data=();
my $dna='';

#enter in the file name within the quotes 
@data=read_fasta("");

$dna=extract_fasta(@data);

print ($dna);

# sub for reading the FASTA file
sub read_fasta{
my($file)=@_;

my @fasta_data=();
unless(open(READ_FASTA,$file)){
print STDERR "Cannot open file\$file\"n\n";
}
@fasta_data=<READ_FASTA>;
close READ_FASTA;
return@fasta_data;

}

# sub for extracting DNA
sub extract_fasta{
my(@fasta)= @_;


my $seq='';

foreach my $line (@fasta){
#blank line
if($line =~ /^\s*$/){
next;
#comment line
} elsif($line =~ /^\s*#/){
next;
#header line
} elsif($line =~ /^>/){
next;

} else{
$seq .=$line;
}
}
$seq =~ s/\s//g;
return $seq;

}
