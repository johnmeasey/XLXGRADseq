# Make a phylogeny from the tab delimited genotypes

First we need to convert the tab file to a nexus file using this script(21_tab_to_interleave.pl):
```perl
#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;



#  This program reads in a tab delimited genotype file generated
#  by the perl program '17_adds_outgroup_to_lots_of_tab_files.pl'
#  and generates an interleaved nexus file that includes degenerate bases and gaps

# it is hardcoded to expect only one base from the reference seq


# run it like this
# 21_tab_to_interleave_nexus.pl input.tab output_interleave.nxs



my $inputfile = $ARGV[0];
my $outputfile = $ARGV[1];

unless (open DATAINPUT, $inputfile) {
	print 'Can not find the input file.\n';
	exit;
}


my @temp;
my @temp1;
my @names;
my %datahash;
my $y;
my $x;
my $watisitnow;
my $count=0;
my $interleave=0;

# Read in datainput file

# Ideally, I'd like to print out the interleaved sections as the vcf file is read
# and keep track of the number of bases, and print this out to screen at the end
# then this number could be added to the nxs file after the script runs


while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split('\t',$line);
	if($temp[0] eq '#CHROM'){
		@names=@temp;
		for ($y = 2; $y <= $#names; $y++ ) {
			$names[$y] =~ s/-//gi; # get rid of dashes in names
			$datahash{$names[$y]}=''; # initialize the hash
		}	
		# print preamble to output file
		unless (open(OUTFILE, ">$outputfile"))  {
			print "I can\'t write to $outputfile\n";
			exit;
		}
		print "Creating output file: $outputfile\n";
		print OUTFILE "#NEXUS\n\n";
		print OUTFILE "BEGIN DATA\;\nDIMENSIONS NTAX=",$#names-1," NCHAR= XXXX\;\n";
		print OUTFILE "FORMAT DATATYPE=DNA  MISSING=? INTERLEAVE GAP=- \;\n";
		print OUTFILE "MATRIX\n";
		print OUTFILE "\n";
	}
	else{	
		# only print ones that are not microsats or indels in the outgroups
		if($interleave>79){  # print this section of the data
			for ($y = 2; $y <= $#names; $y++ ) {
				print OUTFILE $names[$y],"\t\t",$datahash{$names[$y]},"\n";
			}
			print OUTFILE "\n";
			$interleave=0;
			# clear the hash
			for ($y = 2 ; $y <= $#names; $y++ ) {
				$datahash{$names[$y]}='';
			}	
		}
		if(length($temp[2]) == 1){ # all the outgroup seqs are single bp
			$count=$count+1; # this is the count of all positions
			$interleave=$interleave+1; # this is the count of the interleave length
			
			
			# now add data to the hash
			for ($y = 2 ; $y <= 2; $y++ ) {	# first the three outgroups which are haploid
				if($temp[$y] ne '*'){
					$datahash{$names[$y]} = $datahash{$names[$y]}.uc($temp[$y]);
				}
				else{
					$datahash{$names[$y]} = $datahash{$names[$y]}.'-';
				}	
			}
			for ($y = 3 ; $y <= $#temp; $y++ ) { # now the ingroups, which are diploid, usually (except chrX and chrY)
				# for these, we need to use IUPAC codes
				if(($temp[$y] eq 'G/G')||($temp[$y] eq 'C/C')||($temp[$y] eq 'T/T')||($temp[$y] eq 'A/A')||($temp[$y] eq 'G/')||($temp[$y] eq 'C/')||($temp[$y] eq 'T/')||($temp[$y] eq 'A/')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.substr($temp[$y],0,1);
				}
				elsif(($temp[$y] eq './.')||($temp[$y] eq './')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.'-';
				}
				elsif(($temp[$y] eq 'C/T')||($temp[$y] eq 'T/C')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.'Y';
				}
				elsif(($temp[$y] eq 'A/G')||($temp[$y] eq 'G/A')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.'R';
				}
				elsif(($temp[$y] eq 'A/C')||($temp[$y] eq 'C/A')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.'M';
				}
				elsif(($temp[$y] eq 'A/T')||($temp[$y] eq 'T/A')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.'W';
				}
				elsif(($temp[$y] eq 'C/G')||($temp[$y] eq 'G/C')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.'S';
				}
				elsif(($temp[$y] eq 'G/T')||($temp[$y] eq 'T/G')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.'K';
				}
				else{ # this is a microsat, so substitute missing data
					$datahash{$names[$y]} = $datahash{$names[$y]}.'-';
				}
			}	
		}
	}
}		

# print the last line
for ($y = 2; $y <= $#names; $y++ ) {
	print OUTFILE $names[$y],"\t\t",$datahash{$names[$y]},"\n";
}
#print OUTFILE "\n";



print OUTFILE "\;\nEND\;";
print OUTFILE "\n";
print OUTFILE "\n";
#print OUTFILE "BEGIN Mrbayes\;\n";
#print OUTFILE "Prset statefreqpr=dirichlet(1,1,1,1)\;\n";
#print OUTFILE "Lset  nst=6  rates=invgamma\;\n";
#print OUTFILE "mcmc ngen=2000000 savebrlens=yes\;\n";
#print OUTFILE "sumt burnin=10000\;\n";
#print OUTFILE "quit\;\n";

close OUTFILE;
print "The number of sites is $count\n";

# now update the number of bases

my $status;
$status = system("perl -p -i -e 's/XXXX/$count/g' $outputfile");
```
# Analysis with iqtree

`~/2015_SulaRADtag/good_merged_samples/iqtree-1.5.0a-Linux/bin/iqtree -s XLXVXGXM_merged_sorted.bam.vcf.gz.phy -m TEST -nt 1 -pre XLXVXGXM_merged_sorted.bam.vcf.gz.phy_`

# Ultrafast bootstrap
`iqtree -s XLXVXGXM_merged_sorted.bam.vcf.gz.phy -m TEST -bb 1000`





