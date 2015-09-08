#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- program description and required arguments
$scriptname=$0; $scriptname =~ s/.+\///g;
if ($#ARGV != 2 || $ARGV[0] eq "-h")
        {print "\nFilters a FASTQ file to remove sequences containing homopolymer\n";
	print "repeats (HR) longer than the specified threshold\n";
        print "Usage:\t $scriptname sequences crit_length output\n";
        print "Arguments:\n";
        print "\t sequences\t file of short reads to be filtered, fastq format\n";
        print "\t crit_length\t reads containing HRs longer than this will be excluded\n";
        print "\t output\t\t a name for the output file (fastq format)\n";
        print "\n"; exit;
        }

my $seqfile = $ARGV[0];		# raw reads, fastq format
my $critlen = $ARGV[1];		# critical length
my $outfile = $ARGV[2];		# name for output file, fastq format
my @vals = qw{A C G T N};
my @stra = ();
foreach $v (@vals)
	{
	$ssi = $v x $critlen;
	push (@stra, $ssi);
	}

# loop through fastq file and print out only those passing filter
open (IN, $seqfile);
open (OUT, ">$outfile");
my $switch = 0;
while(<IN>)
	{
	chomp;
	$count++;
	if ($count==1) {$ss = substr($_, 0, 8);}
	if ($_ =~ /^$ss/) 
		{
		$thisname = $_;
		$nextup=1;
		$switch=0;
		}
	else
		{
		if ($nextup==1)
			{
			$switch=1;
			foreach $s (@stra) 
				{
				if ($_ =~ /$s/) {$switch=0;}
				}			
			$nextup=0;
			if ($switch==0) {$fail++;}
			elsif ($switch==1) 
				{
				print OUT $thisname, "\n";
				$pass++;
				}
			}
		if ($switch>0) {print OUT $_, "\n";}
		}
	}
close(IN);
print $fail, " reads failed\n";
print $pass, " reads passed\n";
system("date");
