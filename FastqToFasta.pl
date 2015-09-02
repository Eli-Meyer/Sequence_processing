#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check for dependencies and use modules
$mod1="Bio::SeqIO";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;} 
use Bio::SeqIO;
$mod2="Bio::Seq::Quality";
unless(eval("require $mod2")) {print "$mod2 not found. Exiting\n"; exit;} 
use Bio::Seq::Quality;

$scriptname=$0; $scriptname =~ s/.+\///g;

# -- program description and required arguments
unless ($#ARGV == 2)
        {print "Converts a fastq file from Illumina into fasta sequence and quality score files.\n";
        print "Usage:\t $scriptname fastq out_fasta out_qual\n";
        print "Arguments:\n";
        print "\t fastq\t name of fastq input file \n";
        print "\t out_fasta\t name of fasta output file \n";
        print "\t out_qual\t name of qual output file \n";
        print "\n"; exit;
        }

my $wrap = 1000;		# edit to be longer than maximum read length
my $fastqfile = $ARGV[0];
my $inseqs = new Bio::SeqIO(-file=>$fastqfile, -format=>"fastq");
my $osfile = $ARGV[1];
my $oqfile = $ARGV[2];
my $outseqs = new Bio::SeqIO(-file=>">$osfile", -format=>"fasta", -width=>$wrap);
my $outquals = new Bio::SeqIO(-file=>">$oqfile", -format=>"qual", -width=>$wrap);

my %sh; my $scount = 0;
while ($seq = $inseqs->next_seq) 
	{$sh{$seq->display_id} = $seq->seq;
	$outseqs->write_seq($seq);
	$outquals->write_seq($seq);
	$scount++;
	}

print "Converted ", $scount, " reads from FASTQ to FASTA and QUAL.\n";

