#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions
$scriptname=$0; $scriptname =~ s/.+\///g;

# -- check for dependencies
$mod1="File::Which";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;} 
use File::Which;
$dep1 = "cross_match.manyreads";
unless (defined(which($dep1))) {print $dep1, " not found. Exiting.\n"; exit;}

# -- program description and required arguments
unless ($#ARGV == 3)
        {print "\nSearches for reads matching a set of adaptor sequences supplied by the user.\n";
	print "Reads matching adaptors are discarded.\n";
        print "Usage:\t $scriptname sequences adaptors min_score output\n";
        print "Arguments:\n";
        print "\t sequences\t file of short reads to be filtered, fastq format\n";
        print "\t adaptors\t file of adaptor sequences to screen for, fasta format\n";
        print "\t min_score\t score threshold; alignments scoring this high (no. bp) are removed\n";
        print "\t output\t\t a name for the output file (fastq format)\n";
        print "\n"; exit;
        }

my $seqfile = $ARGV[0];		# raw reads, fastq format
my $adfile = $ARGV[1];		# adaptors, fasta format
my $maxn = 0;			# max number of Xs allowed in screened file
my $minmatch = 6;		# min word match for sw alignment
my $minscore = 12;		# min score (matching bases) of alignment to report in cross match log
my $minaln = $ARGV[2];		# min score (no. matching bases) of alignment counted as a match
my $outfile = $ARGV[3];		# name for output file, fastq format

system("date");

# convert fastq to fasta
open (IN, $seqfile);
open (TMP, ">tmp.fasta");
my $switch = 0;
my $count = 0;
while(<IN>)
	{
	chomp;
	$count++;
	if ($count==1) {$ss = substr($_, 0, 4);}
	if ($_ =~ /^$ss/) {$_ =~ s/\@//; print TMP ">", $_, "\n"; $switch = 1; next;}
	if ($_ =~ /^\+/) {$switch = 0; next;}
	if ($switch == 1) {print TMP $_, "\n";}
	}
close(IN);
close(TMP);

# adaptor searching
system("cross_match.manyreads tmp.fasta $adfile -minmatch $minmatch -minscore $minscore -screen >cross_match.log 2>cross_match_errors.log");

# parsing adaptor matches
open(CMT,"cross_match.log");
my %cmh = ();
while (<CMT>)
	{
	chomp;
	unless ($_ =~ /^\s+\d+\s+\d/)	{next;}
	@bits = split(/ +/, $_);
	$scorei = ($bits[7] - $bits[6]) * (100-($bits[2] + $bits[3] + $bits[4]))/100;
	if ($scorei < $minaln) {next;}
	else	{$badh{$bits[5]}++;}
	}
close(CMT);

# finally loop through fastq file and write out trimmed sequences
open (IN, $seqfile);
open (OUT, ">$outfile");
my $switch = 0; $part = 0;
while(<IN>)
	{
	chomp;
	$part++;
	if ($part == 1)
		{
		$_ =~ s/\@//; $_ =~ s/\s.+//;
		if(exists($badh{$_})) {$switch = 0; $bcount++;}
		else {$switch = 1; $goodcount++;}
		$sid = $_;
		$incount++;
		}
	if ($switch == 1) 
		{
		if ($part == 1) {print OUT "@";}
		print OUT $_, "\n";
		}
	if ($part == 4)	{$part = 0;}
	}
close(IN);
close(OUT);

system("rm tmp.fasta*; rm cross_match.log");
print "Output from ", $scriptname, "\n";
print $incount, " reads input\n";
print $bcount, " reads failed\n";
print $goodcount, " reads passed\n";
system("date");

