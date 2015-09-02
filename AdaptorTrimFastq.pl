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
unless ($#ARGV == 4)
        {print "\nFilters a set of short reads in FASTQ format, trimming away regions\n";
	print "matching the specified adaptor sequences\n";
        print "Usage:\t $scriptname sequences adaptors min_bp min_score output\n";
        print "Arguments:\n";
        print "\t sequences\t file of short reads to be filtered, fastq format\n";
        print "\t adaptors\t file of adaptor sequences to screen for, fasta format\n";
        print "\t min_score\t score threshold; alignments scoring this high (no. bp) are removed\n";
        print "\t min_length\t minimum length; reads shorter than this after trimming are discarded\n";
        print "\t output\t\t a name for the output file (fastq format)\n";
        print "\n"; exit;
        }

my $seqfile = $ARGV[0];		# raw reads, fastq format
my $adfile = $ARGV[1];		# adaptors, fasta format
my $maxn = 0;			# max number of Xs allowed in screened file
my $minmatch = 6	;	# min word match for alignment, controls sensitivity
my $minscore = 12	;	# min Smith Waterman score to report in cross_match.log
my $minaln = $ARGV[2];		# min score (no. matching bases) of alignment
my $minlen = $ARGV[3];		# min length for trimmed reads
my $outfile = $ARGV[4];		# name for output file, fastq format

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
	if ($count==2) {$rdlen = length($_);}
	if ($_ =~ /^$ss/) {$_ =~ s/\@//; print TMP ">", $_, "\n"; $switch = 1; next;}
	if ($_ =~ /^\+/) {$switch = 0; next;}
	if ($switch == 1) {print TMP $_, "\n";}
	}
close(IN);
close(TMP);

# adaptor searching
system("cross_match.manyreads tmp.fasta $adfile -minmatch $minmatch -minscore $minscore -screen >cross_match.log 2>cross_match.errors.log");

# identification of good regions
# first build a hash of match positions
open(CMT,"cross_match.log");
my %cmh = ();
while (<CMT>)
	{
	chomp;
	unless ($_ =~ /^\s+\d+\s+\d/)	{next;}
	@bits = split(/ +/, $_);
	$scorei = ($bits[7] - $bits[6]) * (100-($bits[2] + $bits[3] + $bits[4]))/100;
	if ($scorei < $minaln) {next;}
	$cmh{$bits[5]}{$bits[6]} = $bits[7];
#	print join(", ", @bits), "\t", $scorei, "\n";
	}

foreach $rd (sort(keys(%cmh)))
	{
	%rdh = %{$cmh{$rd}};
	@ssa = sort{$a<=>$b}(keys(%rdh));
	$mno = keys(%rdh);
	if ($mno==1)	
		{
#		$keeph{$rd}{0} = $start; $keeph{$rd}{1} = $rdh{$start}; next;
# reads with a single match to adaptor	
		$dista = $ssa[0];
		$distb = $rdlen - $rdh{$ssa[0]};
		if ($distb < $dista) {$pickstart = 1; $pickend = $ssa[0] - 1;}
		else	{$pickstart = $rdh{$ssa[0]} + 1; $pickend = $rdlen;}
#		print $rd, "\t", $ssa[0], "\t", $rdh{$ssa[0]}, "\t", $pickstart, "\t", $pickend, "\n";
		if ($pickend-$pickstart+1 >= $minlen) 
			{
			$trunch{$rd}{0} = $pickstart;
			$trunch{$rd}{1} = $pickend;
			}
		else	{$tooshort{$rd}++;}
		next;	
		}
# reads with > 1 match handled differently
	%gaph = (); $maxno = @ssa;
	for ($sno=0; $sno<$maxno; $sno++)
		{
		$starti = $ssa[$sno];
		$endi = $rdh{$starti};
		if ($sno<$maxno-1)
			{
			$nexti = $ssa[$sno+1];
			$gapi = $nexti - $rdh{$starti};
			$gaph{$starti}{"gap"} = $gapi;
			$gaph{$starti}{"start"} = $rdh{$starti};
			$gaph{$starti}{"end"} = $nexti;
#			print $rd, "\t", $starti, "\t", $rdh{$starti}, "\t", $nexti, "\t", $gapi, "\n";
			}
		else	{
#			print $rd, "\t", $starti, "\t", $rdh{$starti}, "\n";
			}

		}
	@gapa = sort{$gaph{$b}{"gap"}<=>$gaph{$a}{"gap"}}(keys(%gaph));
	$pickstart = $gaph{$gapa[0]}{"start"} + 1;
	$pickend = $gaph{$gapa[0]}{"end"} - 1;
#	print $rd, "\t", "multiple hits", "\t", $pickstart, "\t", $pickend, "\n";
	if ($pickend-$pickstart+1 >= $minlen) 
		{
		$trunch{$rd}{0} = $pickstart;
		$trunch{$rd}{1} = $pickend;
		}
	else	{$tooshort{$rd}++;}
	}

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
		if(exists($tooshort{$_})) {$switch = 0; $bcount++;}
		else {$switch = 1; $goodcount++;}
		$sid = $_;
		$incount++;
		}
	if ($switch == 1) 
		{
		if ($part == 2 || $part == 4)
			{
			if(exists($trunch{$sid}))
				{
				$prtstr = substr($_, $trunch{$sid}{0}-1, $trunch{$sid}{1}-$trunch{$sid}{0}+1);
				}
			else	{$prtstr = $_;}
			}
		else	{$prtstr = $_;}
		if ($part == 1) {print OUT "@";}
		print OUT $prtstr, "\n";
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

