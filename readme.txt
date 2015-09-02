This repository contains a series of scripts used for pre-processing high-throughput DNA sequences
prior to analysis. Some of these require other software packages or modules installed and in your path:
BioPerl			http://www.bioperl.org/wiki/Main_Page
cross_match.manyreads	http://www.phrap.org/phredphrapconsed.html

The following usage information for each script can be reproduced by running each script without arguments.

Any questions or bugs? Please contact Eli Meyer: eli.meyer@science.oregonstate.edu.
-------------------------
AdaptorFilterFastq.pl
-------------------------

Searches for reads matching a set of adaptor sequences supplied by the user.
Reads matching adaptors are discarded.
Usage:	 AdaptorFilterFastq.pl sequences adaptors min_score output
Arguments:
	 sequences	 file of short reads to be filtered, fastq format
	 adaptors	 file of adaptor sequences to screen for, fasta format
	 min_score	 score threshold; alignments scoring this high (no. bp) are removed
	 output		 a name for the output file (fastq format)

-------------------------
AdaptorTrimFastq.pl
-------------------------

Filters a set of short reads in FASTQ format, trimming away regions
matching the specified adaptor sequences
Usage:	 AdaptorTrimFastq.pl sequences adaptors min_bp min_score output
Arguments:
	 sequences	 file of short reads to be filtered, fastq format
	 adaptors	 file of adaptor sequences to screen for, fasta format
	 min_score	 score threshold; alignments scoring this high (no. bp) are removed
	 min_length	 minimum length; reads shorter than this after trimming are discarded
	 output		 a name for the output file (fastq format)

-------------------------
FastqToFasta.pl
-------------------------
Converts a fastq file from Illumina into fasta sequence and quality score files.
Usage:	 FastqToFasta.pl fastq out_fasta out_qual
Arguments:
	 fastq	 name of fastq input file 
	 out_fasta	 name of fasta output file 
	 out_qual	 name of qual output file 

-------------------------
HRFilterFastq.pl
-------------------------

Filters a FASTQ file to remove sequences containing homopolymer
repeats (HR) longer than the specified threshold
Usage:	 HRFilterFastq.pl sequences crit_length output
Arguments:
	 sequences	 file of short reads to be filtered, fastq format
	 crit_length	 reads containing HRs longer than this will be excluded
	 output		 a name for the output file (fastq format)

-------------------------
QualFilterFastq.pl
-------------------------

Removes reads containing too many low quality basecalls from a set of short sequences 
Output:	 high-quality reads in FASTQ format
Usage:	 QualFilterFastq.pl input.fastq low_score min_LQ output.fastq
Arguments:
	 input.fastq	 raw input reads in FASTQ format
	 low score	 quality scores below this are considered low quality (LQ)
	 min_LQ		 reads with more than this many LQ bases are excluded
	 output.fastq	 name for ourput file of HQ reads in FASTQ format

-------------------------
TruncateFastq.pl
-------------------------

Truncates a set of short reads in FASTQ format to keep the region specified
Output:	 FASTQ formatted sequences
Usage:	 script sequences start_position end_position output
Arguments:
	 sequences	 file of short reads to be filtered, fastq format
	 start_position	 beginning of the region to keep, nucleotide position
	 end_position	 end of the region to keep, nucleotide position
	 output		 a name for the output file (fastq format)

