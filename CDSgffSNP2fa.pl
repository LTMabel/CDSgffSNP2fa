#!/usr/bin/perl -w

####################################################################################################
#
#		Sarah B. Kingan
#		University of Rochester
#		16 January 2014
#
#		Title: CDSgffSNP2fa.pl
#
#
#		This program generates fasta sequences using SNP and GFF inputs
#	
#			
#
####################################################################################################

use strict;
use Tabix;
use Bio::Seq;
use Bio::DB::Fasta;
use Data::Dumper;

my $usage = "CDSgffSNP2fa.pl <file.gff> <file.snp.bgz> <samples.txt> <ref.fa>\n";


# get file names from user input:
my $gff_file = shift@ARGV or die $usage;
my $snp_file = shift@ARGV or die $usage;
my $samples_file = shift@ARGV or die $usage;
my $ref_fasta_file = shift@ARGV or die $usage;

# read through gff file and save data in hash
my %gff_hash;
open (GFF, $gff_file);
while (my $line = <GFF>) {
	my @line_array = split("\t", $line);
	my $chrom = $line_array[0];
	my $start = $line_array[3];
	my $end = $line_array[4];
	my $strand = $line_array[6];
	my $FBtr;
	if ($line_array[8] =~ /Parent=((\S)+);/) {
		$FBtr = $1;
	}
	elsif ($line_array[8] =~ /(TCONS_[0-9]{8})/) {
		$FBtr = $1;
	}
	my @FBtr_array = split(",", $FBtr);
	for (my $i = 0; $i< scalar(@FBtr_array); $i++) {
		push(@{$gff_hash{$chrom}{$FBtr_array[$i]}{'start'}},$start);
		push(@{$gff_hash{$chrom}{$FBtr_array[$i]}{'end'}},$end);
		$gff_hash{$chrom}{$FBtr_array[$i]}{'strand'} = $strand;
	}
}

foreach my $chrom (sort keys %gff_hash) {
	foreach my $FBtr (sort keys %{$gff_hash{$chrom}}) {
		my %final_seq_hash;
		open (FBTR, ">$FBtr.fa") or die "Couldn't open file $FBtr.fa, $!";
# loop through each interval and get sequence hash
		for (my $i = 0; $i<scalar(@{$gff_hash{$chrom}{$FBtr}{'start'}}); $i++) {
			my $start = ${$gff_hash{$chrom}{$FBtr}{'start'}}[$i];
			my $end = ${$gff_hash{$chrom}{$FBtr}{'end'}}[$i];
			my %seq_hash = get_seq_hash($chrom, $start, $end);
# loop through each sample and concatenate sequence onto final seq hash
			foreach my $sample (sort keys %seq_hash) {
				my $short_sample;
				if ($sample =~ /(D[a-z]+_(\S)+)\|/) {
					$short_sample = $1;
				}
				elsif ($sample =~ /(\S+)\|/) {
					$short_sample = $1;
				}
				$final_seq_hash{$short_sample} .= $seq_hash{$sample};
			}
		}
		foreach my $sample (sort keys %final_seq_hash) {
			my $seq_obj;
			my $RC_obj;
			my $seq = $final_seq_hash{$sample};
# reverse complement as needed
			if ($gff_hash{$chrom}{$FBtr}{'strand'} eq "-") {
				$seq_obj = Bio::Seq->new(-seq => $seq, -alphabet => 'dna' );
				$RC_obj = $seq_obj->revcom;
				$seq = $RC_obj->seq();
			}
			print FBTR ">", $sample, "\n";
			print FBTR $seq, "\n";
		}
	}
}



###################### SUB ROUTINES ###########################

sub get_seq_hash {
	my ($chrom, $start, $end) = @_;
	my $output = `PBsnp2fa.pl $snp_file $ref_fasta_file $chrom:$start-$end $samples_file`;
	my @output_array = split("\n", $output);
	my %seq_hash;
	my $name;
	foreach my $line (@output_array) {
		if ($line =~ /^>((\S)+)/) {
			$name = $1;
		}
		else {
			$seq_hash{$name} .= $line;
		}
	}
	return %seq_hash;		
}


exit;
