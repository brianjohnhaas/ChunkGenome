#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;

use lib ("$FindBin::Bin/../PerlLib");
use Fasta_reader;



my $fasta_file = $ARGV[0] || die "\n\n\tusage: $0 genome.fasta > N_regions.tsv\n\n";

my $fasta_reader = new Fasta_reader($fasta_file);
while (my $seq_obj = $fasta_reader->next()) {

	my $acc = $seq_obj->get_accession();
	my $sequence = $seq_obj->get_sequence();


	while ($sequence =~ /(N+)/gi) {
		my $start = $-[0];
		my $stop = $+[0];
		
		$start++;

		print "$acc\t$start\t$stop\n";
	}
	
}

exit(0);

