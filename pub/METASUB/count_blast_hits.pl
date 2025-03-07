#!/usr/bin/perl
use warnings;
#==============================================================================
# count_blast_hits.pl / created by Ryo Honda (END-AMR-Asia), 2023-05-27
#==============================================================================
# This perl script outputs read counts of each gene from a blast result file by:
#	$ perl count_blast_hits.pl blast_results.txt
#
# Output format (-outfmt) of blast results should be specified as:
# -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
#------------------------------------------------------------------------------

open(IN1,$ARGV[0]);
my %hash_counts;
my %hash_slen;
my $rpk;
while(<IN1>){
	chomp;
	@col=split("\t",$_);
	if(exists $hash_counts{$col[1]}){
		$hash_counts{$col[1]}++;
	}else{
		$hash_counts{$col[1]}=1;
		$hash_slen{$col[1]}=$col[13];
	}
}
print "sseqid\treads\tslen\tRPK\n";
for my $keys_gene (keys %hash_counts){ 
	$rpk=$hash_counts{$keys_gene}/$hash_slen{$keys_gene}*1000;
	print $keys_gene,"\t",$hash_counts{$keys_gene},"\t",$hash_slen{$keys_gene},"\t",$rpk,"\n";
}
close(IN1);
