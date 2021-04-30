#!/usr/bin/perl -w
use strict;

#cutoff
my $pcutoff=5e-8;


#get all the current genes
my %genes;
open(IN,"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.88_ucsc_gene_annocombo.txt") || die $!;
while(<IN>) {
	tr/\r\n//d;
	next if $_=~/^Gene/;
	
	my @array=split/\t/;
	
	$genes{$array[0]}++;
}
close IN;

#GWAS cata
my %gene2pheno_filtered;
my %gene2pheno_all;

open(IN,"/data/jyin/Databases/gstools-db/GWASCatalog/gwas_catalog_v1.0.2-associations_e100_r2021-04-20.tsv") || die $!;
while(<IN>) {
	tr/\r\n//d;
	next if $_=~/^DATE/;
	my @array=split/\t/;
	
	#upstream, downstream, snpgene
	foreach my $mapped (@array[15..17]) {
		foreach my $gene (split(", ",$mapped)) {
			if($array[27] < $pcutoff) {
				$gene2pheno_filtered{$gene}{$array[7]}++;
			}
			$gene2pheno_all{$gene}{$array[7]}++;
		}
	}

}
close IN;



open(OUT1,">/data/jyin/Databases/gstools-db/GWASCatalog/gwas_catalog_v1.0.2_5e-8_human.txt") || die $!;
print OUT1 "Gene\tAnnotation\n";

open(OUT2,">/data/jyin/Databases/gstools-db/GWASCatalog/gwas_catalog_v1.0.2_all_human.txt") || die $!;
print OUT2 "Gene\tAnnotation\n";

foreach my $gene (sort keys %genes) {
	print OUT1 $gene,"\t";
	print OUT2 $gene,"\t";
	
	if(defined $gene2pheno_filtered{$gene}) {
		print OUT1 join(";",sort keys %{$gene2pheno_filtered{$gene}}),"\n";
	}
	else {
		print OUT1 " \n";
	}
	
	if(defined $gene2pheno_all{$gene}) {
		print OUT2 join(";",sort keys %{$gene2pheno_all{$gene}}),"\n";
	}
	else {
		print OUT2 " \n";
	}
}

close OUT1;
close OUT2;
