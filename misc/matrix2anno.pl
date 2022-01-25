#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);


########
#Prerequisites
########

########
#Interface
########


my $version="0.1";


my $usage="

matrix2anno
version: $version
Usage: gstools matrix2anno [parameters]

Description: convert matrix information to gene annotation


Mandatory Parameters:
    --in|-i           a matrix with each row as a gene and each 
                         column as a comparison. 1 and -1 as up/down-regulated genes
                         E.g. rnaseq-summary_folder/rnaseq-summary_GeneDESigs.txt
	
    --out|-o          annotation file used by gs-fisher
                          first column as gene, second col as annotation, sep by \";\"
	
";


unless (@ARGV) {
	print STDERR $usage;
	exit;
}

my $params=join(" ",@ARGV);
#then call different scripts



########
#Parameters
########
	
my $inputfile;
my $outputfile;
	
my $verbose=1;


GetOptions(
	"input|i=s" => \$inputfile,
	"output|o=s" => \$outputfile,
		
);

########
#Program starts
########

my $linenum=0;
my @cates;
my %gene2cate;

open(IN,$inputfile) || die $!;

open(OUT,">$outputfile") || die $!;

print OUT "Feature\tAnno\n";

while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	
	if($linenum==0) {
		@cates=@array[1..$#array];
	}
	else {
		for(my $num=1;$num<@array;$num++) {
			if($array[$num]>=1) {
				$gene2cate{$array[0]}{$cates[$num-1]."-Up"}++;
				$gene2cate{$array[0]}{$cates[$num-1]}++;
			}
			
			if($array[$num]<=-1) {
				$gene2cate{$array[0]}{$cates[$num-1]."-Down"}++;
				$gene2cate{$array[0]}{$cates[$num-1]}++;
			}
		}

		
	}
	
	if(defined $gene2cate{$array[0]}) {
		print OUT $array[0],"\t",join(";",sort keys %{$gene2cate{$array[0]}}),"\n";
	}
	
	$linenum++;
}
close IN;
close OUT;


