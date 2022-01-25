#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);
#use Data::Dumper qw(Dumper);

########
#Interface
########


my $version="0.1";



my $usage="

gstools list2matrix
version: $version
Usage: gstools list2matrix [parameters]

Description: Convert gene list in by columns into matrix


Mandatory Parameters:
    --in|-i           Input folder(s)					  
    --out|-o          Output folder
	
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
	
my $infile;
my $outfile;


GetOptions(
	"in|i=s" => \$infile,
	"out|o=s"=>\$outfile,
);


my @cates;
my %cate2gene;
my %genes;
my $linenum=0;

open(IN,$infile) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;

	if($linenum==0) {
		@cates=@array;
	}
	else {
		#save gene to cate
		for(my $num =0;$num<@cates;$num++) {
			if(defined $array[$num] && length($array[$num])>0 && $array[$num] ne " ") {
				$cate2gene{$cates[$num]}{$array[$num]}++;
				$genes{$array[$num]}++;
			}
		}
	}
		
	$linenum++;
}
close IN;

open(OUT,">$outfile") || die $!;
print OUT "Gene\t",join("\t",@cates),"\n";
foreach my $gene (sort keys %genes) {
	print OUT $gene,"\t";
	
	#doesn't differentiate 1/-1
	my @marks;
	
	for(my $num =0;$num<@cates;$num++) {
		if(defined $cate2gene{$cates[$num]}{$gene}) {
			push @marks,1;
		}
		else {
			push @marks,0;
		}
	}
	
	print OUT join("\t",@marks),"\n";
}
close OUT;
	

