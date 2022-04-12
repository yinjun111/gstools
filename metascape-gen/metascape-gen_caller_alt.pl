#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);
#use Data::Dumper qw(Dumper);

########
#Interface
########


my $version="0.2";



my $usage="

gstools metascape-gen
version: $version
Usage: gstools metascape-gen [parameters]

Description: Convert significance matrix to gene list format needed by Metascape, e.g. first column as experiment name, 2nd col as genes separated by ,


Mandatory Parameters:
    --in|-i           Input file					  
    --out|-o          Output file

    --cate|-c         Category of changes, both or updown [updown]
                        use updown to generate lists for up and down genes separately
                        use both to generate lists for up and down genes together

    --background|-b   Whether to use all genes from the matrix as background [T]
						
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
my $cateopt="updown";
my $background="T";

GetOptions(
	"in|i=s" => \$infile,
	"out|o=s"=>\$outfile,
	"cate|c=s"=>\$cateopt,
	"background|b=s"=>\$background,
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
		@cates=@array[1..$#array];
	}
	else {
		#save gene to cate
		for(my $num =0;$num<@cates;$num++) {
			
			my $cate=$cates[$num];
			
			if(defined $array[$num+1] && length($array[$num+1])>0 && $array[$num+1] ne " " && $array[$num+1] ne "NA") {
				
				if($array[$num+1]>0) {
					#any number >0
					if($cateopt eq "both") {
						$cate2gene{$cates[$num]}{$array[0]}++;
					}
					else {
						$cate2gene{$cates[$num]."-Up"}{$array[0]}++;
					}
				}
				elsif ($array[$num+1]<0) {
					#any number <0
					if($cateopt eq "both") {
						$cate2gene{$cates[$num]}{$array[0]}++;
					}
					else {
						$cate2gene{$cates[$num]."-Down"}{$array[0]}++;
					}
				}
			}
		}
		
		$genes{$array[0]}++;
	}
		
	$linenum++;
}
close IN;

open(OUT,">$outfile") || die $!;

foreach my $cate (sort keys %cate2gene) {
	print OUT $cate,"\t";
	print OUT join(",",sort keys %{$cate2gene{$cate}}),"\n";
}

#print background
if($background eq "T") {
	print OUT "_BACKGROUND\t",join(",",sort keys %genes),"\n";
}

close OUT;
	

