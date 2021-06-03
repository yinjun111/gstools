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


my $version="0.1a";

#0.1a, change names to first letter uppercase


my $usage="

gs-anno_msigdb
version: $version
Usage: gstools gs-anno-msigdb [parameters]

Description: Annotate Ensembl IDs by msigdb


Mandatory Parameters:
    --in|-i           Ensembl Gene IDs
    --chip|-c         msigdb chip file (/data/jyin/Databases/GSEA/msigdb_v7.4_chip_files_to_download_locally/Human_Gene_Symbol_with_Remapping_MSigDB.v7.4.chip)
    --gmt|-g          msigdb gmt file /data/jyin/Databases/GSEA/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_GMTs/h.all.v7.4.symbols.gmt
    --out|-o          Annotation used by gstools
	
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
my $chipfile; #gene to gs annotation
my $gmtfile; #gene id to symbol
my $outfile;
	
my $verbose=1;
my $runmode="none";

GetOptions(
	"input|i=s" => \$infile,
	"gmt|g=s" => \$gmtfile,
	"chip|c=s" => \$chipfile,
	"out|o=s"=>\$outfile,
	"verbose"=>\$verbose,	
);


print STDERR "\ngstools gs-anno-msigdb $version running ...\n\n" if $verbose;


#the order of dirname and abs_path is important
$outfile = abs_path($outfile);


my $logfile=$outfile;
$logfile=~s/.txt/_gs-anno-msigdb_run.log/;


#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();

print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";

print LOG "\n";


print LOG "\ngstools gs-anno-msigdb $version running ...\n\n";




########
#Process
########

my %genes; #genes in current annotation

open(IN,$infile) || die "ERROR:Error reading $infile.$!\n";
while(<IN>) {
	tr/\r\n//d;
	next if $_=~/^Gene/;
	my @array=split/\t/;
	
	$genes{$array[0]}++;
}
close IN;

print LOG scalar(keys %genes)," genes defined in $infile.\n";
print STDERR scalar(keys %genes)," genes defined in $infile.\n";

#gene 2 symbol form chip
my %gene2symbol;

open(IN,$chipfile) || die "ERROR:Error reading $chipfile.$!\n";
while(<IN>) {
	tr/\r\n//d;
	next if $_=~/^Gene/;
	my @array=split/\t/;
	
	$gene2symbol{$array[0]}=$array[1];
}
close IN;

print LOG scalar(keys %gene2symbol)," genes defined in $chipfile.\n";
print STDERR scalar(keys %gene2symbol)," genes defined in $chipfile.\n";

#symbol 2 anno from GMT

my %symbol2anno;

open(IN,$gmtfile) || die "ERROR:Error reading $gmtfile.$!\n";
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
		
	my $gs;
	
	#upper case first letter of first word #doesn't work very well
	#if($array[0]=~/^([^_]+_)(\w)(.+)/) {
	#	$gs=$1.uc($2).lc($3);
	#}
	#else {
	#	print STDERR "WARNING:Wrong format $gs\n";
		$gs=$array[0];
	#}
	
	for(my $num=2;$num<@array;$num++) {
		$symbol2anno{$array[$num]}{$gs}++;
	}
}
close IN;

print LOG scalar(keys %symbol2anno)," genes defined in $gmtfile.\n";
print STDERR scalar(keys %symbol2anno)," genes defined in $gmtfile.\n";


#output
open(OUT,">$outfile") || die "ERROR:Error reading $outfile.$!\n";
print OUT "Gene\tAnnotation\n";
my $finalnum=0;
foreach my $gene (sort keys %genes) {
	#only output gene defined in both orignial anno and chip anno
	if(defined $gene2symbol{$gene}) {
		if(defined $symbol2anno{$gene2symbol{$gene}}) {
			print OUT $gene,"\t",join(";",sort keys %{$symbol2anno{$gene2symbol{$gene}}}),"\n";
		}
		else {
			print OUT $gene,"\t \n";
		}
		
		$finalnum++;
	}
}

close OUT;
		
print LOG "$finalnum genes defined in both $infile and $chipfile.\n";
print STDERR "$finalnum genes defined in both $infile and $chipfile.\n";


#########
#Functions
#########


sub current_time {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $now = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $mday, $hour, $min, $sec);
	return $now;
}
