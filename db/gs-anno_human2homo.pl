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

gs-anno_human2homo
version: $version
Usage: gstools gs-anno_human2homo [parameters]

Description: Annotate Ensembl IDs by homolog info


Mandatory Parameters:
    --in|-i           Original Annotation
    --homo|-h         Homolog file
    --out|-o          Annotation by homolog
	
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
my $homofile; 
my $outfile;
	
my $verbose=1;


GetOptions(
	"input|i=s" => \$infile,
	"homo|h=s" => \$homofile,
	"out|o=s"=>\$outfile,
);


print STDERR "\ngstools gs-anno_human2homo $version running ...\n\n" if $verbose;


#the order of dirname and abs_path is important
$outfile = abs_path($outfile);


my $logfile=$outfile;
$logfile=~s/.txt/_gs-anno_human2homo_run.log/;


#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();

print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";

print LOG "\n";


print LOG "\ngstools gs-anno_human2homo $version running ...\n\n";




########
#Process
########

#human to anno
my %hgene2anno;
open(IN,$infile) || die "ERROR:Can't read $infile.\n";
while(<IN>) {
	tr/\r\n//d;
	next if $_=~/^Gene/;
	my @array=split/\t/;
	$hgene2anno{$array[0]}=$array[1];
}
close IN;


#other species to human
my %ogene2anno;
my %ogenes;

open(IN,$homofile) || die "ERROR:Can't read $homofile.\n";
while(<IN>) {
	tr/\r\n//d;
	next if $_=~/^Gene/;
	my @array=split/\t/;
	
	#if(defined $array[5] && $array[5] eq "ortholog_one2one") {
	#not limiting to one2one
	if(defined $array[2] ) {
		$ogenes{$array[0]}++;
		if(defined $hgene2anno{$array[2]}) {
			foreach my $term (split(";",$hgene2anno{$array[2]})) {
				next if $term eq " ";
				$ogene2anno{$array[0]}{$term}++;
			}
		}
	}
}
close IN;

#results
#limit to genes with homolog to human
open(OUT,">$outfile") || die "ERROR:Can't write $outfile.\n";
print OUT "Gene\tAnnotation\n";

foreach my $gene (sort keys %ogenes) {
	print OUT $gene,"\t";
	
	if(defined $ogene2anno{$gene}) {
		print OUT join(";",sort keys %{$ogene2anno{$gene}}),"\n";
	}
	else {
		print OUT " \n";
	}
}

close OUT;
	

########
#Function
########
sub current_time {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $now = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $mday, $hour, $min, $sec);
	return $now;
}
