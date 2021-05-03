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

gs-report
version: $version
Usage: gstools gs-report [parameters]

Description: Read rnaseq-summary folder and generate gs-fisher results and summary reports.

Mandatory Parameters:
    --in|-i           rnaseq-summary folder
    --tx|-t           Transcriptome
                        Current support Human.B38.Ensembl88, Mouse.B38.Ensembl88
    --config|-c       Configuration file for gstools-db, by default
                        Human: /data/jyin/Databases/gstools-db/gstools-db-config_human.txt
                        Mouse: /data/jyin/Databases/gstools-db/gstools-db-config_mouse.txt

    --top|-t          # of top terms to show in figures [20]
	
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
	
my $inputfolder;
my $tx;
my $configfile;
my $outputfolder;
my $topnum=20;
my $dev=0;
	
my $verbose=1;


GetOptions(
	"input|i=s" => \$inputfolder,
	"tx|t=s" => \$tx,
	"config|c=s" => \$configfile,
	"top=s"=>\$topnum,	
	"out|o=s"=>\$outputfolder,
	"dev" => \$dev,	
);


###
#Tools needed
####

#omictools
my $mergefiles="/apps/omictools/mergefiles/mergefiles_caller.pl";
my $text2excel="perl /apps/omictools/text2excel/text2excel.pl";


#gstools
my $gstoolsfolder="/apps/gstools/";

#adding --dev switch for better development process
if($dev) {
	$gstoolsfolder="/home/centos/Pipeline/gstools/";
}

my $gs_fisher="$gstoolsfolder/gs-fisher_caller.pl";
my $gs_fisher_summary="$gstoolsfolder/gs-fisher_summary.pl";



#IO
if(!-e $outputfolder) {
	mkdir($outputfolder);
}

$outputfolder = abs_path($outputfolder);
my $outputfoldername = basename($outputfolder);

$inputfolder = abs_path($inputfolder);


my $scriptfolder="$outputfolder/scripts";
if(!-e $scriptfolder) {
	mkdir($scriptfolder);
}

my $scriptfile1="$scriptfolder/gs-report_run.sh";

my $logfile="$outputfolder/gs-report_run.log";


#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();

print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";

print LOG "\n";

print LOG "\ngstools gs-report $version running ...\n\n";
print STDERR "\ngstools gs-report $version running ...\n\n" if $verbose;



########
#Process
########

#find gstools-db config file
my %tx2configfile=(
	"Human.B38.Ensembl88"=>"/data/jyin/Databases/gstools-db/gstools-db-config_human.txt",
	"Mouse.B38.Ensembl88"=>"/data/jyin/Databases/gstools-db/gstools-db-config_mouse.txt",
);


my %tx2ref=(
	"Human.B38.Ensembl88"=> { 
		"star"=>"/data/jyin/Databases/Genomes/Human/hg38/Human.B38.Ensembl88_STAR",
		"rsem"=>"/data/jyin/Databases/Genomes/Human/hg38/Human.B38.Ensembl88_STAR/Human_RSEM",
		"chrsize"=>"/data/jyin/Databases/Genomes/Human/hg38/Human.B38.Ensembl88_STAR/chrNameLength.txt",
		"fasta"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.dna.primary_assembly_ucsc.fa",
		"gtf"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.88_ucsc.gtf",
		"homeranno"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.88_ucsc_homeranno.txt",
		"geneanno"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.88_ucsc_gene_annocombo.txt",
		"txanno"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.88_ucsc_tx_annocombo.txt"},
	"Mouse.B38.Ensembl88"=>{ 
		"star"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mouse.B38.Ensembl88_STAR",
		"rsem"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mouse.B38.Ensembl88_STAR/Mouse_RSEM",
		"chrsize"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mouse.B38.Ensembl88_STAR/chrNameLength.txt",
		"fasta"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.dna.primary_assembly_ucsc.fa",
		"gtf"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.88_ucsc.gtf",
		"homeranno"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.88_ucsc_homeranno.txt",
		"geneanno"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.88_ucsc_gene_annocombo.txt",
		"txanno"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.88_ucsc_tx_annocombo.txt"}
);


unless(defined $configfile && length($configfile) >0) {
	$configfile=$tx2configfile{$tx};
}

#check rnaseq-summary folder

if(!-e "$inputfolder/rnaseq-summary_GeneDESigs.txt") {
	print STDERR "ERROR:Can't read $inputfolder/rnaseq-summary_GeneDESigs.txt. $!\n";
	print LOG "ERROR:Can't read $inputfolder/rnaseq-summary_GeneDESigs.txt. $!\n";
	exit;
}

print STDERR "$inputfolder/rnaseq-summary_GeneDESigs.txt is analyzed by gs-report.\n\n";
print LOG "$inputfolder/rnaseq-summary_GeneDESigs.txt is analyzed by gs-report.\n\n";

#choose db to work on
my %dbfiles;

open(IN,$configfile) || die $!;
while(<IN>) {
	tr/\r\n//d;
	next if $_=~/^Database/;
	my @array=split/\t/;
	
	if($array[5] eq "Y") {
		$dbfiles{$array[3]}=$array[0];
		print STDERR $array[0]," is used for gs-report.\n";
		print LOG $array[0]," is used for gs-report.\n";
	}
}
close IN;


#print scripts

open(S1,">$scriptfile1") || die $!;
foreach my $dbfile (sort keys %dbfiles) {
	#both
	print S1 "perl $gs_fisher -i $inputfolder/rnaseq-summary_GeneDESigs.txt -t matrix -c both -o $outputfolder/$dbfiles{$dbfile} -d $dbfile -a ",$tx2ref{$tx}{"geneanno"},print_dev($dev,";");
	
	print S1 "perl $gs_fisher_summary -i $outputfolder/$dbfiles{$dbfile} -o $outputfolder/$dbfiles{$dbfile}_summary --top $topnum",print_dev($dev,";\n");

	#updown
	print S1 "perl $gs_fisher -i $inputfolder/rnaseq-summary_GeneDESigs.txt -t matrix -c updown -o $outputfolder/$dbfiles{$dbfile}_updown -d $dbfile -a ",$tx2ref{$tx}{"geneanno"},print_dev($dev,";");
	print S1 "perl $gs_fisher_summary -i $outputfolder/$dbfiles{$dbfile}_updown -o $outputfolder/$dbfiles{$dbfile}_updown_summary --top $topnum",print_dev($dev,";\n");

}
close S1;

########
#Functions
########

sub current_time {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $now = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $mday, $hour, $min, $sec);
	return $now;
}

sub print_dev {
	my ($dev_value,$toprint)=@_;
	
	if($dev_value) {
		return " --dev$toprint";
	}
	else {
		return "$toprint";
	}
}