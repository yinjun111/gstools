#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);

my $version="0.1";

#v0.1


my $usage="

gstools
version: $version
Usage: gstools [tool] [parameters]


Parameters:

    ########
    #RNA-Seq
    ########
    rnaseq-process    RNA-seq QC, Align, and RSEM for FASTQ files
    rnaseq-merge      Merge rnaseq-process results for downstream analyses
    rnaseq-de         Perform DE analysis using DESeq2
    rnaseq-summary    Summarize RNA-Seq DE results

    rnaseq-var        RNA-seq variant calling pipeline
    rnaseq-motif      RNA-seq TFBS motif finding pipeline
    rnaseq-motif-summary  RNA-seq TFBS motif finding results summary	
    motif-finder      Transcription factor binding motif prediction


";


unless (@ARGV) {
	print STDERR $usage;
	exit;
}


#functions

my ($command,@params)=@ARGV;

####
#check whether to use --dev version
####
my $params_list=join(" ",@params);

my $dev=0; #developmental version

if($params_list=~/--dev/) {
	$dev=1;
}


my $gstoolsfolder="/apps/gstools/";

#adding --dev switch for better development process
if($dev) {
	$gstoolsfolder="/home/centos/Pipeline/gstools/";
}
else {
	#the tools called will be within the same folder of the script
	$gstoolsfolder=dirname(abs_path($0));
}



#####
#then call different scripts
#####


my $bs_fastq="$gstoolsfolder/bs-fastq/bs-fastq_caller.pl";
my $geo_download="$gstoolsfolder/geo-download/geo-download_caller.pl";


my %commands2program=(
    "bs-fastq"=>$bs_fastq,
	"geo-download"=>$geo_download,
	
    "rnaseq-process"=>$rnaseq_process,
    "rnaseq-merge"=>$rnaseq_merge,
    "rnaseq-de"=>$rnaseq_de,
	"rnaseq-summary"=>$rnaseq_summary,
);



if(defined $commands2program{$command}) {
	system($commands2program{$command}." ".join(" ",@params));
}
else {
	print STDERR "ERORR $command not found in gstools.\n\n";
	system("gstools");
}
