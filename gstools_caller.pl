#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);

my $version="0.2";

#v0.1
#v0.2, added z-score calculation and new dot plot

my $usage="

gstools
version: $version
Usage: gstools [tool] [parameters]

Parameters:

    gs-report    Automatically generate Gene Set Analyses results for rnaseq-summary
    gs-fisher    Perform Fisher's Exact Test for gene lists
    gs-fisher-summary   Summarize gs-fisher results
    list-dbs     List current databases for gs-fisher

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


my $gs_report="perl $gstoolsfolder/gs-report_caller.pl";
my $gs_fisher="perl $gstoolsfolder/gs-fisher_caller.pl";
my $gs_fisher_summary="perl $gstoolsfolder/gs-fisher_summary.pl";
my $list_dbs="sh $gstoolsfolder/list-dbs.sh";


my %commands2program=(
    "gs-report"=>$gs_report,
    "gs-fisher"=>$gs_fisher,
	"gs-fisher-summary"=>$gs_fisher_summary,
	"list-dbs"=>$list_dbs,
);



if(defined $commands2program{$command}) {
	system($commands2program{$command}." ".join(" ",@params));
}
else {
	print STDERR "ERORR $command not found in gstools.\n\n";
	system("gstools");
}
