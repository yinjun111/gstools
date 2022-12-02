#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);

my $version="0.4";

#v0.1
#v0.2, added z-score calculation and new dot plot
#v0.3, changed 2x2 table
#v0.31, add format conversion tools
#v0.4, add pathview/sbgnview

my $usage="

gstools
version: $version
Usage: gstools [tool] [parameters]

Parameters:

    gs-report    Automatically generate Gene Set Analyses results for rnaseq-summary

    gs-fisher    Perform Fisher's Exact Test for gene lists
    gs-fisher-summary   Summarize gs-fisher results

    gsea-gen     Run GSEA analysis
    gsea-gen-summary   Summarize GSEA analysis results

    ipa-gen      Generate files for IPA analysis (not implemented)
    ipa-summary  Summary IPA analysis results

    gs-keggpathview/gs-kegg  Generate pathway figure for KEGG using Pathview
    gs-sbgnview  Generate pathway figure for Reactome/Wikipathways using Sbgnview

    #format conversion tools
    list2matrix  Convert gene lists by columns to matrix format
    matrix2ms    Convert matrix format to Metascape format

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

my $Rscript=find_program("/apps/R-4.0.2/bin/Rscript");


my $gs_report="perl $gstoolsfolder/gs-report/gs-report_caller.pl";
my $gs_fisher="perl $gstoolsfolder/gs-fisher/gs-fisher_caller.pl";
my $gs_fisher_summary="perl $gstoolsfolder/gs-fisher/gs-fisher_summary.pl";
my $list_dbs="sh $gstoolsfolder/misc/list-dbs.sh";

my $gsea_gen="perl $gstoolsfolder/gsea-gen/gsea-gen_caller.pl";
my $gsea_gen_summary="perl $gstoolsfolder/gsea-gen/gsea-gen-summary.pl";

my $ipa_summary="perl $gstoolsfolder/ipa-gen/ipa_summary.pl";

my $gs_kegg="$Rscript $gstoolsfolder/gs-keggpathview/gs_keggpathview.R";

my $list2matrix="perl $gstoolsfolder/misc/list2matrix.pl";
my $matrix2ms="perl $gstoolsfolder/misc/matrix2ms.pl";


my %commands2program=(
    "gs-report"=>$gs_report,
    "gs-fisher"=>$gs_fisher,
	"gs-fisher-summary"=>$gs_fisher_summary,
	"gsea-gen"=>$gsea_gen,	
	"gsea-gen-summary"=>$gsea_gen_summary,
	"ipa-summary"=>$ipa_summary,
	"gs-kegg"=>$gs_kegg,
	"gs-keggpathview"=>$gs_kegg,
	"gsea-gen"=>$gsea_gen,		
	"list2matrix"=>$list2matrix,
	"matrix2ms"=>$matrix2ms,	
	"list-dbs"=>$list_dbs,
);



if(defined $commands2program{$command}) {
	if(@params) {
		system($commands2program{$command}." ".join(" ",@params));
	}
	else {
		if($command=~/gs-kegg/) {
			system($commands2program{$command}." --help");
		}
		else {
			system($commands2program{$command}." ".join(" ",@params));
		}
	}
}
else {
	print STDERR "ERORR $command not found in gstools.\n\n";
	system("gstools");
}


#####
#Functions
#####


sub find_program {
	my $fullprogram=shift @_;
	
	my $program;
	if($fullprogram=~/([^\/]+)$/) {
		$program=$1;
	}
	
	if(-e $fullprogram) {
		return $fullprogram;
	}
	else {
		my $sysout=`$program`;
		if($sysout) {
			my $location=`which $program`;
			return $location;
		}
		else {
			print STDERR "ERROR:$fullprogram or $program not found in your system.\n\n";
			exit;
		}
	}
}
