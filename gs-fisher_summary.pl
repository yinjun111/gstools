#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);


########
#Interface
########


my $version="0.1";

#v0.1

my $usage="

gstools gs-fisher-summary
version: $version
Usage: gstools gs-fisher-summary [parameters]

Description: Merge multiple gs-fisher results


Mandatory Parameters:
    --in|-i           Input folder(s)					  
    --out|-o          Output folder
    --sigonly|-s      Only significant results [T]	
    --top|-t          # of top terms to show in figures [20]	
    --verbose|-v      Verbose
	
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
	
my $infolders;
my $outfolder;

my $sigonly="T";
my $topnum=20;

my $verbose=1;
my $runmode="none";
my $dev=0; #developmental version


GetOptions(
	"in|i=s" => \$infolders,
	"out|o=s"=>\$outfolder,
	"sigonly|s=s"=>\$sigonly,
	"top|t=s"=>\$topnum,	
	"verbose"=>\$verbose,
	"dev" => \$dev,		
);

########
#Prerequisites
########

my $r="/apps/R-4.0.2/bin/R";
my $rscript="/apps/R-4.0.2/bin/Rscript";

#Application version


my $omictoolsfolder="/apps/omictools/";

#adding --dev switch for better development process
if($dev) {
	$omictoolsfolder="/home/centos/Pipeline/omictools/";
#}
#else {
	#the tools called will be within the same folder of the script
	#$omictoolsfolder=get_parent_folder(abs_path(dirname($0)));
}

my $mergefiles="$omictoolsfolder/mergefiles/mergefiles_caller.pl";
my $text2excel="perl $omictoolsfolder/text2excel/text2excel.pl";


#gstools
my $gstoolsfolder="/apps/gstools/";

#adding --dev switch for better development process
if($dev) {
	$gstoolsfolder="/home/centos/Pipeline/gstools/";
#}
#else {
	#the tools called will be within the same folder of the script
	#$omictoolsfolder=get_parent_folder(abs_path(dirname($0)));
}

my $gs_heatmap="$gstoolsfolder/gs_heatmap.R";
my $gs_dotplot="$gstoolsfolder/gs_dotplot.R";


########
#Program running
########

#mkpath for recursive make dir

print STDERR "\nomictools gs-fisher $version running ...\n\n" if $verbose;

if(!-e $outfolder) {
	print STDERR "$outfolder not found. mkdir $outfolder\n";
	mkdir($outfolder); #doesn't do recursive
}

$outfolder=abs_path($outfolder);
my $outfoldername=basename($outfolder);

#the order of dirname and abs_path is important
my $outfile = abs_path($outfolder)."/".$outfoldername.".xlsx";


my $logfile="$outfolder/gs-fisher_summary.log";


#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();

print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";

print LOG "\n";

print LOG "\nomictools gs-fisher_summary $version running ...\n\n";


########
#Process
########

#or,BHP, gene, number

my %file2gs;
my %gss;
my %files;

foreach my $infolder (split(",",$infolders)) {
	foreach my $file (glob("$infolder/*")) {
		if(length($file)>2) {
			#open each file and test
			open(IN,$file) || die $!;
			my $filename=basename($file);
			
			my $linenum=0;
			while(<IN>) {
				tr/\r\n//d;
				if($linenum==0) {
					unless($_=~/^Gene set/) {
						last;
					}
					$files{$filename}++;
					print STDERR $file," detected from gstools.\n";
				}
				else {
					my @array=split/\t/;
					
					$file2gs{$filename}{$array[0]}=[@array];
					
					if($sigonly eq "T") {
						#keeps only significant results
						if($array[7] eq "TRUE") {
							$gss{$array[0]}++;
						}
					}
					else {
						$gss{$array[0]}++;
					}
				}
				
				$linenum++;
			}
			close IN;
		}
	}
}

#2,3,5,6
#number, gene, OddsRatio, bhp, 
my @outfiles;

my $outfile1=$outfile;
$outfile1=~s/\.\w+$/_num.txt/;
open(OUT1,">$outfile1") || die $!;

my $outfile2=$outfile;
$outfile2=~s/\.\w+$/_gene.txt/;
open(OUT2,">$outfile2") || die $!;

my $outfile3=$outfile;
$outfile3=~s/\.\w+$/_or.txt/;
open(OUT3,">$outfile3") || die $!;

my $outfile4=$outfile;
$outfile4=~s/\.\w+$/_q.txt/;
push @outfiles,$outfile4;
open(OUT4,">$outfile4") || die $!;

my $outfile5=$outfile;
$outfile5=~s/\.\w+$/_p.txt/;
push @outfiles,$outfile5;

open(OUT5,">$outfile5") || die $!;

print OUT1 "Gene set\t",join("\t",sort keys %files),"\n";
print OUT2 "Gene set\t",join("\t",sort keys %files),"\n";
print OUT3 "Gene set\t",join("\t",sort keys %files),"\n";
print OUT4 "Gene set\t",join("\t",sort keys %files),"\n";
print OUT5 "Gene set\t",join("\t",sort keys %files),"\n";

foreach my $gs (sort keys %gss) {
	print OUT1 $gs,"\t";
	print OUT2 $gs,"\t";
	print OUT3 $gs,"\t";
	print OUT4 $gs,"\t";
	print OUT5 $gs,"\t";
	
	my (@marks1,@marks2,@marks3,@marks4,@marks5);
	
	foreach my $file (sort keys %file2gs) {
		if(defined $file2gs{$file}{$gs}) {
			push @marks1,$file2gs{$file}{$gs}[2];
			push @marks2,$file2gs{$file}{$gs}[3];
			push @marks3,$file2gs{$file}{$gs}[5];
			push @marks4,$file2gs{$file}{$gs}[6];
			push @marks5,$file2gs{$file}{$gs}[4];
		}
		else {
			push @marks1," ";
			push @marks2," ";
			push @marks3," ";
			push @marks4," ";		
			push @marks5," ";		
		}
	}
	print OUT1 join("\t",@marks1),"\n";
	print OUT2 join("\t",@marks2),"\n";
	print OUT3 join("\t",@marks3),"\n";
	print OUT4 join("\t",@marks4),"\n";	
	print OUT5 join("\t",@marks5),"\n";		
}
close OUT1;
close OUT2;
close OUT3;
close OUT4;
close OUT5;


#Merge files into excel file
my $exceloutfile=$outfile;
$exceloutfile=~s/\.\w+$/\.xlsx/;

print STDERR "Converting outputs to excel file:$exceloutfile.\n\n";
print LOG "Converting outputs to excel file:$exceloutfile.\n\n";

system("$text2excel -i $outfile1,$outfile2,$outfile3,$outfile4,$outfile5 -n Number,Gene,OddsRatio,Q,P -o $exceloutfile --theme theme2");



#draw heatmap and dotplots

#Generate heatmap figures for p/q

foreach my $outfile_sel (@outfiles) {
	my $outfilefig=$outfile_sel;
	$outfilefig=~s/\.\w+$/_heatmap.png/;	

	print STDERR "Generating heatmap for $outfile_sel.\n";
	print LOG "Generating heatmap for $outfile_sel.\n";

	system("$rscript $gs_heatmap -i $outfile_sel -o $outfilefig --top $topnum");
	print LOG "$rscript $gs_heatmap -i $outfile_sel -o $outfilefig --top $topnum\n";
	
}

#Generate dotplot for or+p and or+q
#OddsRatio+p

my $outfiledpfig1=$outfile;
$outfiledpfig1=~s/\.\w+$/_orbhp_dotplot.png/;
	
system("$rscript $gs_dotplot --size $outfile3 --sizename OddsRatio --color $outfile4 --colorname \"'-Log10BHP'\" -o $outfiledpfig1 --top $topnum");
print LOG "$rscript $gs_dotplot --size $outfile3 --sizename OddsRatio --color $outfile4 --colorname \"'-Log10BHP'\" -o $outfiledpfig1 --top $topnum\n";


my $outfiledpfig2=$outfile;
$outfiledpfig2=~s/\.\w+$/_orrawp_dotplot.png/;

system("$rscript $gs_dotplot --size $outfile3 --sizename OddsRatio --color $outfile5 --colorname \"'-Log10RawP'\" -o $outfiledpfig2 --top $topnum");
print LOG "$rscript $gs_dotplot --size $outfile3 --sizename OddsRatio --color $outfile5 --colorname \"'-Log10RawP'\" -o $outfiledpfig2 --top $topnum\n";

close LOG;

########
#Function
########
sub current_time {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $now = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $mday, $hour, $min, $sec);
	return $now;
}
