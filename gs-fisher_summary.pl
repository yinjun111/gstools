#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);


#get an input folder 



########
#Prerequisites
########

#my $r="/apps/R-3.4.1/bin/R";
#my $rscript="/apps/R-3.4.1/bin/Rscript";

#Application version
my $mergefiles="/apps/sbptools/mergefiles/mergefiles_caller.pl";
my $text2excel="perl /apps/sbptools/text2excel/text2excel.pl";

#dev version
#my $mergefiles="perl /home/jyin/Projects/Pipeline/sbptools/mergefiles/mergefiles_caller.pl";

########
#Interface
########


my $version="0.4";

#v0.2, add --type to accept different input types
#v0.3, add text2excel to merge outputs
#v0.3a, add pvalue


my $usage="

gstools gs-fisher-summary
version: $version
Usage: gstools gs-fisher-summary [parameters]

Description: Merge multiple gs-fisher results


Mandatory Parameters:
    --in|-i           Input folder(s)					  
    --out|-o          Output folder
    --sigonly|-s      Only significant results [T]	
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

my $verbose=1;
my $runmode="none";

GetOptions(
	"in|i=s" => \$infolders,
	"out|o=s"=>\$outfolder,
	"sigonly|s=s"=>\$sigonly,
	"verbose"=>\$verbose,	
);


#mkpath for recursive make dir

print STDERR "\nsbptools gs-fisher $version running ...\n\n" if $verbose;

if(!-e $outfolder) {
	print STDERR "$outfolder not found. mkdir $outfolder\n";
	mkdir($outfolder); #doesn't do recursive
}

$outfolder=abs_path($outfolder);
my $outfoldername=basename($outfolder);

#the order of dirname and abs_path is important
my $outfile = abs_path($outfolder)."/".$outfoldername.".xlsx";


my $logfile=$outfile;
$logfile=~s/\.\w+$/_gs-fisher_summary.log/;


#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();

print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";

print LOG "\n";


print LOG "\nsbptools gs-fisher $version running ...\n\n";




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
#number, gene, OR, bhp, 
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
open(OUT4,">$outfile4") || die $!;

my $outfile5=$outfile;
$outfile5=~s/\.\w+$/_p.txt/;
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

system("$text2excel -i $outfile1,$outfile2,$outfile3,$outfile4,$outfile5 -n Number,Gene,OR,Q,P -o $exceloutfile --theme theme2");


close LOG;

########
#Function
########
sub current_time {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $now = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $mday, $hour, $min, $sec);
	return $now;
}
