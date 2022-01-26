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

#v0.2, change default dotplot to z+p/q


my $usage="

gstools gs-fisher-summary
version: $version
Usage: gstools gs-fisher-summary [parameters]

Description: Merge multiple gs-fisher results


Mandatory Parameters:
    --in|-i           Input folder(s)					  
    --out|-o          Output folder
    --sigonly|-s      Only significant results [T]	

    --comparisons     (Optional) Name of comparisons to be tested in gs analyses
                         By default, all the comparisons in the gs-fisher folder

    --suffix|-x       Remove suffix from file names [_gs-fisher.txt]
    --top|-t          # of top terms to show in figures [20]
    --sortby|-y       Sort by which comparison [avg]
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
my $sortby="avg";
my $comparisons;
my $suffix="_gs-fisher.txt";

my $verbose=1;
my $runmode="none";
my $dev=0; #developmental version


GetOptions(
	"in|i=s" => \$infolders,
	"out|o=s"=>\$outfolder,
	"sigonly|s=s"=>\$sigonly,
	"comparisons=s" => \$comparisons,	
	"sortby=s" => \$sortby,
	"suffix|x=s"=>\$suffix,	
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

my $mergefiles="perl $omictoolsfolder/mergefiles/mergefiles_caller.pl";
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

my $gs_fisher="perl $gstoolsfolder/gs-fisher/gs-fisher_caller.pl";
my $gs_fisher_summary="perl $gstoolsfolder/gs-fisher/gs-fisher_summary.pl";


########
#Program running
########

#mkpath for recursive make dir

print STDERR "\ngstools gs-fisher $version running ...\n\n" if $verbose;

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

print LOG "\ngstools gs-fisher_summary $version running ...\n\n";


########
#Process
########

#pre-selected comparisons
my @selcomparisons;
my %selcomparisons_hash;
if(defined $comparisons && length($comparisons)>0) { 
	#keep original order
	#need to convert non-char in $comp 
	foreach my $comp (split(",",$comparisons)) {
		$comp=~s/\W/_/g;
		push @selcomparisons,$comp;
	}
	
	%selcomparisons_hash=map {$_,1} @selcomparisons;
}

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
			
			my $samplename=$filename;
			
			if(defined $suffix && length($suffix)>0) {
				$samplename=~s/$suffix//;
			}
			
			my $linenum=0;
			while(<IN>) {
				tr/\r\n//d;
				if($linenum==0) {
					unless($_=~/^Gene set/) {
						last;
					}
					$files{$samplename}++;
					print STDERR $file," detected from gstools.\n";
					print LOG $file," detected from gstools.\n";
					
					#check --comparisons
					if(@selcomparisons) {
						unless(defined $selcomparisons_hash{$samplename}) {
							print STDERR "$samplename is not defined in --comparisons $comparisons. Skip...\n";
							print LOG "$samplename is not defined in --comparisons $comparisons. Skip...\n";
							next;
						}
					}
				}
				else {
					my @array=split/\t/;
					
					$file2gs{$samplename}{$array[0]}=[@array];
					
					if($sigonly eq "T") {
						#keeps only significant results
						if($array[8] eq "TRUE") {
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

#check comparisons vs detected samplename

my @usedcomparisons;

if(@selcomparisons) {
	foreach my $comp (@selcomparisons) {
		unless(defined $file2gs{$comp}) {
			print STDERR "ERROR:$comp was not found in --in $infolders.\n";
			print LOG "ERROR:$comp was not found in --in $infolders.\n";
			exit;
		}
	}
	@usedcomparisons=@selcomparisons;
}
else {
	@usedcomparisons=sort keys %file2gs;
}


#error message

if(@usedcomparisons ==0) {
	print STDERR "\nERROR:No gs-fisher results are detected in $infolders\n\n";
	print LOG "\nERROR:No gs-fisher results are detected in $infolders\n\n";
	exit;
}


#2,3,5,6
#number, gene, OddsRatio, bhp, 
my @outfiles_forheatmap;

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
push @outfiles_forheatmap,$outfile4;
open(OUT4,">$outfile4") || die $!;

my $outfile5=$outfile;
$outfile5=~s/\.\w+$/_p.txt/;
push @outfiles_forheatmap,$outfile5;
open(OUT5,">$outfile5") || die $!;

my $outfile6=$outfile;
$outfile6=~s/\.\w+$/_z.txt/;
open(OUT6,">$outfile6") || die $!;




print OUT1 "Gene set\t",join("\t",@usedcomparisons),"\n";
print OUT2 "Gene set\t",join("\t",@usedcomparisons),"\n";
print OUT3 "Gene set\t",join("\t",@usedcomparisons),"\n";
print OUT4 "Gene set\t",join("\t",@usedcomparisons),"\n";
print OUT5 "Gene set\t",join("\t",@usedcomparisons),"\n";
print OUT6 "Gene set\t",join("\t",@usedcomparisons),"\n";

foreach my $gs (sort keys %gss) {
	print OUT1 $gs,"\t";
	print OUT2 $gs,"\t";
	print OUT3 $gs,"\t";
	print OUT4 $gs,"\t";
	print OUT5 $gs,"\t";
	print OUT6 $gs,"\t";
	
	my (@marks1,@marks2,@marks3,@marks4,@marks5,@marks6);
	
	foreach my $file (@usedcomparisons) {
		if(defined $file2gs{$file}{$gs}) {
			push @marks1,$file2gs{$file}{$gs}[2]; #num
			push @marks2,$file2gs{$file}{$gs}[3]; #gene
			push @marks3,$file2gs{$file}{$gs}[6]; #or
			push @marks4,$file2gs{$file}{$gs}[7]; #q
			push @marks5,$file2gs{$file}{$gs}[5]; #p
			push @marks6,$file2gs{$file}{$gs}[4]; #z
		}
		else {
			push @marks1," ";
			push @marks2," ";
			push @marks3," ";
			push @marks4," ";		
			push @marks5," ";
			push @marks6," ";			
		}
	}
	print OUT1 join("\t",@marks1),"\n";
	print OUT2 join("\t",@marks2),"\n";
	print OUT3 join("\t",@marks3),"\n";
	print OUT4 join("\t",@marks4),"\n";	
	print OUT5 join("\t",@marks5),"\n";
	print OUT6 join("\t",@marks6),"\n";	
}
close OUT1;
close OUT2;
close OUT3;
close OUT4;
close OUT5;
close OUT6;

#Merge files into excel file
my $exceloutfile=$outfile;
$exceloutfile=~s/\.\w+$/\.xlsx/;

print STDERR "Converting outputs to excel file:$exceloutfile.\n\n";
print LOG "Converting outputs to excel file:$exceloutfile.\n\n";

system("$text2excel -i $outfile1,$outfile2,$outfile3,$outfile4,$outfile5,$outfile6 -n Number,Gene,OddsRatio,Q,P,Z -o $exceloutfile --theme theme2");



#draw heatmap and dotplots

#Generate heatmap figures for p/q

foreach my $outfile_sel (@outfiles_forheatmap) {
	my $outfilefig=$outfile_sel;
	$outfilefig=~s/\.\w+$/_heatmap.png/;	

	print STDERR "Generating heatmap for $outfile_sel.\n";
	print LOG "Generating heatmap for $outfile_sel.\n";

	system("$rscript $gs_heatmap -i $outfile_sel -o $outfilefig --top $topnum --sortby $sortby");
	print LOG "$rscript $gs_heatmap -i $outfile_sel -o $outfilefig --top $topnum --sortby $sortby\n";
	
}

#Generate dotplot for or+p and or+q
#OddsRatio+p

my $outfiledpfig1=$outfile;

#size by OR, color by bhp	
#$outfiledpfig1=~s/\.\w+$/_orbhp_dotplot.png/;
#system("$rscript $gs_dotplot --size $outfile3 --sizename OddsRatio --color $outfile4 --colorname \"'-Log10BHP'\" -o $outfiledpfig1 --top $topnum --sortby $sortby");
#print LOG "$rscript $gs_dotplot --size $outfile3 --sizename OddsRatio --color $outfile4 --colorname \"'-Log10BHP'\" -o $outfiledpfig1 --top $topnum --sortby $sortby\n";

#size by bhp, color by z
$outfiledpfig1=~s/\.\w+$/_z_bhp_dotplot.png/;
system("$rscript $gs_dotplot --size $outfile4 --sizename \"'-Log10BHP'\" --color $outfile6 --colorname Zscore -o $outfiledpfig1 --top $topnum --sortby $sortby");
print LOG "$rscript $gs_dotplot --size $outfile4 --sizename \"'-Log10BHP'\" --color $outfile6 --colorname Zscore -o $outfiledpfig1 --top $topnum --sortby $sortby\n";


my $outfiledpfig2=$outfile;


#size by OR, color by rawp	
#$outfiledpfig2=~s/\.\w+$/_orrawp_dotplot.png/;
#system("$rscript $gs_dotplot --size $outfile3 --sizename OddsRatio --color $outfile5 --colorname \"'-Log10RawP'\" -o $outfiledpfig2 --top $topnum --sortby $sortby");
#print LOG "$rscript $gs_dotplot --size $outfile3 --sizename OddsRatio --color $outfile5 --colorname \"'-Log10RawP'\" -o $outfiledpfig2 --top $topnum --sortby $sortby\n";

#size by bhp, color by z
$outfiledpfig2=~s/\.\w+$/_z_rawp_dotplot.png/;
system("$rscript $gs_dotplot --size $outfile5 --sizename \"'-Log10RawP'\" --color $outfile6 --colorname Zscore -o $outfiledpfig2 --top $topnum --sortby $sortby");
print LOG "$rscript $gs_dotplot --size $outfile5 --sizename \"'-Log10RawP'\" --color $outfile6 --colorname Zscore -o $outfiledpfig2 --top $topnum --sortby $sortby\n";


close LOG;

########
#Function
########
sub current_time {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $now = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $mday, $hour, $min, $sec);
	return $now;
}
