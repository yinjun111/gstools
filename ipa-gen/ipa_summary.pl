#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);


#get an input folder 


########
#Interface
########


my $version="0.3";

#v0.2, add UR support
#v0.3, dotplot/heatmap


my $usage="

ipa-summary
version: $version
Usage: gstools ipa-summary [parameters]

Description: Merge multiple IPA results


Mandatory Parameters:
    --in|-i           Input folder(s)					  
    --out|-o          Output folder
    --db|-d           IPA database CP or UR [CP]
    --nochem          Whether to remove IPA UR Chems, default False [F]
    --sigonly|-s      Only significant results [T]	
    --suffix|-x       Remove suffix from file names [_CP.txt]	
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

my $db="CP";
my $nochem="F";
my $sigonly="T";
my $topnum=20;
my $sortby="avg";
my $comparisons;
my $suffix="_CP.txt";

my $verbose=1;
my $runmode="none";

my $pcutoff=0.05;

my $dev=0;

GetOptions(
	"in|i=s" => \$infolders,
	"out|o=s"=>\$outfolder,
	"db|d=s"=>\$db,
	"nochem=s"=>\$nochem,	
	"comparisons=s" => \$comparisons,	
	"sigonly|s=s"=>\$sigonly,
	"suffix|x=s"=>\$suffix,		
	"sortby=s" => \$sortby,
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



##################

#mkpath for recursive make dir

print STDERR "\nipa-summary $version running ...\n\n" if $verbose;

if(!-e $outfolder) {
	print STDERR "$outfolder not found. mkdir $outfolder\n";
	mkdir($outfolder); #doesn't do recursive
}

$outfolder=abs_path($outfolder);
my $outfoldername=basename($outfolder);

#the order of dirname and abs_path is important
my $outfile = abs_path($outfolder)."/".$outfoldername.".xlsx";


my $logfile=$outfile;
$logfile=~s/\.\w+$/_ipa_summary.log/;


#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();

print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";

print LOG "\n";


print LOG "\nipa-summary $version running ...\n\n";




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




#z,BHP, gene, number

my $pcol;
my $zcol;
my $genecol;

if($db eq "CP") {
	$pcol=1;
	$zcol=3;
	$genecol=4;
}
#elsif($db eq "UR") {
#	$pcol=6;
#	$zcol=4;
#	$genecol=7;
#}
#else {
#	print STDERR "ERROR:--db $db not supported.\n\n";
#}


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
			
			#print STDERR $filename,"\n";

			#$files{$filename}++;
			$files{$samplename}++;
			print STDERR $file," detected from ipa.\n";
			print LOG $file," detected from ipa.\n";
			
			#check --comparisons
			if(@selcomparisons) {
				unless(defined $selcomparisons_hash{$samplename}) {
					print STDERR "$samplename is not defined in --comparisons $comparisons. Skip...\n";
					print LOG "$samplename is not defined in --comparisons $comparisons. Skip...\n";
					next;
				}
			}
			
			my $linenum=0;
			
			my %title2colnum;
			
			while(<IN>) {
				tr/\r\n//d;
				#skip copyright and title
				
				my @array=split/\t/;
				
				if($db eq "UR" && $linenum == 2) {
					#Expr Log Ratio	Molecule Type	Predicted Activation State	Activation z-score	Flags	p-value of overlap	Target Molecules in Dataset
					for(my $num=0;$num<@array;$num++) {
						my $term=$array[$num];
						if($term eq "p-value of overlap") {
							$pcol=$num;
						}
						
						if($term eq "Activation z-score") {
							$zcol=$num;
						}

						if($term eq "Target Molecules in Dataset") {
							$genecol=$num;
						}
					}
					
					print STDERR $filename,"\t",join(",",$pcol,$zcol,$genecol),"\n";
				}
				
				if($linenum>=3) {					
					
					$file2gs{$filename}{$array[0]}=[@array[$pcol,$zcol,$genecol]];
					
					if($db eq "CP") {
						if($sigonly eq "T") {
							#keeps only significant results
							#if($array[7] eq "TRUE") {
							
							if(10**-$array[$pcol] <$pcutoff) {
								$gss{$array[0]}++;
							}
						}
						else {
							$gss{$array[0]}++;
						}
					}
					elsif($db eq "UR") {
						if($sigonly eq "T") {
							
							if($array[$pcol] >$pcutoff) {
								next;
							}
						}
						
						if($nochem eq "T") {
							if($array[2]=~/^chem/) {
								next;
							}
						}

						$gss{$array[0]}++;
					}
				}
				
				$linenum++;
			}
			close IN;
		}
	}
}

#2,4,5
#p, z, genes
my $outfile1=$outfile;
$outfile1=~s/\.\w+$/_p.txt/;
open(OUT1,">$outfile1") || die $!;

my $outfile2=$outfile;
$outfile2=~s/\.\w+$/_z.txt/;
open(OUT2,">$outfile2") || die $!;

my $outfile3=$outfile;
$outfile3=~s/\.\w+$/_genes.txt/;
open(OUT3,">$outfile3") || die $!;

print OUT1 "Gene set\t",join("\t",sort keys %files),"\n";
print OUT2 "Gene set\t",join("\t",sort keys %files),"\n";
print OUT3 "Gene set\t",join("\t",sort keys %files),"\n";

foreach my $gs (sort keys %gss) {
	print OUT1 $gs,"\t";
	print OUT2 $gs,"\t";
	print OUT3 $gs,"\t";
	
	#p, z, gene
	my (@marks1,@marks2,@marks3);
	
	foreach my $file (sort keys %file2gs) {
		if(defined $file2gs{$file}{$gs}) {
			if($db eq "CP") {
				push @marks1,10**-$file2gs{$file}{$gs}[0]; #conversion of p value
			}
			else {
				push @marks1,$file2gs{$file}{$gs}[0]; #conversion of p value			
			}
			
			push @marks2,$file2gs{$file}{$gs}[1];
			push @marks3,$file2gs{$file}{$gs}[2];
		}
		else {
			push @marks1," ";
			push @marks2," ";
			push @marks3," ";		
		}
	}
	print OUT1 join("\t",@marks1),"\n";
	print OUT2 join("\t",@marks2),"\n";
	print OUT3 join("\t",@marks3),"\n";
	
}
close OUT1;
close OUT2;
close OUT3;



#Merge files into excel file
my $exceloutfile=$outfile;
$exceloutfile=~s/\.\w+$/\.xlsx/;

print STDERR "Converting outputs to excel file:$exceloutfile.\n\n";
print LOG "Converting outputs to excel file:$exceloutfile.\n\n";

system("$text2excel -i $outfile1,$outfile2,$outfile3 -n P,Z,Genes -o $exceloutfile --theme theme2");



my $outfiledpfig1=$outfile;


#size by bhp, color by z
$outfiledpfig1=~s/\.\w+$/_z_p_dotplot.png/;
system("$rscript $gs_dotplot --size $outfile1 --sizename \"'-Log10P'\" --color $outfile2 --colorname Zscore -o $outfiledpfig1 --top $topnum --sortby $sortby");
print LOG "$rscript $gs_dotplot --size $outfile1 --sizename \"'-Log10P'\" --color $outfile2 --colorname Zscore -o $outfiledpfig1 --top $topnum --sortby $sortby\n";


close LOG;

########
#Function
########
sub current_time {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $now = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $mday, $hour, $min, $sec);
	return $now;
}
