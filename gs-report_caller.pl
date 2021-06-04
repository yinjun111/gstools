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


my $version="0.2";

#v0.2, changed to matrix input
#v0.21, keep original order in input


my $usage="

gs-report
version: $version
Usage: gstools gs-report [parameters]

Description: Read rnaseq-summary folder and generate gs-fisher results and summary reports.

Mandatory Parameters:
    --in|-i           a matrix with each row as a gene and each 
                         column as a comparison. 1 and -1 as up/down-regulated genes
                         E.g. rnaseq-summary_folder/rnaseq-summary_GeneDESigs.txt

    --tx|-t           Transcriptome
                        Current support Human.B38.Ensembl88, Mouse.B38.Ensembl88

    --config|-c       Configuration file for gstools-db, by default
                        Human: /data/jyin/Databases/gstools-db/gstools-db-config_human.txt
                        Mouse: /data/jyin/Databases/gstools-db/gstools-db-config_mouse.txt

    --comparisons     (Optional) Name of comparisons to be tested in gs analyses
                         By default, all the comparisons in the rnaseq-summary folder

    --top             No. of top terms to show in figures [20]
    --sortby          Column name used to sort the top terms, or averge of all columns [avg]
    --qcutoff         Corrected P cutoff [0.05]
	
    --out|-o          Output folder
	
    --runmode|-r      Where to run the scripts, local, cluster or none [none]
	
	
    #Parallel computing controls
    --task            Number of tasks to be paralleled. By default 20 tasks. [20]
    --ncpus           No. of cpus for each task [2]
    --mem|-m          Memory usage for each process, e.g. 100mb, 100gb
	
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
	
my $inputfile;
my $tx;
my $configfile;
my $comparisons;
my $outputfolder;
my $topnum=20;
my $sortby="avg";
my $qcutoff=0.05;
my $dev=0;

my $task=20;
my $ncpus=2;
my $mem;
my $runmode="none";
	
my $verbose=1;


GetOptions(
	"input|i=s" => \$inputfile,
	"tx|t=s" => \$tx,
	"config|c=s" => \$configfile,
	"comparisons=s" => \$comparisons,
	"top=s"=>\$topnum,	
	"sortby=s" => \$sortby,
	"qcutoff=f"=>\$qcutoff,
	"out|o=s"=>\$outputfolder,

	"task=s" => \$task,
	"ncpus=s" => \$ncpus,
	"mem=s" => \$mem,	

	"runmode|r=s" => \$runmode,			
	"dev" => \$dev,	
);


###
#Tools needed
####


my $omictoolsfolder="/apps/omictools/";

#adding --dev switch for better development process
#if($dev) {
#	$omictoolsfolder="/home/jyin/Projects/Pipeline/omictools/";
#}
#else {
#	#the tools called will be within the same folder of the script
#	$omictoolsfolder=get_parent_folder(abs_path(dirname($0)));
#}

my $mergefiles="$omictoolsfolder/mergefiles/mergefiles_caller.pl";
my $parallel_job="$omictoolsfolder/parallel-job/parallel-job_caller.pl";
my $text2excel="$omictoolsfolder/text2excel/text2excel.pl";


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

$inputfile = abs_path($inputfile);

my $scriptlocalrun="$outputfolder/gs-report_local_submission.sh";
my $scriptclusterrun="$outputfolder/gs-report_cluster_submission.sh";

my $scriptfolder="$outputfolder/scripts";
if(!-e $scriptfolder) {
	mkdir($scriptfolder);
}

my $scriptfile1="$scriptfolder/gs-report_run.sh";

my $logfile="$outputfolder/gs-report_run.log";


#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();
my $timestamp=build_timestamp($now,"long");

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

if(!-e $inputfile) {
	print STDERR "ERROR:Can't read $inputfile. $!\n";
	print LOG "ERROR:Can't read $inputfile. $!\n";
	exit;
}

print STDERR "$inputfile is analyzed by gs-report.\n\n";
print LOG "$inputfile is analyzed by gs-report.\n\n";


#generate new DE sigs based on --comparisons
#pre-selected comparisons
my @selcomparisons;
my @selcomparisons_updown;
my %selcomparisons_hash;
if(defined $comparisons && length($comparisons)>0) { 
	#keep original order
	@selcomparisons=split(",",$comparisons);
	#@selcomparisons_updown=((map {$_."_up"} @selcomparisons),(map {$_."_down"} @selcomparisons));
	%selcomparisons_hash=map {$_,1} @selcomparisons;
}

my $newsigfile="$outputfolder/gs-report_inputs.txt";

my %comp2col;
my @selcols;
my @usedcomparisons;

if(@selcomparisons) {
	open(OUT,">$newsigfile") || die $!;
	open(IN,"$inputfile") || die $!;
	while(<IN>) {
		tr/\r\n//d;
		my @array=split/\t/;
		
		if($_=~/^Gene/) {
			for(my $num=1;$num<@array;$num++) {
				$comp2col{$array[$num]}=$num;
			}
			
			foreach my $comp (@selcomparisons) {
				if(defined $comp2col{$comp}) {
					push @selcols,$comp2col{$comp};
					print STDERR $comp," is identifed at $comp2col{$comp}(th) column.\n";
					print LOG $comp," is identifed at $comp2col{$comp}(th) column.\n";
				}
				else {
					print STDERR "ERROR:$comp is not defined in $inputfile\n\n";
					print LOG "ERROR:$comp is not defined in $inputfile\n\n";
					exit;
				}
			}
			
			print STDERR "Columns ", join(",",@selcols)," are used by --comparisons $comparisons.\n";
			print LOG "Columns ", join(",",@selcols)," are used by --comparisons $comparisons.\n";			
			
		}
		
		print OUT join("\t",@array[0,@selcols]),"\n";
	}
	close IN;
	close OUT;

}
else {
	system("cp $inputfile $newsigfile");
	
	#keep original order
	open(IN,"$inputfile") || die $!;
	while(<IN>) {
		tr/\r\n//d;
		my @array=split/\t/;
		
		if($_=~/^Gene/) {
			@usedcomparisons=@array[1..$#array];
			last;
		}
	}
	close IN;
	
	print STDERR "No --comparisons defined. Use existing order: ",join(",",@usedcomparisons),"\n";
	print LOG "No --comparisons defined. Use existing order: ",join(",",@usedcomparisons),"\n";
	
}



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
	print S1 "perl $gs_fisher -i $newsigfile -t matrix -c both -o $outputfolder/$dbfiles{$dbfile} -d $dbfile --top $topnum --siglevel $qcutoff -a ",$tx2ref{$tx}{"geneanno"},print_dev($dev,";");
	
	if(@selcomparisons) {
		print S1 "perl $gs_fisher_summary -i $outputfolder/$dbfiles{$dbfile} -o $outputfolder/$dbfiles{$dbfile}_summary --comparisons $comparisons --top $topnum --sortby $sortby",print_dev($dev,";\n");
	}
	else {
		print S1 "perl $gs_fisher_summary -i $outputfolder/$dbfiles{$dbfile} -o $outputfolder/$dbfiles{$dbfile}_summary --comparisons ",join(",",@usedcomparisons)," --top $topnum --sortby $sortby",print_dev($dev,";\n");
	}
}
close S1;




#######
#Run mode
#######

open(LOUT,">$scriptlocalrun") || die "ERROR:can't write to $scriptlocalrun. $!";
open(SOUT,">$scriptclusterrun") || die "ERROR:can't write to $scriptclusterrun. $!";


my @scripts_all=($scriptfile1);


#print out command for local and cluster parallel runs
my $jobnumber=0;
my $jobname="gs-report-$timestamp";

if($task eq "auto") {
	$jobnumber=0;
}
else {
	$jobnumber=$task;
}

my @local_runs;
my @script_names;

foreach my $script (@scripts_all) {
	push @local_runs,"cat $script | parallel -j $jobnumber";

	if($script=~/([^\/]+)\.\w+$/) {
		#push @script_names,$1."_".basename_short($outputfolder);
		push @script_names,$1;
	}
}

my $localcommand="screen -S $jobname -dm bash -c \"source ~/.bashrc;".join(";",@local_runs).";\"";

print LOUT $localcommand,"\n";
close LOUT;

#print out command for cluster parallel runs

my $clustercommand="perl $parallel_job -i ".join(",", @scripts_all)." -o $scriptfolder -n ".join(",",@script_names)." --tandem -t $task --ncpus $ncpus -r ";

if(defined $mem && length($mem)>0) {
	$clustercommand.=" -m $mem";
}

print SOUT $clustercommand,"\n";
close SOUT;



if($runmode eq "none") {
	print STDERR "\nTo run locally, in shell type: sh $scriptlocalrun\n";
	print STDERR "To run in cluster, in shell type: sh $scriptclusterrun\n";
	
	print LOG "\nTo run locally, in shell type: sh $scriptlocalrun\n";
	print LOG "To run in cluster, in shell type: sh $scriptclusterrun\n";
}
elsif($runmode eq "local") {
	#local mode
	#implemented for Falco
	
	system("sh $scriptlocalrun");
	print LOG "sh $scriptlocalrun;\n\n";

	print STDERR "Starting local paralleled processing using $jobnumber tasks. To monitor process, use \"screen -r $jobname\".\n\n";
	print LOG "Starting local paralleled processing using $jobnumber tasks. To monitor process, use \"screen -r $jobname\".\n\n";
	
}
elsif($runmode eq "cluster") {
	#cluster mode
	#implement for Firefly
	
	system("sh $scriptclusterrun");
	print LOG "sh $scriptclusterrun;\n\n";

	print STDERR "Starting cluster paralleled processing using $jobnumber tasks. To monitor process, use \"qstat\".\n\n";

}

close LOG;



########
#Functions
########

sub current_time {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $now = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $mday, $hour, $min, $sec);
	return $now;
}

sub build_timestamp {
	my ($now,$opt)=@_;
	
	if($opt eq "long") {
		$now=~tr/ /_/;
		$now=~tr/://d;
	}
	else {
		$now=substr($now,0,10);
	}
	
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