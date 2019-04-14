#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename);


########
#Prerequisites
########

my $r="/apps/R-3.4.1/bin/R";
my $rscript="/apps/R-3.4.1/bin/Rscript";

#Application version
my $mergefiles="/apps/sbptools/mergefiles/mergefiles_caller.pl";

#dev version
#my $mergefiles="perl /home/jyin/Projects/Pipeline/sbptools/mergefiles/mergefiles_caller.pl";

########
#Interface
########


my $version="0.1";


my $usage="

gstools gs-fisher
version: $version
Usage: gstools gs-fisher [parameters]

Description: Perform fisher's exact test for select genes for a list of gene sets(gs)


Mandatory Parameters:
    --in|-i           Input file
    --out|-o          Output file
    --anno|-a         Annotation file
    --backgroud|-b    Background file (optional)

    --qmethod         Multiple testing correction method [BH]
    --qcutoff         Corrected P cutoff [0.05]

    --runmode|-r      Where to run the scripts, local, server or none [none]
    --verbose|-v      Verbose
	
";


#R parameters
my $rparams="
  -i, --in IN			Expr input file
  -a, --anno ANNO			Annotation file
  -o, --out OUT			Output file
  -f, --formula FORMULA			DESeq formula
  -t, --treat TREAT			treatment name
  -r, --ref REF			reference name
  --fccutoff FCCUTOFF			Log2 FC cutoff [default: 1]
  -q, --qcutoff QCUTOFF			qcutoff [default: 0.05]
  -p, --pmethod PMETHOD			Method used in DESeq2 [default: Wald]
  --qmethod QMETHOD			FDR Method [default: BH]
  --filter FILTER			Count filter for all the exps [default: 10]
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
my $annofile;
my $bkfile;
my $outfile;

my $pmethod="fisher";
my $qmethod="BH";
my $alternative="greater"; #greater, less
my $qcutoff=0.05;
my $bkmodel="total"; #total, selected
my $sortby="pvalue"; #pvalue, depth, name

	
my $verbose=1;
my $runmode="none";

GetOptions(
	"input=s" => \$infile,
	"annotation=s" => \$annofile,
	#"pmethod=s"=>\$pmethod,
	"qmethod=s"=>\$qmethod,
	"qcutoff=f"=>$qcutoff,
	"alternative=s"=>\$alternative,
	"background=s"=>\$bkfile,
	"sort=s"=> \$sortby,
	"out|o=s"=>\$outfile,
	"verbose"=>\$verbose,	
);


$outfile = abs_path($outfile);


my $logfile=$outfile;
$logfile=~s/.txt/_log.txt/;


#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();

print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";

print LOG "\n";

print STDERR "\nsbptools gs-fisher $version running ...\n\n" if $verbose;
print LOG "\nsbptools gs-fisher $version running ...\n\n";




########
#Process
########


my $filefolder="temp";

my $rand=randstring();


my %selgenes;
my %gene_id;
my %gs2gene;
my %gene_anno;

#####
#BK file
#####

my $blinenum=0;
my %bkgenes;
if(length($bkfile)>0) {
	open(IN,$bkfile) || die "ERROR:Can't read $bkfile.$!\n";
	
	while(<IN>) {
		tr/\r\n//d;
		unless($blinenum==0) {
			my @array=split/\t/;
			$bkgenes{$array[0]}++;
			$gene_anno{$array[0]}++;
		}
		$blinenum++;
	}
	close IN;
	
	print STDERR scalar(keys %gene_anno)," features found in the $bkfile file.\n" if $verbose;
}
else {
	print STDERR "No background file defined.\n" if $verbose;
}
	
#########
#Annotation File
#########

my $alinenum=0;
open(ANNO,$annofile) || die "ERROR:Can't read $annofile.$!\n";
while(<ANNO>) {
	tr/\r\n\"\'//d;
	
	unless($alinenum==0) {
		my @array=split/\t/;
		
		if(length($bkfile)>0) {
			unless(defined $bkgenes{$array[0]}) {
				next;
			}
		}
		
		#genes in the annotation
		$gene_anno{$array[0]}++;
		
		if(defined $array[1] && $array[1] ne " ") {
	
			foreach my $gs (split(";",$array[1])) {
				if($gs eq "-" || $gs eq " ") {
					next;
				}
				
				$gs2gene{$gs}{$array[0]}++;
			}
		}

	}
	
	$alinenum++;
}
close ANNO;


if(length($bkfile)==0) {
	%bkgenes=%gene_anno; #use all genes in the anno as bk
}



######
#Input
######

#support multiple columns of genes

open(IN,$infile) || die "ERROR:Can't read $infile.$!\n";
my @contents=<IN>;
my $title=shift @contents;
$title=~tr/\r\n//d;
my @comparison_names=split(/\t/,$title);
print STDERR "Sample list(s) ",join(",",@comparison_names)," are identified from $infile.\n" if $verbose;


foreach (@contents) {
	tr/\r\n\"//d;
	my @array=split/\t/;
	my $num=-1;
	foreach my $id (@array) {
		$num++;
		
		if(defined $id && length($id)>0 && $id ne " ") {
			if(defined $gene_anno{$id}) {
				$selgenes{$comparison_names[$num]}{$id}++;
			}
		}
	}
}



#######
#Calculation
#######

my $all_num=keys %bkgenes;

my $outfolder=abs_path(dirname($outfile));

if(!-e "$outfolder/$filefolder") {
	mkdir("$outfolder/$filefolder");
}

my $tempfile="$outfolder/$filefolder/gs-fisher_tem_$rand.txt";

#if(@comparison_names<=1) {
#	my $exp_id=$comparison_names[0];
#}
#else {
	#should only use this after DE_sigs
	
	foreach my $exp_id (@comparison_names) {
		open(TMP,">$tempfile") || die $!;
		my %gs2gene_sel; 
		my %changed_gene; #n
		
		foreach my $gs (sort keys %gs2gene) {
			foreach my $id (sort keys %{$selgenes{$exp_id}}) {
				next unless defined $gene_anno{$id};
				
				$changed_gene{$id}++;
				
				if(defined $gs2gene{$gs}{$id}) {
					$gs2gene_sel{$gs}{$id}++;
				}
			}
		}
		
		my $total_num=keys %changed_gene; #n
		
		print TMP "Total\t$all_num\t$total_num\t-\n"; #N, n,
		
		my $k_num;
		
		foreach my $gs (sort keys %gs2gene) {
			if(defined $gs2gene_sel{$gs}) {
				my $num2=keys %{$gs2gene_sel{$gs}}; #m
				my $num1=keys %{$gs2gene{$gs}}; #M
				
				print TMP $gs,"\t",$num1,"\t",$num2,"\t",join(";",sort keys %{$gs2gene_sel{$gs}}),"\n"; #M, m
			
			}
		}
		close TMP;
		
		my $newid=$exp_id;
		$newid=~s/\W/_/g;
		
		my $file_out=substr($outfile,0,length($outfile)-4)."_$newid.txt";
		
		print STDERR "Result File,$file_out\t" if $verbose;
		
		path_r($tempfile,$file_out);
		
		print STDERR "done!\n" if $verbose;
		
		##unlink($tempfile); #masked for test
	}
#}



			
###################################
#R scripts here
###################################
sub max {
	my @nums=@_;
	my $max=$nums[0];
	foreach (@nums) {
		if($_>$max) {
			$max=$_;
		}
	}
	return $max;
}

sub path_r {
	my ($file_in,$file_out)=@_;
	open(RPROG,"|$r --no-restore --no-save --slave") || die $!;

	select RPROG;
	print <<"CODE";
	
###This script reads the pathway list from microarray result.
###The P value is calculated by Hypergeometric distribution, a package in R.
###The Q value is calculated by a R package integrated below, which was from CRAN.
###By default, FDR is set to 0.05.

path_stat <- function(file_in,file_out,pmethod=c("hyper","fisher"),qmethod="BH",Q=0.05,alternative="greater",sortby="pvalue") {
	table<-read.table(file_in,header=FALSE,sep="\t") #row.names=1
	m<-table[1,3]	#n
	n<-table[1,2]-table[1,3]	#N-n
	k<-table[,2]	#M
	x<-table[,3]	#m
	num<-2;
	r_enrich<-c()
	
	#fisher
	if(pmethod=="fisher") {
		p_hyper<-c()
		while(num <= length(x)) {
			data<-matrix(c(x[num],m-x[num],k[num]-x[num],n-k[num]+x[num]),nr=2)
			p_fisher<-fisher.test(data,alternative=alternative)
			
			p_hyper[num]<-p_fisher\$p.value

			if(is.infinite(p_fisher\$estimate)) {
				#an aritificial OR to replace inf
				r_enrich[num]<-(x[num]+0.1)/(k[num]-x[num]+0.1)/(m-x[num]+0.1)*(n-k[num+x[num]+0.1)
			}
			else {
				r_enrich[num]<-p_fisher\$estimate
			}
			
			num<-num+1
		}
		realp<-p_hyper[2:length(p_hyper)]
		realr<-r_enrich[2:length(r_enrich)]
	}

	#Enrichment factor
	#num<-2
	#r_enrich<-c()
	#while(num <= length(x)) {
	#	r_enrich[num]<- x[num]*(n+m)/(k[num]*m)
	#	num<-num+1
	#}
	#realr<-r_enrich[2:length(r_enrich)]
	
	
    if(qmethod == "BH") {
		qvalue<-p.adjust(realp,method="BH");
		sig<-qvalue
		sig[qvalue<Q]="TRUE";
		sig[qvalue>=Q]="FALSE";
		#qvalue<-c("-",qvalue)
		#sig<-c('-',sig);
	} else {
		qvalue<-p.adjust(realp,method=qmethod);
		sig<-qvalue
		sig[qvalue<Q]="TRUE";
		sig[qvalue>=Q]="FALSE";
		#qvalue<-c("-",qvalue)
		#sig<-c('-',sig);	
	}
	
    #pvalue<-c('-',realp)
	#rvalue<-c('-',realr)
	
	totalline<-c("Total",table[1,2],table[1,3],"-","-","-","-","-")
	#content<-data.frame(cbind(table,pvalue,rvalue,qvalue,sig))
	content<-cbind(as.matrix(table[2:nrow(table),]),realp,realr,qvalue,sig)
	
	#sort function
	if(nrow(content)>1) {
		if(sortby=="pvalue") {
			content<-content[order(realp),]
		} else if(sortby=="name"){
			content<-content[order(content[,1]),]
		}
	}
	
	result<-rbind(totalline,as.matrix(content))
	
	colnames(result)<-c("Gene set",'No. of features in the gene set','No. of selected features in the gene set','Feature ID',paste(c('p_value',pmethod,alternative),collapse="_"),'Odds Ratio',paste(c('q_value',qmethod),collapse="_"),paste(c("significance",Q),collapse="_"))

	write.table(data.frame(content),file=file_out,sep="\\t",row.names=FALSE,quote=F)
	#cat(log)
	#pathout<-"path_stat.log"
	#write(log,file=pathout,append = TRUE)
}

path_stat("$file_in","$file_out",Q=$qcutoff,qmethod="$qmethod",pmethod="$pmethod",alternative="$alternative",sortby="$sortby");

#warnings()

q()
CODE
close RPROG;

}

sub randstring {
	my @chars=("A".."Z");
	my $string;
	
	$string.=$chars[rand(@chars)] for 1..6;
	
	return $string;
}


