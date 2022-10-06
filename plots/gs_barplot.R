

library("argparser",quietly =T)

version="0.2"

#v0.2, change colors

description=paste0("gs_barplot\nversion ",version,"\n","Usage:\nDescription: Generate dot image two data matrixes\n")


#####
#Add arguments
#####

parser <- arg_parser(description=description)


parser <- add_argument(parser, arg="--in",short="-i", type="character", help = "Input file from gs-fisher")
parser <- add_argument(parser, arg="--out",short="-o", type="character", help = "Output file, png format")
parser <- add_argument(parser, arg="--top",short="-t", type="character", default ="20",help = "# of top terms to show")
parser <- add_argument(parser, arg="--siglevel",short="-s", type="character", default ="0.05",help = "significance level")


#need to add size/color transformation options

args = parse_args(parser)

#print(args)

######
#Drawing function
#####

library(Cairo)
#library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(scales)

#cols <- brewer.pal(9, "Reds")
#pal <- colorRampPalette(cols)


siglevel<-as.numeric(args$siglevel)

data<-read.table(args$"in",header=T,row.names=1,sep="\t",quote="",comment.char="")

#remove first row
data<-data[-1,]

rownames(data)<-make.names(substr(rownames(data),1,50),unique=T)


#values
data.z<-as.numeric(as.vector(unlist(data[,4])))
data.p<-as.numeric(as.vector(unlist(data[,5])))
data.p[is.na(data.p)]<-1
data.p<- -log10(data.p+1e-32) #conversion
data.q<-as.numeric(as.vector(unlist(data[,7])))
data.q[is.na(data.q)]<-1
data.q<- -log10(data.q+1e-32) #conversion

#decide top number. at least 5 are shown
if(args$top =="all") {
	topnum=nrow(data)
} else {
	topnum<-as.numeric(args$top)
}

topnum.p=max(min(topnum,sum(data.p>-log10(siglevel))),5)
topnum.q=max(min(topnum,sum(data.q>-log10(siglevel))),5)


#####
#Function
#####

gs_bar<-function(z,p,gs,zname,pname,signame="p",siglevel=0.05) {
  
  #create data frame
  
  gs<-make.names(gs,unique=T)
  
  data.df<-data.frame(factor(gs,levels=rev(gs)),		
                      as.numeric(unlist(z)),
                      as.numeric(unlist(p)))
  
  colnames(data.df)<-c("Gene Set",zname, pname)			  
  
  zname.rev<-paste("`",zname,"`",sep="")
  pname.rev<-paste("`",pname,"`",sep="")
  
  plot1<-ggplot(data.df, aes_string(x="`Gene Set`", y=pname.rev, fill=zname.rev)) + 
    geom_bar(stat='identity',width=0.9)+
    theme_classic() +
    scale_fill_gradient2(low = "#0571b0",  mid="grey",high = "#ca0020", space = "Lab", limit = c(-4, 4),oob=squish) + 
    coord_flip()+
    geom_hline(yintercept = -log10(siglevel), linetype="dashed",color="black") + 
    annotate(geom="text", x=1, y=-log10(0.05), label=paste(signame,"<",siglevel,sep=""),color="black",size=4,hjust = 0) +
    theme(axis.text.x = element_text(angle = 0,size = 12 ),axis.text.y = element_text(angle = 0,size = 12 ))
  #panel.border = element_rect(colour = "black", fill=NA, size=1)
  
  print(plot1)  
}

#figure output, qvalue
figure1=sub(".png$","_bhp.png",args$out,perl=T)

cat("Generating ",figure1,"\n")

#auto calculate length

CairoPNG(file=figure1,res = 300,width = 3.5+floor(max(nchar(rownames(data)[1:topnum.q])/10)),height = topnum.q*0.3,units = "in")

gs_bar(z=data.z[1:topnum.q],p=data.q[1:topnum.q],gs=rownames(data)[1:topnum.q],zname="ZScore",pname="-Log10BHP",siglevel=siglevel,signame="BHP")

dev.off()


#figure output, pvalue

figure2=sub(".png$","_rawp.png",args$out,perl=T)

cat("Generating ",figure2,"\n")

CairoPNG(file=figure2,res = 300,width = 3.5+floor(max(nchar(rownames(data)[1:topnum.p])/10)),height = topnum.p*0.3,units = "in")

gs_bar(z=data.z[1:topnum.p],p=data.p[1:topnum.p],gs=rownames(data)[1:topnum.p],zname="ZScore",pname="-Log10RawP",siglevel=siglevel,signame="P")

dev.off()