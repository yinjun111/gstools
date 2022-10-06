

library("argparser",quietly =T)

version="0.1"

#

description=paste0("gs_waterfallplot\nversion ",version,"\n","Usage:\nDescription: Generate dot image two data matrixes\n")


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

#topnum.p=max(min(topnum,sum(data.p>-log10(siglevel))),5)
#topnum.q=max(min(topnum,sum(data.q>-log10(siglevel))),5)


#convert data set to the top list
data.q.sign<-data.q*sign(data.z)

data.q.sign.order<-data.q.sign[order(data.q.sign,decreasing = T)]
data.z.order<-data.z[order(data.q.sign,decreasing = T)]
data.names<-rownames(data)[order(data.q.sign,decreasing = T)]

if(length(data.q.sign.order)>=topnum) {
	data.q.sign.order.sel<-data.q.sign.order[c(1:(topnum/2),(length(data.q.sign.order)-topnum/2+1):length(data.q.sign.order))]
	data.z.order.sel<-data.z.order[c(1:(topnum/2),(length(data.q.sign.order)-topnum/2+1):length(data.q.sign.order))]
	data.names.sel<-data.names[c(1:(topnum/2),(length(data.q.sign.order)-topnum/2+1):length(data.q.sign.order))]
} else {
	data.q.sign.order.sel<-data.q.sign.order
	data.z.order.sel<-data.z.order
	data.names.sel<-data.names
}



#####
#Function
#####

gs_waterfall<-function(z,p,gs,zname,pname,signame="p",siglevel=0.05) {
  
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
    #geom_hline(yintercept = -log10(siglevel), linetype="dashed",color="black") + 
    geom_segment(aes(y=log10(0.05),x=sum(p>0)+0.5,yend=log10(0.05),xend=0), linetype="dashed",color="black",size=0.3)+
    geom_segment(aes(y=-log10(0.05),x=sum(p>0)+0.5,yend=-log10(0.05),xend=length(p)+0.5), linetype="dashed",color="black",size=0.3)+
    annotate(geom="text", x=1, y=0, label=paste(signame,"<",siglevel,sep=""),color="black",size=3,hjust = 0) +
    theme(axis.text.x = element_text(angle = 0,size = 12 ),axis.text.y = element_text(angle = 0,size = 12 ))
  #panel.border = element_rect(colour = "black", fill=NA, size=1)
  
  print(plot1)  
}

#figure output, qvalue
figure1=sub(".png$","_bhp.png",args$out,perl=T)

cat("Generating ",figure1,"\n")

#auto calculate length

CairoPNG(file=figure1,res = 300,width = 3.5+floor(max(nchar(data.names.sel)/10)),height = topnum*0.3,units = "in")

gs_waterfall(z=data.z.order.sel,p=data.q.sign.order.sel,gs=data.names.sel,zname="ZScore",pname="-Log10BHP, by ZScore sign",siglevel=siglevel,signame="BHP")

dev.off()
