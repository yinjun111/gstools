

library("argparser",quietly =T)

version="0.2"

#v0.2, make z+p/q as default

description=paste0("gs_dotplot\nversion ",version,"\n","Usage:\nDescription: Generate dot image two data matrixes\n")


#####
#Add arguments
#####

parser <- arg_parser(description=description)

parser <- add_argument(parser, arg="--size",short="-s", type="character", help = " Input file for dot size, e.g. P/Q value")
parser <- add_argument(parser, arg="--sizename", type="character", help = "Name for dot size, e.g. -Log10(BHP)")
parser <- add_argument(parser, arg="--color",short="-c", type="character", help = "Input file for dot color, e.g. Zscore")
parser <- add_argument(parser, arg="--colorname",type="character", help = "Name for dot color, e.g. Zscore")
parser <- add_argument(parser, arg="--out",short="-o", type="character", help = "Output file, png format")
parser <- add_argument(parser, arg="--top",short="-t", type="character", help = "# of top terms to show")
parser <- add_argument(parser, arg="--sortby",short="-y", type="character",default ="avg", help = "Sort by which column")
parser <- add_argument(parser, arg="--colortransform",short="-c", type="character",default ="none", help = "Transformation for color")
parser <- add_argument(parser, arg="--sizetransform",short="-c", type="character",default ="mLog10", help = "Transformation for size")
parser <- add_argument(parser, arg="--colorna",short="-c", type="character",default ="0", help = "NA conversion for color")
parser <- add_argument(parser, arg="--sizena",short="-c", type="character",default ="1", help = "NA conversion for size")


#need to add size/color transformation options

args = parse_args(parser)

print(args)


######
#Drawing function
#####

library(Cairo)
#library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(scales)

cols <- brewer.pal(9, "Reds")
pal <- colorRampPalette(cols)

data.size<-read.table(args$"size",header=T,row.names=1,sep="\t",quote="",comment.char="")
data.color<-read.table(args$"color",header=T,row.names=1,sep="\t",quote="",comment.char="")

#remove non-word chars
colnames(data.size)==gsub("\\W","_",colnames(data.size))
colnames(data.color)==gsub("\\W","_",colnames(data.color))
args$sortby=gsub("\\W","_",args$sortby)

#conversion for size #p/q

#na conversion
data.size[is.na(data.size)]<-as.numeric(args$"sizena")
data.color[is.na(data.color)]<-as.numeric(args$"colorna")


#conversion for size #p
if(args$"sizetransform"=="mLog10") {
	data.size[data.size==0] <- 1e-32
	data.size<- -log10(data.size)
} else if(args$"sizetransform"=="none") {
	data.size<-data.size
}

#conversion for color #z
if(args$"colortransform"=="mLog10") {
	data.color[data.color==0] <- 1e-32
	data.color<- -log10(data.color)
} else if(args$"colortransform"=="none") {
	data.color<-data.color
}

#top terms #by size
topnum=as.numeric(args$top)
topnum=min(topnum,nrow(data.size))

#select by data size
if(args$sortby=="avg") {
	data.color.top<-data.color[order(apply(data.size,1,sum),decreasing = T)[1:topnum],]
	data.size.top<-data.size[order(apply(data.size,1,sum),decreasing = T)[1:topnum],]
} else {
	data.color.top<-data.color[order(data.size[,args$sortby],decreasing = T)[1:topnum],]
	data.size.top<-data.size[order(data.size[,args$sortby],decreasing = T)[1:topnum],]
}

rownames(data.size.top)<-make.names(substr(rownames(data.size.top),1,50),unique=T)


#df
#sizename<-make.names(args$sizename)
#colorname<-make.names(args$colorname)
sizename<-args$sizename
colorname<-args$colorname

#for -Log10(BHP)...
sizename.rev<- paste("`",sizename,"`",sep="")

#create data frame
data.df<-data.frame(factor(rep(colnames(data.color),each=topnum),levels=colnames(data.color)),
    factor(rep(rownames(data.color.top),ncol(data.color.top)),levels=rev(rownames(data.color.top))),
	as.numeric(unlist(data.size.top)),
	as.numeric(unlist(data.color.top)))

#filter OR==0
data.df<-data.df[data.df[,3]!=0,]

colnames(data.df)<-c("Comparison","Gene Set", sizename, colorname)			  

#figure output

#unclustered version
figure1<-args$out

#auto calculate length

CairoPNG(file=figure1,res = 300,width = 2+floor(max(nchar(rownames(data.color.top))/10))+ncol(data.color)*0.5,height = 3+0.4*topnum,units = "in")

#plot1<-ggplot(data.df, aes_string(x="Comparison", y="`Gene Set`" , size=sizename.rev, color=colorname)) + geom_point(alpha = 1)+theme_classic() +scale_color_gradient2(low = "blue",  mid="grey",high = "red", space = "Lab", limit = c(0, max(data.df[[colorname]])))+scale_size(range = c(0, 10))+ theme(axis.text.x = element_text(angle = 90,size = 10 ),axis.text.y = element_text(angle = 0,size = 10 ))
plot1<-ggplot(data.df, aes_string(x="Comparison", y="`Gene Set`" , size=sizename.rev, color=colorname)) + geom_point(alpha = 1)+theme_classic() +scale_color_gradient2(low = "blue",  mid="grey",high = "red", space = "Lab", limit = c(-4, 4),oob=squish)+scale_size(range = c(0, 10))+ theme(axis.text.x = element_text(angle = 90,size = 10 ),axis.text.y = element_text(angle = 0,size = 10 ))

print(plot1)

dev.off()

#pvalue only version
#figure2=sub(".png$","_ponly.png",figure1,perl=T)

#CairoPNG(file=figure2,res = 300,width = 2+floor(max(nchar(rownames(data.color.top))/10))+ncol(data.color)*0.5,height = 3+0.4*topnum,units = "in")

#plot1<-ggplot(data.df, aes_string(x="Comparison", y="`Gene Set`" , size=sizename.rev, color=sizename.rev)) + geom_point(alpha = 1)+theme_classic() +scale_color_gradient2(low = "blue",  mid="grey",high = "red", space = "Lab", limit = c(0, max(data.df[[colorname]])))+scale_size(range = c(0, 10))+ theme(axis.text.x = element_text(angle = 90,size = 10 ),axis.text.y = element_text(angle = 0,size = 10 ))

#print(plot1)

#dev.off()



#clustered version
#figure3=sub(".png$","_clustered.png",figure1,perl=T)

#CairoPNG(file=figure2,res = 300,width = 15,height = 10,units = "in")
#gs_heatmap(data,cluster_cols = T,topnum = 30)
#dev.off()
