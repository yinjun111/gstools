

library("argparser",quietly =T)

version="0.1"

description=paste0("gs_heatmap\nversion ",version,"\n","Usage:\nDescription: Generate heatmap image for p/q matrix\n")


#####
#Add arguments
#####

parser <- arg_parser(description=description)

parser <- add_argument(parser, arg="--in",short="-i", type="character", help = "Input file, a p/q value matrix")
parser <- add_argument(parser, arg="--out",short="-o", type="character", help = "Output file, png format")
parser <- add_argument(parser, arg="--sortby",short="-y", type="character", help = "Sort by which column")
parser <- add_argument(parser, arg="--top",short="-t", type="character", help = "# of top terms to show")

args = parse_args(parser)

print(args)



######
#Drawing function
#####

library(Cairo)
library(pheatmap)
library(RColorBrewer)

cols <- brewer.pal(9, "Reds")
pal <- colorRampPalette(cols)


gs_heatmap<-function(val,topnum=30,na.val=1,zero.val=0.0001,order.by="decreasing",sortby="avg",conversion="-log10",cluster_cols =F) {
  

  val[is.na(val)]<- na.val
  val[val==0]<-zero.val
  
  if(conversion == "-log10") {
    val<--log10(val)
  }
  
  topnum=min(topnum,nrow(val))
  
  
  if(order.by=="decreasing") {
	if(sortby=="avg") {
		val.top<-val[order(apply(val,1,sum),decreasing = T)[1:topnum],]
	} else {
		#by comparison name
		val.top<-val[order(val[,sortby],decreasing = T)[1:topnum],]
	}
  } else {
	if(sortby=="avg") {
		val.top<-val[order(apply(val,1,sum),decreasing = F)[1:topnum],]
	} else {
		#by comparison name
		val.top<-val[order(val[,sortby],decreasing = F)[1:topnum],]	
	}
  }
  
  #cell height and width will ensure pheatmap has border
  pheatmap(val.top,color=c("white",pal(10)),cluster_cols=cluster_cols,fontsize = 8,border_color = "grey40",cellheight=10,cellwidth=30)
}


######
#Input/output
#####

data<-read.table(args$"in",header=T,row.names=1,sep="\t",quote="",comment.char="")

#remove non-word chars
colnames(data)==gsub("\\W","_",colnames(data))
args$sortby=gsub("\\W","_",args$sortby)

#top terms 
topnum=as.numeric(args$top)
topnum=min(topnum,nrow(data))

figure1<-args$out

figure2=sub(".png$","_clustered.png",figure1,perl=T)


CairoPNG(file=figure1,res = 300,width = 15,height = 10,units = "in")
gs_heatmap(data,cluster_cols = F,topnum = topnum,sortby=args$sortby)
dev.off()

CairoPNG(file=figure2,res = 300,width = 15,height = 10,units = "in")
gs_heatmap(data,cluster_cols = T,topnum = topnum,sortby=args$sortby)
dev.off()
