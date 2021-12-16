library(ggplot2)
library(ggrepel)

options(warn=-1)
args<-commandArgs(T)
edges<-read.table(args[1],header=F,sep="\t")
colnames(edges)<-c("x","xend","count","Junction_type","pos","x-lab")
edges_in=subset(edges,edges[,4]=="Inclusion reads")
edges_ex=subset(edges,edges[,4]=="Exclusion reads")
edges_new=subset(edges,edges[,4]=="New junction reads")
pos<- sort(unique(as.numeric(unlist(strsplit(unlist(as.vector(paste(edges[,5],collapse=","))),split=',')))))
min_pos=min(pos)
max_pos=max(pos)
data<-read.table(args[2],header=F,sep="\t")
data2<-subset(data,data[,2]>=min_pos)
data1<-subset(data2,data2[,2]<=max_pos)
data<-subset(data1,data1[,3]>0)
y_max=max(data[,3])
pdf(args[3])
if(nrow(edges_new)>0){
edges[,4]<-factor(edges[,4],c("Inclusion reads","Exclusion reads","New junction reads"))
x_lab=edges[1,6]
ggplot(data)+geom_bar(aes(x = data[,2], y =data[,3]),stat = "identity",fill="grey")+geom_curve(aes(x = x, y = 0, xend = xend, yend = 0),data = edges_ex,color="#FF5733", curvature = -0.4,alpha = 1)+geom_curve(aes(x = x, y = 0, xend = xend, yend = 0),data = edges_new,color="#7CB342", curvature = -0.4,alpha = 1)+geom_text_repel(aes(x = (xend-x)/2+x, y = (y_max-0)/(max_pos-min_pos)*0.4*(xend-x), label = count,color=Junction_type),data = edges,size = 4,fontface = "bold")+theme_bw()+theme(panel.grid=element_blank(),plot.margin = unit(c(0,5,1,1),"cm"),aspect.ratio=0.5,axis.text.y=element_text(color="black"),axis.text.x = element_text(size=9,angle = 60, hjust = 1, vjust = 1,color="black"),title=element_text(size=9))+ylab("Deapth")+xlab(x_lab)+scale_x_continuous(limits=c(min_pos,max_pos),breaks=pos)+guides(color=guide_legend(title=NULL))+scale_y_continuous(limits=c(0,y_max))+scale_color_manual(values = c("#0D47A1","#641E16","#1B5E20"))
}else{
edges[,4]<-factor(edges[,4],c("Inclusion reads","Exclusion reads"))
x_lab=edges[1,6]
ggplot(data)+geom_bar(aes(x = data[,2], y =data[,3]),stat = "identity",fill="grey")+geom_curve(aes(x = x, y = 0, xend = xend, yend = 0),data = edges_ex,color="#FF5733", curvature = -0.4,alpha = 1)+geom_text_repel(aes(x = (xend-x)/2+x, y = (y_max-0)/(max_pos-min_pos)*0.4*(xend-x), label = count,color=Junction_type),data = edges,size = 4,fontface = "bold")+theme_bw()+theme(panel.grid=element_blank(),plot.margin = unit(c(0,5,1,1),"cm"),aspect.ratio=0.5,axis.text.x = element_text(size=9,angle = 60, hjust = 1, vjust = 1,color="black"),axis.text.y=element_text(color="black"),title=element_text(size=9))+ylab("Depth")+xlab(x_lab)+scale_x_continuous(limits=c(min_pos,max_pos),breaks=pos)+guides(color=guide_legend(title=NULL))+scale_y_continuous(limits=c(0,y_max))+scale_color_manual(values = c("#0D47A1","#641E16"))
}
dev.off()
