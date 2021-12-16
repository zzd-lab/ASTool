#compare2view.r
library(ggplot2)
#library(gridExtra)
#library(cowplot)
#library(ggpubr)
#library(reshape2)

args<-commandArgs(T)

depthfile<-unlist(strsplit(args[1],split=","))
designfile<-unlist(strsplit(args[3],split=","))
psi_values<-unlist(strsplit(args[2],split=","))
sample_names<-unlist(strsplit(args[7],split=","))
file_index=0
out<-data.frame(V1=0,V2=0,V3=0,sample="",Treatment="")
for (d_file in depthfile){
	file_index=file_index+1
	sample_name<-paste("PSI =",psi_values[file_index],"\n(",sample_names[file_index],")")
	data_d_file<-read.table(d_file,header=F,sep="\t")
	data_d_file2<-data.frame(data_d_file,sample=sample_name,Treatment=designfile[file_index])
	out<-rbind(out,data_d_file2)
}
out <- out[-1,]

pos<-as.numeric(unlist(strsplit(args[5],split="_")))
min_pos=min(pos)
max_pos=max(pos)

data2<-subset(out,out[,2]>=min_pos)
data1<-subset(data2,data2[,2]<=max_pos)
data<-subset(data1,data1[,3]>0)
x_lab<-args[6]
x_lab<-gsub(","," ",x_lab)
options(warn=-1)
pdf(args[4])
#ggarrange(p2, p1, heights = c(2,3), ncol = 1, nrow = 2,align = "v")
ggplot(data) + geom_bar(aes(x = data[,2], y =data[,3],fill=Treatment),stat = "identity")+ theme_bw()+theme(panel.grid=element_blank(),title=element_text(size=9), legend.title=element_text(size=9),axis.text.y=element_text(color="black"),axis.text.x = element_text(size=9,angle = 60, hjust = 1, vjust = 1,color="black"),plot.margin = unit(c(1,5,9,1),"cm"),strip.background = element_rect(fill = 'white', colour = 'white'))+ylab("Depth")+xlab(x_lab)+scale_x_continuous(limits=c(min_pos,max_pos),breaks=pos)+ scale_fill_brewer(palette="Set2")+ facet_wrap(~sample,nrow=2)
dev.off()