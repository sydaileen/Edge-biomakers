# Compare the classification performance of node/edge features 
# selected by SFFS and SVM
# Usage: Rscript demo.R data labl
# data: m by n matrix of node data, where m is the number
#       of node features, n is the number of samples; 
#		rownames are feature names; colnames are sample names
# labl: n by 1 matrix labeling the sample classes: 1,2,3,...,K;
#		rownames are sample names,colnames as "grp";
# Created by Yidi Sun (sunyidi@picb.ac.cn)
# 22:30 Feb.28th, 2017

args=commandArgs(TRUE)
file1=args[1]
file2=args[2]

# read data(nodes) and group information
data=read.table(file1,header=TRUE,row.names=1)
labl=read.table(file2,header=FALSE,row.names=1)	

data=data[apply(data,1,function(x){ !(sum(is.na(x))/length(x)>0.5) }), ]
impute=function(data){
	for (i in 1:nrow(data)){
		tmp=data[i,]
		tmp[is.na(tmp)]=mean(tmp,na.rm=T)
		data[i,]=tmp
	}
	data
}
data=impute(as.matrix(data))

## Differential expressed node features

diff_node=function(data,labl,ffparam=0.05){	
	feature=rownames(data)
	sample=colnames(data)

	nf=dim(data)[1]
	ns=dim(data)[2]
	nc=length(unique(labl$grp))

	pval=matrix(0,dim(data)[1],nc*(nc-1)/2)
	l=0
	for (i in 1:(nc-1)){
		for (j in (i+1):nc){
			l=l+1
			lablindx1=which(labl$grp==i)
			datatmp1=data[,lablindx1]
			lablindx2=which(labl$grp==j)
			datatmp2=data[,lablindx2]
			datatmp=cbind(datatmp1,datatmp2)
			tmpn1=dim(datatmp1)[2]
			tmpn=dim(datatmp)[2]
			pval[,l]=apply(datatmp,1,function(x){t.test(x[1:tmpn1],x[(tmpn1+1):tmpn])$p.value})
		}
	}
	idx=c()
	for (i in 1:nrow(pval)){
		if(any(pval[i,]<ffparam)){
			idx=c(idx,i)
		}
	}
	data_node=data[idx,]
	data_node
}

# data_node=diff_node(data,labl)
snp_train_node=diff_node(snp_train,sam_train)
pro_train_node=diff_node(pro_train,sam_train)
car_train_node=diff_node(car_train,sam_train)
all_train_node=diff_node(all_train,sam_train)

## SVM cross-validation for classification

svmCV=function(data,labl,kfold=10){
	db=cbind(t(data),labl)
	colnames(db)[dim(db)[2]]="grp"
	sv=svm(grp~.,data=db,cross=kfold,type='C-classification',kernel='linear') 
	acc=sv$accuracies
	acc
}

# kfold=100
# acc_node=svmCV(data_node,labl,kfold=kfold)
# acc_edge=svmCV(data_edge,labl)

## SFFS algorithm for feature selection
library(e1071)

maxIter=100
kfold=10

## cross validation for node feature sets

result_node=sffs(data_node,labl,maxIter,kfold)
indx_node=result_node[[1]]
acc_node=result_node[[2]]

snp_train_node=sffs(snp_train,sam_train,maxIter,kfold)
pro_train_node=sffs(pro_train,sam_train,maxIter,kfold)
car_train_node=sffs(car_train,sam_train,maxIter,kfold)
lip_train_node=sffs(lip_train,sam_train,maxIter,kfold)
aa_train_node=sffs(aa_train,sam_train,maxIter,kfold)

all_train1_node=sffs(all_train1,sam_train,maxIter,kfold)
all_train2_node=sffs(all_train2,sam_train,maxIter,kfold)
all_train3_node=sffs(all_train3,sam_train,maxIter,kfold)
all_train4_node=sffs(all_train4,sam_train,maxIter,kfold)

all_train_node_snp=sffs(all_train_snp,sam_train,maxIter,kfold)
all_train_node_pro=sffs(all_train_pro,sam_train,maxIter,kfold)
all_train_node_car=sffs(all_train_car,sam_train,maxIter,kfold)
all_train_node_lip=sffs(all_train_lip,sam_train,maxIter,kfold)
all_train_node_aa=sffs(all_train_aa,sam_train,maxIter,kfold)

## transforming node data into edge features

# with deafault parameters
EdgeStruct=node2edge(data,labl)   
edata=EdgeStruct[[1]]

# with self-defined parameters
ffparam=0.05
pccparam=0.5
EdgeStruct=node2edge(data,labl,ffparam,pccparam)
edata=EdgeStruct[[1]]

snpEdge=node2edge(snp_train,sam_train,pccparam = 0.5, ffparam = 0.17)
snp_edge=snpEdge[[1]]
snp_edge1=getEdge(snpEdge,snp_test,sam_test)
proEdge=node2edge(pro_train,sam_train,pccparam = 0.5, ffparam = 0.15)
pro_edge=proEdge[[1]]
pro_edge1=getEdge(proEdge,pro_test,sam_test)
proTMTEdge=node2edge(proTMT_train,sam_train,pccparam = 0.5, ffparam = 0.2)
proTMT_edge=proTMTEdge[[1]]
proTMT_edge1=getEdge(proTMTEdge,proTMT_test,sam_test)
carEdge=node2edge(car_train,sam_train,pccparam = 0.5, ffparam = 0.3)
car_edge=carEdge[[1]]
car_edge1=getEdge(carEdge,car_test,sam_test)
lipEdge=node2edge(lip_train,sam_train,pccparam = 0.5, ffparam = 0.3)
lip_edge=lipEdge[[1]]
lip_edge1=getEdge(lipEdge,lip_test,sam_test)
lipXLEdge=node2edge(lipXL_train,sam_train,pccparam = 0.5, ffparam = 0.3)
lipXL_edge=lipXLEdge[[1]]
lipXL_edge1=getEdge(lipXLEdge,lipXL_test,sam_test)
aaEdge=node2edge(aa_train,sam_train,pccparam=0.5,ffparam=0.3)
aa_edge=aaEdge[[1]]
aa_edge1=getEdge(aaEdge,aa_test,sam_test)

all1Edge=node2edge(all_train1,sam_train,pccparam = 0.5, ffparam = 0.1)
all1_edge=all1Edge[[1]]
all1_edge1=getEdge(all1Edge,all_test1,sam_test)
all2Edge=node2edge(all_train2,sam_train,pccparam = 0.5, ffparam = 0.1)
all2_edge=all2Edge[[1]]
all2_edge1=getEdge(all2Edge,all_test2,sam_test)
all3Edge=node2edge(all_train3,sam_train,pccparam = 0.5, ffparam = 0.1)
all3_edge=all3Edge[[1]]
all3_edge1=getEdge(all3Edge,all_test3,sam_test)
all4Edge=node2edge(all_train4,sam_train,pccparam = 0.5, ffparam = 0.1)
all4_edge=all4Edge[[1]]
all4_edge1=getEdge(all4Edge,all_test4,sam_test)

allEdge=node2edge(all_train,sam_train,pccparam = 0.5, ffparam = 0.05)
all_edge=allEdge[[1]]
all_edge1=getEdge(allEdge,all_test,sam_test)

snp_train_edge=sffs(snp_edge,sam_train,maxIter,kfold)
pro_train_edge=sffs(pro_edge,sam_train,maxIter,kfold)
car_train_edge=sffs(car_edge,sam_train,maxIter,kfold)
lip_train_edge=sffs(lip_edge,sam_train,maxIter,kfold)
lipXL_train_edge=sffs(lipXL_edge,sam_train,maxIter,kfold)
aa_train_edge=sffs(aa_edge,sam_train,maxIter,kfold)

all_train1_edge=sffs(all1_edge,sam_train,maxIter,kfold)
all_train2_edge=sffs(all2_edge,sam_train,maxIter,kfold)
all_train3_edge=sffs(all3_edge,sam_train,maxIter,kfold)
all_train4_edge=sffs(all4_edge,sam_train,maxIter,kfold)
all_train_edge=sffs(all_edge,sam_train,maxIter,kfold)

all_train_edge_snp=sffs(all_edge_snp,sam_train,maxIter,kfold)
all_train_edge_pro=sffs(all_edge_pro,sam_train,maxIter,kfold)
all_train_edge_car=sffs(all_edge_car,sam_train,maxIter,kfold)
all_train_edge_lip=sffs(all_edge_lip,sam_train,maxIter,kfold)
all_train_edge_aa=sffs(all_edge_aa,sam_train,maxIter,kfold)

## selecting edge features by SFFS algorithm
## cross validation for selected edge feature sets

# result_edge=sffs(edata,labl,maxIter,kfold)
# indx_edge=result_edge[[1]]
# acc_edge=result_edge[[2]]

# ## cross validation for node+edge feature sets

# result_all=sffs(rbind(edata,data_node),labl,maxIter,kfold)
# indx_all=result_all[[1]]
# acc_all=result_all[[2]]

## Cross validation in independent validation set

# result:result of sffs feature selection
# testSet: independent validation dataset with the same rows as the training set
# testlabl: independent validation dataset labls indicating the sample groups

testSVM=function(result,testSet,testlabl,name,kfold=10){
	colclass<-rep("numeric",2)
	colnames<-c("Features","Accuracy")
	acc<-read.table(text="",colClasses = colclass,col.names=colnames)
	
	for (k in 1:length(result[[2]])){
		index=result[[1]][[k]]
		datatmp=testSet[index,]
		dbtmp=cbind(t(datatmp),testlabl)
		colnames(dbtmp)[dim(dbtmp)[2]]="grp"
		svtmp=svm(grp~.,data=dbtmp,cross=kfold,type='C-classification',kernel='sigmoid') 
		acc[k,1]=length(index)
		acc[k,2]=mean(svtmp$accuracies)
	}
	acc$class=rep(name,nrow(acc))
	acc
}

# test_acc=testSVM(result_node,test,sam,10)
snp_test_node=testSVM(snp_train_node,snp_test,sam_test,"snp")
pro_test_node=testSVM(pro_train_node,pro_test,sam_test,"pro")
car_test_node=testSVM(car_train_node,car_test,sam_test,"car")
lip_test_node=testSVM(lip_train_node,lip_test,sam_test,"lip")
aa_test_node=testSVM(aa_train_node,aa_test,sam_test,"aa")
all_test_node=testSVM(all_train_node,all_test,sam_test,"All")

snp_test_edge=testSVM(snp_train_edge,snp_edge1,sam_test,"snp")
pro_test_edge=testSVM(pro_train_edge,pro_edge1,sam_test,"pro")
car_test_edge=testSVM(car_train_edge,car_edge1,sam_test,"car")
lip_test_edge=testSVM(lip_train_edge,lip_edge1,sam_test,"lip")
aa_test_edge=testSVM(aa_train_edge,aa_edge1,sam_test,"aa")
all_test_edge=testSVM(all_train_edge,all_edge1,sam_test,"All")

## Format sffs result file

format_acc=function(result,name){
	colclass<-rep("numeric",2)
	colnames<-c("Features","Accuracy")
	acctable<-read.table(text="",colClasses = colclass,col.names=colnames)
	for (k in 1:length(result[[2]])){
		acctable[k,1]=length(result[[1]][[k]])
		acctable[k,2]=result[[2]][k]
	}
	acctable$class=rep(name,nrow(acctable))
	acctable
}

snp_train_node_acc=format_acc(snp_train_node,"snp")
pro_train_node_acc=format_acc(pro_train_node,"pro")
car_train_node_acc=format_acc(car_train_node,"car")
lip_train_node_acc=format_acc(lip_train_node,"lip")
aa_train_node_acc=format_acc(aa_train_node,"aa")
all_train_node_acc=format_acc(all_train_node,"All")
all_train1_node_acc=format_acc(all_train1_node,"All1")
all_train2_node_acc=format_acc(all_train2_node,"All2")

all_train_node_snp_acc=format_acc(all_train_node_snp,"noSNP")
all_train_node_pro_acc=format_acc(all_train_node_pro,"noPRO")
all_train_node_car_acc=format_acc(all_train_node_car,"noCAR")
all_train_node_lip_acc=format_acc(all_train_node_lip,"noLIP")
all_train_node_aa_acc=format_acc(all_train_node_aa,"noAA")


snp_train_edge_acc=format_acc(snp_train_edge,"snp")
pro_train_edge_acc=format_acc(pro_train_edge,"pro")
car_train_edge_acc=format_acc(car_train_edge,"car")
lip_train_edge_acc=format_acc(lip_train_edge,"lip")
aa_train_edge_acc=format_acc(aa_train_edge,"aa")
all_train_edge_acc=format_acc(all_train_edge,"All")
all_train1_edge_acc=format_acc(all_train1_edge,"All1")
all_train2_edge_acc=format_acc(all_train2_edge,"All2")

all_train_edge_snp_acc=format_acc(all_train_edge_snp,"noSNP")
all_train_edge_pro_acc=format_acc(all_train_edge_pro,"noPRO")
all_train_edge_car_acc=format_acc(all_train_edge_car,"noCAR")
all_train_edge_lip_acc=format_acc(all_train_edge_lip,"noLIP")
all_train_edge_aa_acc=format_acc(all_train_edge_aa,"noAA")


## ggplot2 line plot

ggplot(data=acc, aes(x=featue, y=accuracy, group=type, shape=type, colour=type)) + 
    geom_line() +     # Set linetype by sex
    geom_point(size=3) +         # Use larger points, fill with white
#    expand_limits(y=0) +                       # 设置x y轴的起止范围,这里是y�?0开�?
#    scale_colour_hue(name="Sex of payer",      # Set legend title
#                     l=30)  +                  # Use darker colors (lightness=30)
#    scale_shape_manual(name="Sex of payer",
#                       values=c(22,21)) +      # Use points with a fill color
#    scale_linetype_discrete(name="Sex of payer") +

    xlab("Number of Features") + ylab("Total bill") + # Set axis labels
    ggtitle("Average bill for 2 people") +     # Set title
    theme_bw() +                          # 设置主题
    theme(legend.position=c(.7, .4))           # 设置图例的位�?

lineplot=function(acc,setname){
	colnames(acc)[2]="accuracy"
	acc$class=factor(acc$class, levels=unique(acc$class))
	f=ggplot(data=acc, aes(x=Features, y=accuracy, group=class, 
		shape=class, colour=class)) + geom_line() + geom_point(size=1) + 
		xlab("Number of Features") + ylab("Accuracy") + 
		ggtitle(paste("Cross validation on ",setname,sep="")) + 
		theme_bw()+theme(plot.title = element_text(hjust = 0.45))+
		xlim(0,30)+ylim(20,100)
	f
}

f1=lineplot(train_nodes,"Training Set")
f11=lineplot(test_nodes,"validation Set")
f12=lineplot(test2_nodes,"validation Set II")
f2=lineplot(train_edges,"Training Set")
f3=lineplot(train_all,"Training Set")
f21=lineplot(test_edges,"validation Set")
f22=lineplot(test2_edges,"validation Set II")
f31=lineplot(test1_all,"validation Set I")
f32=lineplot(test2_all,"validation Set II")

##define multiplot function
    multiplot <- function(..., plotlist=NULL, cols) {
        require(grid)

        # Make a list from the ... arguments and plotlist
        plots <- c(list(...), plotlist)

        numPlots = length(plots)

        # Make the panel
        plotCols = cols                          # Number of columns of plots
        plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols

        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
        vplayout <- function(x, y)
            viewport(layout.pos.row = x, layout.pos.col = y)

        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            curRow = ceiling(i/plotCols)
            curCol = (i-1) %% plotCols + 1
            print(plots[[i]], vp = vplayout(curRow, curCol ))
        }

    }

## plot in one frame
myplot=multiplot(f1,f2,f3,f11,f21,f31,f12,f22,f32, cols=3)


### 3D PCA plot

library(rgl)
pcnode=princomp(t(tmpdb), cor=TRUE, scores=TRUE)
pcedge=princomp(t(tmpdbedge), cor=TRUE, scores=TRUE)
plot3d(pcnode$scores[,1:3],col=groupsam,xlab="",ylab="",zlab="")
rgl.postscript("nodefeature3d.pdf","pdf")
plot3d(pcedge$scores[,1:3],col=groupsam,x="",y="",z="")
rgl.postscript("edgefeature3d.pdf","pdf")

grpmean=function(db,grp){
	data=cbind(t(db),grp)
	mean1=apply(data[data[,"grp"]==1,-ncol(data)],2,mean)
	mean2=apply(data[data[,"grp"]==2,-ncol(data)],2,mean)
	mean3=apply(data[data[,"grp"]==3,-ncol(data)],2,mean)
	mean=data.frame(cbind(mean1,mean2,mean3))
	mean
}

proname=function(db,proid){
	for (i in 1:nrow(db)){
		tmp=strsplit(rownames(db)[i],"_")[[1]]
		idx1=which(proid[,1]==tmp[1])
		idx2=which(proid[,1]==tmp[2])
		if(length(idx1)>0){
			db[i,"Source"]=proid[idx1,2][1]
		}else{
			db[i,"Source"]=tmp[1]
		}
		if(length(idx2)>0){
			db[i,"Target"]=proid[idx2,2][1]
		}else{
			db[i,"Target"]=tmp[2]
		}
		db[i,"Group"]=tmp[3]
	}
	db
}

