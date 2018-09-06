# Transforming node data to edge data
# Usage: (edata,EdgeData) = node2edge(file,labl,pccparam,ffparam)
# file: m by n matrix of node data, where m is the number
#       of node features, n is the number of samples; 
#		rownames are feature names; colnames are sample names
# labl: n by 1 matrix labeling the sample classes: 1,2,3,...,K;
#		rownames are sample names;
# pccparam: double; criterion for selecting differentially 
#           correlated feature pairs.(default 0.5) 
#           For a certain feature pair, there will be K PCCs 
#           with respect to K classes. Among these PCCs, if 
#           there are two PCCs between which the absolute 
#           difference is larger than the pccparam, then 
#           select the feature pair
# ffparam: double; criterion for filtering differential 
#          edge features.(default 0.05) 
#          For a certain edge feature, there will be K(K-1)/2
#          pvalues (one-against-one ttest) among K classes of 
#          samples. If all of the pvalues are smaller than 
#          ffparam, then select the edge feature
# EdgeData: a list of five items
# 	edata: l by n matrix of edge data, where l is the number
#        of edge features, n is the number of samples;
#		 rownames are edge feature names taking the form:
#		 "<name of node1>_<name of node2>_<class i>", 
#		 where name of node1 is prior to name of node2 
#        in lexicographic order
#	mu: m by K matrix including the mean value of 
#		each gene in each class K
#	std: m by K matrix including the standard error
#		of each gene in each class K
#	efeature: a vector including all the feature names 
#	of the edge feature
#	etable: l by 3 cell array; column 1 and 2 are 
#                     strings of nodes for each edge; 
#                     column 3 are labels of class whose 
#                     means and standard deviations are 
#                     used to transform the data

# Created by Yidi Sun (sunyidi@picb.ac.cn)

node2edge=function(data,labl,pccparam=0.5,ffparam=0.05){

#	data=read.table(file,header=TRUE,row.names=1)
	nc=length(levels(labl$grp))
	ns=dim(data)[2]
	nf=dim(data)[1]

	# Generalized z-score transformation
	zdata=matrix(0, nc*nf,ns)
	zdata=as.data.frame(zdata)
	mu=matrix(0, nf,nc)
	std=matrix(0, nf,nc)
	for (i in 1:nc){
		lablindx=which(labl$grp==levels(labl$grp)[i])
		datatmp=data[,lablindx]
		mu[,i]=apply(datatmp,1,mean)
		std[,i]=apply(datatmp,1,sd)
		zdata[(1:nf)+(i-1)*nf,]=(data - (matrix(1,1,ns) %x% mu[,i])) / (matrix(1,1,ns) %x% std[,i])
	}

	# Calculate of PCCs of each feature pair
	pcc=array(0,dim=c(nf,nf,nc))
	for (i in 1:nc){
		lablindx=which(labl$grp==levels(labl$grp)[i])
		datatmp=data[,lablindx]
		pcc[,,i]=cor(t(datatmp))
	}
	pcc[is.na(pcc)]=0
	# Absolute difference of PCCs
	adpcc=array(0,dim=c(nf,nf,nc*(nc-1)/2))
	l=0
	for (i in 1:(nc-1)){
		for (j in (i+1):nc){
			l=l+1
			adpcc[,,l]=abs(pcc[,,i]-pcc[,,j])
		}
	}
	edgeIdx1=c()
	edgeIdx2=c()
	for (m in 1:(nf-1)){
		for (n in (m+1):nf){
			if(any(adpcc[m,n,]>pccparam)){
				edgeIdx1=c(edgeIdx1,m)
				edgeIdx2=c(edgeIdx2,n)
			}
		}
	}
	edgeIdx_tmp=cbind(edgeIdx1,edgeIdx2)
	ne=dim(edgeIdx_tmp)[1]
	edgeIdx=matrix(0,nc*ne,2)
	feature=rownames(data)
	edge_tmp=cbind(feature[edgeIdx1],feature[edgeIdx2])
	edge=matrix(, nc*ne, 3)
	for (i in 1:nc){
		edgeIdx[((1:ne)+(i-1)*ne),]=edgeIdx_tmp+(i-1)*nf
		edge[((1:ne)+(i-1)*ne),1:2]=edge_tmp
		edge[((1:ne)+(i-1)*ne),3]=i*rep(1,ne)
	}

	edata=zdata[edgeIdx[,1],]*zdata[edgeIdx[,2],]
	colnames(edata)=colnames(data)
	edata=na.omit(edata)

	# Filtering differential expressed edge features
	pval=matrix(0,dim(edata)[1],nc*(nc-1)/2)
	l=0
	for (i in 1:(nc-1)){
		for (j in (i+1):nc){
			l=l+1
			lablindx1=which(labl$grp==levels(labl$grp)[i])
			edatatmp1=edata[,lablindx1]
			lablindx2=which(labl$grp==levels(labl$grp)[j])
			edatatmp2=edata[,lablindx2]
			edatatmp=cbind(edatatmp1,edatatmp2)
			tmpn1=dim(edatatmp1)[2]
			tmpn=dim(edatatmp)[2]
			pval[,l]=apply(edatatmp,1,function(x){t.test(x[1:tmpn1],x[(tmpn1+1):tmpn])$p.value})
		}
	}
	idx=c()
	for (i in 1:nrow(pval)){
		if(all(pval[i,]<ffparam)){
			idx=c(idx,i)
		}
	}
	edata=edata[idx,]
	edge=edge[idx,]

	rownames(edata)=paste(edge[,1],edge[,2],edge[,3],sep="_")
	efeature=sort(unique(c(edge[,1],edge[,2])))
	EdgeData=list(edata,mu,std,efeature,edge)
}