# Transforming node data to edge data
# Usage: edge_data = getEdge(EdgeData,file,labl)
# file: m by n matrix of node data, where m is the number
#       of node features, n is the number of test samples; 
#		rownames are feature names; colnames are sample names
# labl: n by 1 matrix labeling the sample classes: 1,2,3,...,K;
#		rownames are sample names;
# EdgeData: The structured edge list from function "node2edge"
# edge_data: l by n matrix of edge data, where l is the number
#        of edge features from EdgeData

# Created by Yidi Sun (sunyidi@picb.ac.cn)

getEdge=function(EdgeData,data,labl){

	mu=EdgeData[[2]]
	std=EdgeData[[3]]
	efeature=EdgeData[[4]]
	etable=EdgeData[[5]]
	idx=c()
	for (i in 1:length(efeature)){
		idx=c(idx,which(rownames(data)==efeature[i]))
	}
	data=data[idx,]
	mu=mu[idx,]
	std=std[idx,]
	nc=length(levels(labl$grp))
	ns=dim(data)[2]
	nf=dim(data)[1]

	# Generalized z-score transformation
	zdata=matrix(0, nc*nf,ns)
	zdata=as.data.frame(zdata)
	for (i in 1:nc){
		lablindx=which(labl$grp==levels(labl$grp)[i])
		datatmp=data[,lablindx]
		zdata[(1:nf)+(i-1)*nf,]=(data - (matrix(1,1,ns) %x% mu[,i])) / (matrix(1,1,ns) %x% std[,i])
	}
	edgeIdx1=c()
	edgeIdx2=c()
	for (i in 1:nrow(etable)){
		idx1=which(rownames(data)==etable[i,1])
		idx2=which(rownames(data)==etable[i,2])
		idx_nc=as.numeric(etable[i,3])
		idx1=idx1+(idx_nc-1)*nf
		idx2=idx2+(idx_nc-1)*nf
		edgeIdx1=c(edgeIdx1,idx1)
		edgeIdx2=c(edgeIdx2,idx2)
	}
	edata=zdata[edgeIdx1,]*zdata[edgeIdx2,]
	rownames(edata)=paste(etable[,1],etable[,2],etable[,3],sep="_")
	edata
}