shared_genes<-intersect(rownames(bulk_dataset),rownames(Seurat_tmp))
sc_exprs <- as.matrix(Seurat_tmp@assays$RNA@data)
correlation_matrix<-cor(bulk_dataset[shared_genes,],sc_exprs[shared_genes,])
data = list(x = correlation_matrix, y = phenotype)
fit = SGL(data, leidengroup, type = "logit",alpha =0.5)
lam<-fit[['lambdas']]
cvfit<-cvSGL(data,leidengroup,type='logit',nfold = 5,alpha =0.5,lambdas = lam)
error<-cvfit$lldiff
m<-min(error)
h<-which(error==m)
a<-fit[["beta"]]
b<-a[,h]
LP_SGL_pos<-colnames(sc_dataset)[which(b>0)]
LP_SGL_neg<-colnames(sc_dataset)[which(b<0)]
Background<-colnames(sc_dataset)[which(b==0)]

