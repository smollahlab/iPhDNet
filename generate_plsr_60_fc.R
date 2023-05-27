# ----------------------------------------------------------------------------
# Author: Shamim Mollah  
# Created: 10-12-2016
# Last updated: 1-12-2017
#
# Histone prediction using PLSR
#-----------------------------------------------------------------------------

setwd("/Users/smollah")
rm(list = ls())
histon_data<-read.table("60histon_24hr_final.txt", sep="\t", header=T, row.names=1)
phospho_data<-read.table("peptide_24hr_final.txt", sep="\t", header=T, row.names=1)

library(plsdof)

XNormZ<-scale(as.matrix(phospho_data[,]))
somePDFPath = "measured_predicted.pdf"
pdf(file=somePDFPath) 
par(mfrow=c(2,3))

pls.cv<-function (X, y, k = 10, groups = NULL, m = ncol(X), use.kernel = FALSE, 
                  compute.covariance = FALSE,method.cor="pearson") 
{
    n <- nrow(X)
    p <- ncol(X)
    if (is.null(groups) == FALSE) {
        f = as.factor(groups)
        k = length(levels(f))
        my.names = levels(f)
        all.folds <- split(1:n, f)
    }
    if (is.null(groups) == TRUE) {
        f <- rep(1:k, length = n)
        my.names <- 1:k
        all.folds <- split(sample(1:n), f)
    }
    
    ntrain = vector(length = k)
    for (i in 1:k) {
        ntrain[i] = n - length(all.folds[[i]])
    }
    ntrain.min = min(ntrain)
    m = min(m, ntrain.min - 1, p)
    cv.error.matrix = matrix(0, k, m + 1)
    rownames(cv.error.matrix) = my.names
    colnames(cv.error.matrix) = 0:m
    cor.error.matrix<-cv.error.matrix
    for (i in seq(k)) {
        omit <- all.folds[[i]]
        Xtrain = X[-omit, , drop = FALSE]
        ytrain = y[-omit]
        Xtest = X[omit, , drop = FALSE]
        ytest = y[omit]
        pls.object <- pls.model(Xtrain, ytrain, m = m, Xtest = Xtest, 
                                ytest = ytest, compute.DoF = FALSE, use.kernel = use.kernel,method.cor=method.cor)
        cv.error.matrix[i, ] <- pls.object$mse
        cor.error.matrix[i, ] <- pls.object$cor
    }
    cv.error = apply(cv.error.matrix, 2, mean)
    cor.error<-apply(cor.error.matrix,2,mean)
    # change by smollah 11-13-2016
    
    k=2
    lst=""
    iter=length(pls.object$RSS) -1
    for (i in 1:iter)  {
        if ((pls.object$RSS[i] - pls.object$RSS[k]) <= 0.05) {
            lst=append(lst,i)
        }
        k=k+1  
    }
    m.opt = as.numeric(lst[2]) 
    # plot RSS to choose the optimal component
    plot(pls.object$RSS,main=m.opt)
    #m.opt <- which.min(cv.error) - 1
    m.opt.cor<-which.max(cor.error) - 1
    if (compute.covariance == TRUE) {
        use.kernel = FALSE
    }
    pls.object <- pls.model(X, y, m = max(m.opt, m.opt.cor,1), use.kernel = use.kernel, 
                            compute.DoF = compute.covariance, compute.jacobian = compute.covariance)
    intercept <- pls.object$intercept[m.opt + 1]
    
    coefficients <- pls.object$coefficients[, m.opt + 1]
    covariance <- pls.object$covariance
    # edited by SMollah on 10-31-2016
    DoF <- pls.object$DoF
    intercept.cor <- pls.object$intercept[m.opt.cor + 1]
    coefficients.cor <- pls.object$coefficients[, m.opt.cor + 1]
    if (compute.covariance == TRUE) {
        #covariancve.cor<-covariance[m.opt.cor + 1, , ]
        covariance <- covariance[m.opt + 1, , ]
    }
    outlist = list(cv.error.matrix = cv.error.matrix, cor.error.matrix=cor.error.matrix,cv.error = cv.error, cor.error=cor.error,
                   m.opt = m.opt, m.opt.cor=m.opt.cor,covariance = covariance, DoF = DoF, intercept = intercept, intercept.cor=intercept.cor,
                   coefficients = coefficients,coefficients.cor=coefficients.cor)
    class(outlist) = "plsdof"
    return(outlist)
}

# m = maximal number of Partial Least Squares components. Default is m=min(ncol(X),nrow(X)-1)
comp=min(ncol(XNormZ),nrow(XNormZ)-1)

# total histone marks are 60, therefore, 60 different response values (y)
histone_num = 60      # 60 histones in GCP 
peptide_num = 96      # phosphoproteins in P100 
ii=1
cc=0
zz <- file("all.Rout", open="wt")
sink(zz, type="message")
for (i in 1:histone_num ) {
	set.seed(1234)
    yCent<-scale(as.vector(histon_data[,i]), scale = FALSE)
        # compute PLS coefficients for all the components (m) and plot Degrees of Freedom
        mypls1<-pls.model(XNormZ,yCent, m=comp, compute.DoF=TRUE)
        # add naive estimate
        #plot(0:comp,mypls1$DoF,pch="*",cex=3,xlab="principle components",ylab="DoF",ylim=c(0,max(mypls1$DoF)+1))
        # add naive estimate
        #lines(0:comp,0:comp,lwd=3)
        #mypls2<-pls.ic(XNormZ,yCent,criterion="bic")
        mypls3<-pls.cv(XNormZ,yCent,compute.covariance=TRUE,m=comp)
        my.vcov<-vcov(mypls3)
        my.sd<-sqrt(diag(my.vcov)) # standard deviation of the regression coefficients
        
        str=paste("p_list_1234",sep="")
        str2=paste("non_p_list_1234",sep="")
        
        index= mypls3$m.opt +1
        myvec = mypls3$coefficients
        mat=XNormZ%*%myvec
        
        plot(yCent,mat,xlab="measured",ylab="predicted(mycalc)", ylim=c(-2,2), xlim=c(-2,2),main=names(histon_data)[i])
        # add naive estimate
        lines(-2:2,-2:2,lwd=3)
        
        plot(yCent,mypls1$Yhat[,index],xlab="measured",ylab="predicted(Yhat)", ylim=c(-2,2), xlim=c(-2,2),main=names(histon_data)[i])
        
        # add naive estimate
        lines(-2:2,-2:2,lwd=3)
        
        
        for (k in 1:peptide_num ) {
            pval=dt(mypls3$coefficients[k]/my.sd[k], (mypls3$m.opt))
            if (pval < 0.01) {
                if (mypls3$coefficients[k] < 0){
                    direction = "neg"
                } 
                else {
                    direction = "pos" 
                } # end of 2nd if else
                l_str= paste(names(histon_data)[i], names(phospho_data)[k],(abs(mypls3$coefficients[k])),pval,direction)
                write.table(l_str, str, append=TRUE, sep="\t")
            } # end of 1st if
            else {
                l_str2= paste(names(histon_data)[i], names(phospho_data)[k],10*(abs(mypls3$coefficients[k])),pval)
                write.table(l_str2, str2, append=TRUE, sep="\t")
                
            }
        } 
    
    # end of 2nd for
} # end of 1st for
sink()
## Display the log file
readLines("all.Rout")
dev.off()
