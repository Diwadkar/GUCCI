## Function
refine.cci <- function(bt,se,weight,omega) {
  ## b = A*beta + e;

  # Set beta and SE
  b <- as.vector(bt)
  beta.out <- rep(NA,length(b)-1)
  se <- as.vector(se)
  se.out <- beta.out

  #Indices for NAs and non-NAs
  ix.rm <- which(is.na(b))
  ix.keep <- which(!is.na(b))

  if(is.na(b[1])) {
    ##we cannot refine if the bk-eQTL is NA;
    return(list(beta=b[-1],
                se=se[-1],
                beta.diff=b[(length(b)-1):length(b)],
                se.diff=sqrt(sum((se[c(length(b)-1,length(b))])^2))))

  }

  #If there is data left
  if(length(ix.keep)>0){

    #If NAs present, update stats by removing NAs
    if(length(ix.rm)>0) {

      b <- unlist(b[-ix.rm])
      se <- unlist(se[-ix.rm])
      ## updating the weights
      weight <- unlist(weight[-ix.rm]/sum(weight[-ix.rm]))
      omega <- matrix(omega[-ix.rm,-ix.rm],nrow=nrow(omega)-length(ix.rm))
    }
    #print(length(se))
    #print(dim(omega))

    ##Parameter estimation: WLSQ
    omega <- ginv(diag(se)%*%omega%*%diag(se))

    A <- diag(length(b)-1)
    A <- matrix(unlist(rbind(A,weight)),ncol=length(b)-1)

    #Generate beta and covariance matrix
    beta <- ginv(t(A)%*%omega%*%A)%*%t(A)%*%omega%*%unlist(b)
    cov.beta <- ginv(t(A)%*%omega%*%A)

    ix.out <- ix.keep[which(!ix.keep%in%c(1))]-1

    beta.out[ix.out] <- beta
    se.out[ix.out] <- sqrt(diag(cov.beta))

    #Add covariances from covariance matrix
    cov.out <- matrix(NA,nrow=length(beta.out),ncol=length(beta.out))
    cov.out[ix.out,ix.out] <- cov.beta

    ##Effect size differences : Bg-Bl ; sqrt(seg^2 + sel^2 - 2 cov(B1, B2))
    beta.diff <- beta.out[length(beta.out)]-beta.out[length(beta.out)-1]
    se.diff <- sqrt(cov.out[length(beta.out),length(beta.out)]-2*cov.out[length(beta.out),length(beta.out)-1]+cov.out[length(beta.out)-1,length(beta.out)-1])

  } else {

    #In case all entries are NAs
    b <- as.vector(bt)
    ix.out <- ix.rm[which(!ix.rm%in%c(1))]-1
    beta.out[ix.out] <- rep(NA,length(b)-1)
    se.out[ix.out] <- rep(NA,length(b)-1)
    beta.diff <- NA
    se.diff <- NA

  }

  return(list(beta=beta.out,
              se=se.out,
              beta.diff=beta.diff,
              se.diff=se.diff))

}