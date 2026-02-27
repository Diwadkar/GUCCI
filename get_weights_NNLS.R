library(rareGWAMA);
library(pracma)

eqtl <- as.data.frame(fread(file='/storage/group/dxl46/default/private/dajiang/projects/cci/merged.sc.bk.eqtl.genename.txt',header=T));

#Filter beta and SE are remove any NAs
betas <- eqtl[,8:22];
ses <- eqtl[,23:37];
ix.keep <- which(rowSums(is.na(betas))==0);
betas.nona <- betas[ix.keep,];
ses.nona <- ses[ix.keep,];

# Set constraints
Aeq <- matrix(rep(1, ncol(betas.nona)-1), nrow= 1)
beq <- c(1)

## Lower and upper bounds of the parameters, i.e [0, 1]
lb <- rep(0, ncol(betas.nona)-1)
ub <- rep(1, ncol(betas.nona)-1)

## Use function that solves linearly constrained linear least-squares problems
frac <- lsqlincon(matrix(unlist(betas.nona[,-1]),nrow=nrow(betas.nona),ncol=ncol(betas.nona)-1),
                  unlist(betas.nona[,1]),
                  Aeq= Aeq,
                  beq= beq,
                  lb= lb,
                  ub= ub)


## Generate dataframe with weights
frac <- t(as.data.frame(frac))
colnames(frac) <- gsub("BETA_","",colnames(betas)[-1])
weight <- frac
write.table(weight, "merged.sc.bk.eqtl.genename.full.weights.txt", sep="\t", quote=F, row.names=F)
