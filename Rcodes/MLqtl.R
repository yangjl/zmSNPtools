# Jinliang Yang
# 8.16.2012
# purpose: QTL functions

library(glmnet)

mrk2RNA <- function(idx=1, binsize=10000000, fgsinfo=canonical,
                    rc=fgsrc, mrk=seg_marker, alpha=1, dif=1){
  #################################################################
  # idx=1: the gene to fit
  # binsize=10000000: determined the cis region
  # fgsinfo=canonical: FGS annotation infomation
  # rc=fgsrc: reads count data
  # mrk=seg_marker: segment marker genotype data
  # alpha=1: determine which method to use, 1 for L1 and 0 for L2
  # dif=0: determine the different regulation over cis and trans region 
  #  	 dif=0(cis=1, trans=2); dif=-1(cis=1, trans=3), dif=1(cis=1, trans=1)			
  #################################################################
  
  geneid <- names(rc)[idx];
  geneinfo <- subset(fgsinfo, transcript_id==geneid)
  
  pos <- geneinfo$Mid;
  mrk$penalty <- dif;
  mrk[mrk$Chr==geneinfo$chromosome & mrk$Mid >=pos-binsize & mrk$Mid <= pos+binsize,]$penalty <- 1;
  
  tmrk <- as.data.frame(t(mrk[,-1:-7]));
  names(tmrk) <- mrk$Segmarker
  
  rc$ID <- row.names(rc);
  fullset <- merge(as.data.frame(rc[, c("ID", geneid)]), tmrk, by.x="ID", by="row.names")
  
  set.seed(1234*idx);
  idx2 <- sample(1:nrow(fullset), floor(nrow(fullset)*0.8))
  trainset <- fullset[idx2,]
  testset <- fullset[-idx2,]
  
  ####################### REGRESSION WITH REGULARIZATION
  
  penalty <- mrk$penalty;
  penalty <- 2-penalty;
  
  #cv.fit <- cv.glmnet(as.matrix(trainset[,-1:-2]), as.vector(trainset[, trait]), family="gaussian",
  #		nlambda=100, alpha=alpha, maxit=100000)
  
  fit <- glmnet(as.matrix(trainset[,-1:-2]), as.vector(trainset[, 2]), family="gaussian",
                nlambda=100, alpha=alpha, penalty.factor=penalty)
  
  
  
  pred <- predict(fit, newx=as.matrix(testset[,-1:-2]))
  pred <- as.data.frame(pred)
  pred$obs <- as.vector(testset[,2])
  corl <- cor(x=pred, as.vector(testset[,2]));
  corl <- corl[c(-1, -nrow(corl)),]
  minbeta <- names(corl)[which.max(abs(corl))];
  
  coef <- coef(fit)
  coef <- coef[, minbeta]
  coef <- coef[coef!=0]
  
  return(list(corl[which.max(abs(corl))], coef));	
  
}

################
RNA2pheno <- function(fgsinfo=canonical, dif=1, perm=1000, alpha=0,
                      rc=fgsrc, mrk=seg_marker, pheno=trait15, trait="SDW", penalty=res){
  
  fullset <- merge(pheno[, c("ID",trait)], rc, by.x="ID", by.y="row.names");
  fullset <- fullset[!is.na(fullset[,2]),]
  
  
  set.seed(1234);
  idx2 <- sample(1:nrow(fullset), floor(nrow(fullset)*0.8))
  trainset <- fullset[idx2,]
  testset <- fullset[-idx2,]
  
  ####################### REGRESSION WITH REGULARIZATION
  #if(abs(penalty)){
  penalty <- round((dif- abs(penalty)),0);
  
  #cv.fit <- cv.glmnet(as.matrix(trainset[,-1:-2]), as.vector(trainset[, trait]), family="gaussian",
  #  	nlambda=100, alpha=0, maxit=100000)
  
  fit <- glmnet(as.matrix(trainset[,-1:-2]), as.vector(trainset[, trait]), family="gaussian",
                nlambda=100, alpha=alpha, penalty.factor=penalty)
  
  pred <- predict(fit, newx=as.matrix(testset[,-1:-2]))
  pred <- as.data.frame(pred)
  pred$obs <- as.vector(testset[,2])
  corl <- cor(x=pred, as.vector(testset[,2]));
  corl <- corl[c(-1, -nrow(corl)),]
  minbeta <- names(corl)[which.max(abs(corl))];
  
  coef <- coef(fit)
  coef <- coef[, minbeta]
  coef <- coef[coef!=0]
  
  #####
  h2 <- lm(as.formula(paste(trait, "~", paste(names(coef)[-1], collapse = "+"), sep = "")),
           data=fullset)
  
  
  return(list(corl[which.max(abs(corl))], coef, h2));
  
}
#############################################################################