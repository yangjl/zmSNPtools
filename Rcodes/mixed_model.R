#############################
library("nlme")
### an R function use nlme package to fit mixed model.
mixed_model <- function(data = df, model = KRN ~ Pedigree, random = ~1 | Farm/Rep, 
                        trait = "KRN") {
  
  trait <- as.character(model)[2]
  data <- data[!is.na(data[, trait]), ]
  data[, trait] <- as.numeric(as.character(data[, trait]))
  
  lmeout1 <- lme(model, data = data, random = random)
  ped.hat1 <- lmeout1$coef$fixed
  ped.hat1[-1] <- ped.hat1[-1] + ped.hat1[1]
  
  fix <- as.character(model)[3]
  names(ped.hat1)[1] <- as.character(data[order(data[, fix]), fix][1])
  names(ped.hat1) <- gsub(fix, "", names(ped.hat1))
  tped <- data.frame(Genotype = names(ped.hat1), trait = ped.hat1)
  names(tped)[2] <- trait
  return(tped)
}

