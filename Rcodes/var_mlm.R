
var_mlm <- function(df = pheno, nullmodel="~ (1 | GCA1.all) + (1 | GCA2.all)", fullmodel="~ (1 | GCA1.all) + (1 | GCA2.all) + (1 | SCA.all)"){
    
    # Structures to hold results
    r2marg <- data.frame(trait=vector(), r2marg=vector(), r2cond=vector())
    sig <- data.frame(trait=vector(), chisq=vector(), df=vector(), p.value=vector())
    
    # Construct MLMs for each of the 61 time points
    nms <- names(moon)
    idx1 <- which(nms == "V1")
    idx2 <- which(nms == "V61")
    
    for (i in idx1:idx2) {
        trait <- names(df)[2]
        
        # Construct the null model
        null.formula <- as.formula(paste(trait, nullmodel, sep=" "))
        model.null <- lmer(null.formula, data = df, REML = FALSE)
        
        # Construct the model
        full.formula <- as.formula(paste(trait, fullmodel, sep=" "))
        model.full <- lmer(full.formula, data = df, REML = FALSE)
        
        # Retrieve the p-value
        inter1a <- anova(model.null, model.full)[2, 6:8]
        inter1b <- data.frame(trait=trait, chisq=inter1a[1, 1], df=inter1a[1, 2], p.value=inter1a[1, 3])
        sig <- rbind(sig, inter1b)
        
        # Calculate the marginal R-squared
        inter2a <- rsquared.glmm(model.full)
        inter2b <- data.frame(trait=trait, r2marg=inter2a[1, 4], r2cond=inter2a[1, 5])
        r2marg <- rbind(r2marg, inter2b)
    }
    
    # Calculate q-values
    sig$q.value <- p.adjust(sig[, "p.value"], method="fdr")
    
    time.points <- seq(0, 180, 3)
    
    res <- data.frame(time=time.points, sig=sig, r2=r2marg)
    return(res)
}

