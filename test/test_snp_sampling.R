# Jinliang Yang
# 3/26/2013
# Purpose: SNP formatting

obj <- load("~/Documents/XBSA/XBSA_proj/cache/xsnp0.RData")
#Note: decimal number generated during the indel merging
for(i in 8:19){
  xsnp0[,i] <- round(xsnp0[,i], 0)
}

### test version
test <- xsnp0[1:1000, c(1,2,3,8,9)]
test$refcum <- cumsum(test$XBSAlib1_REF_COUNT)
test$altcum <- sum(test$XBSAlib1_REF_COUNT) + cumsum(test$XBSAlib1_ALT_COUNT)
write.table(test, "~/Documents/Codes/zmSNPtools/test/test_sampling.txt", sep="\t",
            row.names=FALSE, quote=FALSE)


head(test)
res <- read.table("~/Documents/Codes/zmSNPtools/test/test.txt", header=T)
test <- merge(test, res, by="snpid")

test$ref <- test$XBSAlib1_REF_COUNT - test$ref_count
test$alt <- test$XBSAlib1_ALT_COUNT - test$alt_count
hist(test$ref, breaks=30)
hist(test$alt, breaks=30)

sum(test$ref)
sum(test$alt)
sum(test$XBSAlib1_REF_COUNT)
sum(test$XBSAlib1_ALT_COUNT)