x<-read.table("/Users/cjfiscus/Desktop/Working/SNP_rand_10k.hapmap") # opens file as data frame

xm=as.matrix(x[,c(1,12:length(x))]) # turns data frame into matrix, but skips columns 2-11

y<-t(as.matrix(xm)) # rotates the matrix

letters=c("AA","CC","GG","TT","NN") # define lettered genotypes (NN is missing)

numbers=c(11,22,33,44,99) # define new genotypes (99 for missing)

for( i in 1:5){ y[which(y==letters[i])]=numbers[[i]]} # for loop that replaces all AA with 11, CC with 22, etc.

write.table(y,"/Users/cjfiscus/Desktop/Working/reformatted.txt",quote=F,row.names=F,col.names=F,sep="\t") # writes new file to “reformatted.txt"
