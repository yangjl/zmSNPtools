canonical$exonSize <- sapply(1:nrow(canonical), 
                             function(i){
                               ac <- subset(exon, txid %in% canonical$transcript_id[i])
                               print(i);
                               return(sum(ac$size));
                             })