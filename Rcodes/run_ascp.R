

# ascp root: vog.hin.mln.ibcn.ptf@ptfnona:
#Remainder of path:
#    /sra/sra-instant/reads/ByRun/sra/{SRR|ERR|DRR}/<first 6 characters of accession>/<accession>/<accession>.sra
#Where
#{SRR|ERR|DRR} should be either ‘SRR’, ‘ERR’, or ‘DRR’ and should match the prefix of the target .sra file
#ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -QT -l 100m
# anonftpftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR161/SRR1610960/SRR1610960.sra

run_ascp <- function(sraid = "SRR1610960", maxspeed="200m", outdir="."){
    
    id1 <- substr(sraid, start=1, stop=3)
    id2 <- substr(sraid, start=1, stop=6)
    #http://www.ncbi.nlm.nih.gov/books/NBK158899/#SRA_download.downloading_sra_data_using
    cmd <- paste0("ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -QT -l ", maxspeed,
                  " anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/",
                  id1, "/", id2, "/", sraid, ".sra ", outdir)
    return(cmd)
}



