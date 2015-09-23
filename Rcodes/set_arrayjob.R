# Jinliang Yang
# function to prepare slurm script

set_arrayjob <- function(shid="largedata/GenSel/CL_test.sh",
                         shcode="sh largedata/myscript.sh",
                         arrayjobs="1-700",
                         wd=NULL, jobid="myjob", email=NULL){
    
    #message(sprintf("###>>> cp from Introgression, tailored for pvpDiallel"))  
  
    ##### setup working directory
    if(is.null(wd)){
       wd <- getwd()
    }
    dir.create("slurm-log", showWarnings = FALSE)
    sbath <- paste0(wd, "/slurm-log/")
    sbatho <- paste0(sbath, "testout-%j.txt")
    sbathe <- paste0(sbath, "err-%j.txt")
   
    #### parameters pass to slurm script
    cat(paste("#!/bin/bash -l"),
        #-D sets your project directory.
        #-o sets where standard output (of your batch script) goes.
        #-e sets where standard error (of your batch script) goes.
        #-J sets the job name.
        paste("#SBATCH -D", wd, sep=" "),
        paste("#SBATCH -o", sbatho, sep=" "),
        paste("#SBATCH -e", sbathe, sep=" "),
        paste("#SBATCH -J", jobid, sep=" "),
        paste0("#SBATCH --array=", arrayjobs),
        paste0("#SBATCH --mail-user=", email),
        paste("#SBATCH --mail-type=END"),
        paste("#SBATCH --mail-type=FAIL #email if fails"),
            
        
        "set -e",
        "set -u",
        "",
        #"module load gmap/2014-05-15",
        file=shid, sep="\n", append=FALSE);
    
    #### attach some sh scripts
    cat(shcode, file=shid, sep="\n", append=TRUE)
    message(paste("###>>> In this path: cd ", wd, sep=""), "\n",
            paste("###>>> [ note: --ntasks=INT, number of cup ]"),"\n",
            paste("###>>> [ note: --mem=16000, 16G memory ]"),"\n",
            paste("###>>> RUN: sbatch -p bigmemh", shid),
            "")
    
}

