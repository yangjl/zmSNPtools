# Jinliang Yang
# function to 

setUpslurm <- function(slurmsh="largedata/GenSel/CL_test.sh",
                       oneline=TRUE,
                       codesh="myscript.sh",
                       wd=NULL,
                       sbatho="/home/jolyang/Documents/pvpDiallel/slurm-log/testout-%j.txt",
                       sbathe="/home/jolyang/Documents/pvpDiallel/slurm-log/error-%j.txt",
                       sbathJ="jobid"){
    
    
    ##### setup working directory
    if(is.null(wd)){
       wd <- getwd()
    }
    
    #### parameters pass to slurm script
    cat(paste("#!/bin/bash"),
        #-D sets your project directory.
        #-o sets where standard output (of your batch script) goes.
        #-e sets where standard error (of your batch script) goes.
        #-J sets the job name.
        paste("#SBATCH -D", wd, sep=" "),
        paste("#SBATCH -o", sbatho, sep=" "),
        paste("#SBATCH -e", sbathe, sep=" "),
        paste("#SBATCH -J", sbathJ, sep=" "),
        "set -e",
        "set -u",
        "",
        file=slurmsh, sep="\n", append=FALSE);
    
    #### attach some sh scripts
    if(oneline){
       cat(codesh, file=slurmsh, sep="\n", append=TRUE) 
    }else{
        cat(paste("sh", codesh),
            file=slurmsh, sep="\n", append=TRUE)
    }
    
    #### warning and message
    cat("",
        paste("python /home/jolyang/bin/send_email.py -s", slurmsh),
        file=slurmsh, sep="\n", append=TRUE);
    
    message(paste("###>>> In this path: cd ", wd, sep=""), "\n",
            paste("###>>> note --ntask=x, 8GB of memory per CPU"),"\n",
            paste("###>>> RUN: sbatch -p bigmemh --mem 24000", slurmsh),
            "")
    
}

