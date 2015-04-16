
sendEmail <- function(subject="Job Done", message="Job DONE!"){
  sys.arg <- paste("python ~/bin/send_email.py -s ", "'", subject, "'",
                   " -m ", "'", message, "'", sep=" ")
  system(sys.arg)
}