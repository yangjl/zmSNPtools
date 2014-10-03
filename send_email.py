#usr/bin/python
from optparse import OptionParser

def send_email(subject="Subject", txt="message"):
    import smtplib
    gmail_user = "yangjl0930@gmail.com"
    gmail_pwd = "yangjl@810903"
    FROM = 'yangjl0930@gmail.com'
    TO = ['yangjl0930@gmail.com'] #must be a list
    SUBJECT = subject
    TEXT = txt

    # Prepare actual message
    message = """\From: %s\nTo: %s\nSubject: %s\n\n%s
    """ % (FROM, ", ".join(TO), SUBJECT, TEXT)
    try:
        #server = smtplib.SMTP(SERVER)
        server = smtplib.SMTP("smtp.gmail.com", 587) #or port 465 doesn't seem to work!
        server.ehlo()
        server.starttls()
        server.login(gmail_user, gmail_pwd)
        server.sendmail(FROM, TO, message)
        #server.quit()
        server.close()
        print 'successfully sent the mail'
    except:
        print "failed to send mail"

def main():
    usage = "Usage: %prog -s subject -m message"
    parser = OptionParser(usage=usage)
    parser.add_option("-s", "--subject", dest="subject", help="subject of the email", default="Job Done")
    parser.add_option("-m", "--message", dest="message", help="message of the email", default="Job Done")

    (options, args) = parser.parse_args()
    send_email(subject=options.subject, txt=options.message)

if __name__ == '__main__':
    main()

