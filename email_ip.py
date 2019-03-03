import smtplib
import private as pr
from email.mime.image import MIMEImage
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication
from email.mime.text import MIMEText
import os.path

carriers = {'tmobile': '@tmomail.net', 'att': '@txt.att.net'}
def SendText(content,emailcontent=None,plotfiles=[],lstfile=None,numbers=pr.numbers,emails=pr.emails):
    try:
        mailserver = smtplib.SMTP('smtp.gmail.com', 587)
        mailserver.ehlo()
        mailserver.starttls()
        mailserver.ehlo()
        mailserver.login(pr.username, pr.password)
        for number,carrier in numbers:
            mailserver.sendmail(pr.username, number + carriers[carrier], content)
            print("Sent {} to {}.".format(content,number))
        msg = MIMEMultipart()
        msg['Subject'] = 'GW alert'
        msg['From'] = pr.username
        msg.preamble = 'Gravitational Wave Alert'
        if emailcontent == None:
            msg.attach(MIMEText(content))
        else:
            msg.attach(MIMEText(emailcontent))
        if lstfile != None:
            with open(lstfile) as f:
                attach = MIMEApplication(f.read())
                attach.add_header('Content-Disposition', 'attachment', filename = lstfile)
                msg.attach(attach)
        for plotfile in plotfiles:
            if not os.path.exists(plotfile):
                print("Couldn't find " + plotfile)
                continue
            with open(plotfile,'rb') as f:
                attach = MIMEApplication(f.read(),'pdf')
                attach.add_header('Content-Disposition', 'attachment', filename = plotfile)
                msg.attach(attach)
        for em in emails:
            mailserver.sendmail(pr.username,em,msg.as_string())
            print("Sent alert to {}.".format(em))
        mailserver.close()
    except:
        print('Error sending email!')

def main():
    SendText('TEST',plotfiles=['LSTs_MS181101ab.pdf','MOLL_GWHET_M2052.pdf'],numbers=[],emails=['majoburo@gmail.com'])

if __name__=='__main__':
    main()
