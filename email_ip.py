import smtplib
import email.utils
from email.mime.image import MIMEImage
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication
from email.mime.text import MIMEText
import os.path
import datetime
import importlib
carriers = {'gphi':'@msg.fi.google.com','verizon':'@vtext.com','tmobile': '@tmomail.net', 'att': '@txt.att.net'}
def SendText(content,emailcontent=None,plotfiles=[],lstfile=None,numbers=None,emails=None,recipients='recipients.py'):
    try:
        pr = importlib.import_module(recipients.split('.')[0])
        if numbers == None:
            numbers = pr.numbers
        if emails == None:
            emails = pr.emails
        mailserver = smtplib.SMTP('smtp.gmail.com', 587)
        mailserver.ehlo()
        mailserver.starttls()
        mailserver.ehlo()
        mailserver.login(pr.username, pr.password)
        msg = MIMEMultipart()
        msg['Subject'] = 'GW ALERT'
        msg['From'] = pr.username
        #msg['Date'] = email.utils.localtime()
        #msg.preamble = 'Gravitational Wave Alert'
        msg.attach(MIMEText(content))
        for number,carrier in numbers:
            msg['To'] = number+carriers[carrier]
            mailserver.sendmail(pr.username, number + carriers[carrier], msg.as_string())
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
    SendText('TEST  ',plotfiles=['LSTs_MS181101ab.pdf','MOLL_GWHET_M2052.pdf'],numbers=[('5125763501','tmobile')],emails=['majoburo@gmail.com'])

if __name__=='__main__':
    main()
