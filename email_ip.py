import smtplib
import private as pr


carriers = {'tmobile': '@tmomail.net', 'att': '@txt.att.net'}
def SendText(content,numbers=pr.numbers,emails=pr.emails):
    try:
        mailserver = smtplib.SMTP('smtp.gmail.com', 587)
        mailserver.ehlo()
        mailserver.starttls()
        mailserver.ehlo()
        mailserver.login(pr.username, pr.password)
        for number,carrier in numbers:
            mailserver.sendmail(pr.username, number + carriers[carrier], content)
            print("Sent {} to {}.".format(content,number))
        contentsub = 'Subject: {}\n\n{}'.format('GW alert', content)
        for em in emails:
            mailserver.sendmail(pr.username, em, contentsub)
            print("Sent {} to {}.".format(content,em))
        mailserver.close()
    except:
        print('Error sending email!')

def main():
    SendText('TEST',numbers=[],emails=['majoburo@gmail.com'])

if __name__=='__main__':
    main()
