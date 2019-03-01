import smtplib


username = 'virusgwhet@gmail.com' 
password = 'xXYk4S8byXvG88u'
carriers = {'tmobile': '@tmomail.net', 'att': '@txt.att.net'}
            
def SendText(content,numbers=numbers,emails=emails):
    try:
        mailserver = smtplib.SMTP('smtp.gmail.com', 587)
        mailserver.ehlo()
        mailserver.starttls()
        mailserver.ehlo()
        mailserver.login(username, password)
        for number,carrier in numbers:
            mailserver.sendmail(username, number + carriers[carrier], content)
            print("Sent {} to {}.".format(content,number))
        contentsub = 'Subject: {}\n\n{}'.format('GW alert', content)
        for em in emails:
            
            mailserver.sendmail(username, em, contentsub)
            print("Sent {} to {}.".format(content,em))
        mailserver.close()
    except:
        print('Error sending email!')

