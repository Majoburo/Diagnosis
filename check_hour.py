from astropy.time import Time
import email_ip
import numpy

now = Time.now().jd
last = np.loadtxt('time_last.txt')
if now - last > 1/24.:
    email_ip.SendText("DIDN'T SEND GCN TEST ALERT WITHIN HOUR. CHECK DIAGNOSIS.",numbers=[('5125763501','tmobile')],emails=['majoburo@gmail.com'])
