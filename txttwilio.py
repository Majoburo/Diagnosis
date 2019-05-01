from twilio.rest import Client

def sentxt(account_sid, auth_token, phone, message):
# Your Account Sid and Auth Token from twilio.com/console
# DANGER! This is insecure. See http://twil.io/secure
    client = Client(account_sid, auth_token)

    message = client.messages \
                    .create(
                         body=message,
                         from_='+17372042250',
                         to='+1'+phone
                     )

    #print(message.sid)
    return
