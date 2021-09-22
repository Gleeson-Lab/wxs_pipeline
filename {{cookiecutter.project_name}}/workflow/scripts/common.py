import subprocess
from email.mime.text import MIMEText

def email_user_with_error(email, subject, message):
    print(f"sending email from {{cookiecutter.project_name}} to {email} "
          f"with {subject=} and {message=}")
    msg = MIMEText(message)
    msg["From"] = "{{cookiecutter.project_name}}"
    msg["To"] = email
    msg["Subject"] = subject
    p = subprocess.Popen(["/usr/sbin/sendmail", "-t", "-oi"],
                         stdin=subprocess.PIPE, universal_newlines=True)
    p.communicate(msg.as_string())
