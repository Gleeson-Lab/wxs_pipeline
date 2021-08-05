import os

user = os.getlogin()

logs = "{{cookiecutter.logs}}".replace("$USER", user)
output = "{{cookiecutter.output}}".replace("$USER", user)
resources = "{{cookiecutter.resources}}".replace("$USER", user)
scratch = "{{cookiecutter.resources}}".replace("$USER", user)

for link, directory in (
    ("logs", logs),
    ("output", output),
    ("resources", resources),
    ("scratch", scratch),
):
    if not os.path.isdir(directory):
        os.makedirs(directory)
        os.symlink(directory, f"../{link}")
