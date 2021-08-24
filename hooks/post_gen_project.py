import os
import yaml

genome_build = "{{cookiecutter.genome_build.lower()}}"
with open("config/config.yaml") as c:
    config = yaml.load(c, yaml.Loader)
if genome_build not in config["reference"]:
    print(f"genome_build ({genome_build}) not found in config.yaml.  "
          f"Valid options are {','.join(config['reference'].keys())}")
    sys.exit(1)

user = os.getlogin()
logs = "{{cookiecutter.logs}}".replace("$USER", user)
output = "{{cookiecutter.output}}".replace("$USER", user)
resources = "{{cookiecutter.resources}}".replace("$USER", user)
scratch = "{{cookiecutter.scratch}}".replace("$USER", user)

for link, directory in (
    ("logs", logs),
    ("output", output),
    ("resources", resources),
    ("scratch", scratch),
):
    if not os.path.isdir(directory):
        os.makedirs(directory)
    os.symlink(directory, link)
