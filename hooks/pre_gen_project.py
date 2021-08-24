import sys

project_name = "{{cookiecutter.project_name}}"
if " " in project_name:
    print("project_name is used as the name of the directory, so should not have any spaces")
    sys.exit(1)

genome_build = "{{cookiecutter.genome_build.lower()}}"
import yaml
with open("../{{cookiecutter.project_name}}/config/config.yaml") as c:
    config = yaml.load(c, yaml.Loader)
if genome_build not in config["reference"]:
    print(f"genome_build ({genome_build}) not found in config.yaml.  "
          f"Valid options are {','.join(config['reference'].keys())}")
    sys.exit(1)
