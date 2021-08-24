import sys

project_name = "{{cookiecutter.project_name}}"
if " " in project_name:
    print("project_name is used as the name of the directory, so should not have any spaces")
    sys.exit(1)
