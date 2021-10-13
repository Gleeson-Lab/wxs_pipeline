import sys

project_name = "{{cookiecutter.project_name}}"
if " " in project_name:
    print(
        "project_name is used as the name of the directory, so should not have any spaces"
    )
    sys.exit(1)

max_contamination = "{{cookiecutter.max_intraspecies_contamination}}"
try:
    max_contamination = float(max_contamination)
    if max_contamination < 0 or max_contamination > 1:
        raise ValueError
except ValueError:
    print(
        f"max_intraspecies_contamination ({max_contamination}) must be in the range [0,1]."
    )
    sys.exit(1)

min_relatedness = "{{cookiecutter.min_identity_relatedness}}"
try:
    min_relatedness = float(min_relatedness)
    if min_relatedness < 0 or min_relatedness > 1:
        raise ValueError
except ValueError:
    print(f"min_identity_relatedness ({min_relatedness}) must be in the range [0,1].")
    sys.exit(1)

ploidy = "{{cookiecutter.ploidy}}"
try:
    ploidy = int(ploidy)
    if ploidy < 1:
        raise ValueError
except ValueError:
    print(f"{ploidy=} must be an integer at least 1")
    sys.exit(1)
