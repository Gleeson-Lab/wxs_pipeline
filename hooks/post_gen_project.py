import os
import sys

# adjusts max gaussians to pass for VQSR
# https://gatk.broadinstitute.org/hc/en-us/articles/4402736812443-Which-training-sets-arguments-should-I-use-for-running-VQSR-
sequencing_type = "{{cookiecutter.sequencing_type}}"
if sequencing_type == "genome":
    snp_gaussians = "*max_gaussians_snp_genome"
    indel_gaussians = "*max_gaussians_indel_genome"
elif sequencing_type in ("exome", "amplicon"):
    snp_gaussians = "*max_gaussians_snp_other"
    indel_gaussians = "*max_gaussians_indel_other"
else:
    print(f"Unrecognized sequencing_type ({sequencing_type}.")
    sys.exit(1)
with open("config/config.yaml") as config_fh:
    config_content = config_fh.read()
config_content = config_content.replace("$((SNP_GAUSSIANS))", snp_gaussians)
config_content = config_content.replace("$((INDEL_GAUSSIANS))", indel_gaussians)
with open("config/config.yaml", "w") as config_fh:
    config_fh.write(config_content)

# create symlinks to locations data is stored in scratch space
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
    if os.path.isfile(link) or os.path.islink(link):
        os.remove(link)
    elif os.path.isdir(link):
        os.rmdir(link)
    if not os.path.isdir(directory):
        os.makedirs(directory)
    os.symlink(directory, link)
