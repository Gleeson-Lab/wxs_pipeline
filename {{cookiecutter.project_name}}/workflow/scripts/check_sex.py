"""
Confirm the declared sex of the sample with the sex inferred from the pipeline.
"""
import os
from common import email_user_with_error

SEX_CHECK_FAILURE_MESSAGE = "Failed sex check."
# checked so as to not email more than once
already_failed = os.path.isfile(snakemake.params.failure_file)
if snakemake.log:
    log_fh = open(str(snakemake.log), "w")
else:
    import sys

    log_fh = sys.stdout
try:
    declared_sex_by_sample = {}
    with open(snakemake.input["original_ped"]) as ped_fh:
        for line in ped_fh:
            _, individual_id, _, _, sex, _ = line.strip().split("\t")
            declared_sex_by_sample[individual_id] = sex
    sex_mismatches = []
    with open(snakemake.input["predicted_ped"]) as ped_fh:
        for line in ped_fh:
            _, individual_id, _, _, sex, _ = line.strip().split("\t")
            if declared_sex_by_sample[individual_id] != sex:
                sex_mismatches.append(
                    (individual_id, declared_sex_by_sample[individual_id], sex)
                )
    if sex_mismatches:
        messages = [
            "The following samples had a mismatch between the "
            "declared sex and predicted sex:"
        ]
        for individual_id, declared_sex, predicted_sex in sex_mismatches:
            messages.append(f"{individual_id} {declared_sex=}, {predicted_sex=}")
    else:
        messages = ["No sex mismatches detected."]
    message = "\n".join(messages) + "\n"
    log_fh.write(message)
finally:
    if snakemake.log:
        log_fh.close()
if sex_mismatches:
    with open(snakemake.params.failure_file, "w") as failure_fh:
        failure_fh.write(message)
    if not already_failed:
        email_user_with_error(
            snakemake.params.email,
            subject=SEX_CHECK_FAILURE_MESSAGE,
            message=message,
        )
    raise ValueError(SEX_CHECK_FAILURE_MESSAGE)
# success
if already_failed:
    os.remove(snakemake.params.failure_file)
