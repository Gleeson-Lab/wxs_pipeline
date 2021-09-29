"""
Check the pairwise relatedness of each read group for the sample as reported by
somalier.
"""

import os
from common import email_user_with_error

RELATEDNESS_FAILURE_MESSAGE = (
    f"Failed read group relatedness check on {snakemake.wildcards.sample}."
)
# checked so as to not email more than once
already_failed =  os.path.isfile(snakemake.params.failure_file)
if snakemake.log:
    log_fh = open(str(snakemake.log), "w")
else:
    import sys

    log_fh = sys.stdout
try:
    rg_pairs_with_relatedness_problem = []
    log_lines = ["read_group_a\tread_group_b\trelatedness\tflag"]
    comparisons = {}
    with open(snakemake.input[0]) as somalier_fh:
        header = somalier_fh.readline().strip().split("\t")
        for line in somalier_fh:
            d = dict(zip(header, line.strip().split("\t")))
            # sample name format is rg_sample; all we want is the read group
            a = int(d["#sample_a"].split("_")[0])
            b = int(d["sample_b"].split("_")[0])
            if b < a:
                a, b = b, a
            relatedness = float(d["relatedness"])
            comparisons[(a, b)] = relatedness
    for read_group_a, read_group_b in sorted(comparisons):
        relatedness = comparisons[(read_group_a, read_group_b)]
        if relatedness < snakemake.params.min_relatedness:
            rg_pairs_with_relatedness_problem.append(
                f"{read_group_a},{read_group_b}"
            )
            flag = "error"
        else:
            flag = "ok"
        log_lines.append(f"{read_group_a}\t{read_group_b}\t{relatedness}\t{flag}")
    log_fh.write("\n".join(log_lines) + "\n")
finally:
    if snakemake.log:
        log_fh.close()

if rg_pairs_with_relatedness_problem:
    message = (
        f"The following read group pairs of {snakemake.wildcards.sample}: "
        f"{'; '.join(rg_pairs_with_relatedness_problem)} had less than "
        f"{snakemake.params.min_relatedness} relatedness.\n"
    )
    with open(snakemake.params.failure_file, "w") as failure_fh:
        failure_fh.write(message)
    if not already_failed:
        email_user_with_error(
            snakemake.params.email, subject=RELATEDNESS_FAILURE_MESSAGE, message=message
        )
    raise ValueError(RELATEDNESS_FAILURE_MESSAGE)
# success
if already_failed:
    os.remove(snakemake.params.failure_file)
