"""
Check the FREEMIX field of the verifybamid2 output, which is the estimated
amount of contamination.
"""
import os
from common import email_user_with_error

CONTAMINATION_FAILURE_MESSAGE = (
    f"Failed intra-species contamination check on {snakemake.wildcards.sample}."
)
# checked so as to not email more than once
already_failed = os.path.isfile(snakemake.params.failure_file)
if snakemake.log:
    log_fh = open(str(snakemake.log), "w")
else:
    import sys

    log_fh = sys.stdout
try:
    CONTAMINATION_FIELD = "FREEMIX"
    fns_with_excess_contamination = []
    log_lines = ["file_name\tcontamination\tflag"]
    for fn in sorted(snakemake.input):
        with open(fn) as fh:
            header = fh.readline().strip().split("\t")
            data = fh.readline().strip().split("\t")
        d = dict(zip(header, data))
        contamination = float(d[CONTAMINATION_FIELD])
        fn = os.path.basename(fn)
        if contamination > snakemake.params.max_contamination:
            fns_with_excess_contamination.append(fn)
            flag = "error"
        else:
            flag = "ok"
        log_lines.append(f"{fn}\t{contamination}\t{flag}")
    log_fh.write("\n".join(log_lines) + "\n")
finally:
    if snakemake.log:
        log_fh.close()
if fns_with_excess_contamination:
    message = (
        f"The following: {', '.join(fns_with_excess_contamination)} "
        f"had greater than {snakemake.params.max_contamination} intra-species "
        "contamination.\n"
    )
    with open(snakemake.params.failure_file, "w") as failure_fh:
        failure_fh.write(message)
    if not already_failed:
        email_user_with_error(
            snakemake.params.email, subject=CONTAMINATION_FAILURE_MESSAGE, message=message
        )
    raise ValueError(CONTAMINATION_FAILURE_MESSAGE)
# success
if already_failed:
    os.remove(snakemake.params.failure_file)
