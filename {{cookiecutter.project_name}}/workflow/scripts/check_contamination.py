"""
Check the FREEMIX field of the verifybamid2 output, which is the estimated
amount of contamination.
"""
CONTAMINATION_FIELD = "FREEMIX"
fns_with_excess_contamination = []
for fn in snakemake.input:
    with open(fn) as fh:
        header = fh.readline().strip().split("\t")
        data = fh.readline().strip().split("\t")
    d = dict(zip(header, data))
    contamination = float(d[CONTAMINATION_FIELD])
    print(f"{fn} has {contamination} contamination.")
    if contamination > snakemake.params.max_contamination:
        fns_with_excess_contamination.append(fn)
if fns_with_excess_contamination:
    prinf(f"the following: {fns_with_excess_contamination} had greater than "
          f"{snakemake.params.max_contamination}.")
