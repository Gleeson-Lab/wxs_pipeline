# wxs_pipeline
This pipeline takes arbitrary numbers of W[EG]S samples in BAM, or FASTQ format to create single-sample BAMs, gVCFs, VCFs, and a joint-genotyped VCF at the end, using the [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-).  The pipeline is currently configured to only support hg38 or GRCh37.  Please contact <brcopeland@ucsd.edu> if you need help with something else.

# Requirements:
1. [Install conda](https://docs.conda.io/en/latest/miniconda.html).
2. [Install snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).
3. [Install cookiecutter](https://cookiecutter.readthedocs.io/en/latest/installation.html#install-cookiecutter).
4. (Recommended): [Initialize snakemake cluster profile](https://github.com/Gleeson-Lab/pbs-torque).
5. [Create an SSH key](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent#generating-a-new-ssh-key) and [configure github to use it](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account).
6. Copy this repository with `cookiecutter git+ssh://git@github.com/Gleeson-Lab/wxs_pipeline.git`
---

# Configuration:
Please see [deploying the pipeline](../../wiki/Deploying-the-Pipeline) for more information.
