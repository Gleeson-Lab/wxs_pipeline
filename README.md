# wxs_pipeline
This pipeline takes arbitrary numbers of W[EG]S samples in BAM and/or FASTQ format to create single-sample BAMs, gVCFs, and a joint-genotyped VCF at the end, using the [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-).  The pipeline is currently configured to only support hg38 or GRCh37.  Please contact <brcopeland@ucsd.edu> if you need help with something else.

# Requirements:
1. [Install conda](https://docs.conda.io/en/latest/miniconda.html).
2. [Install snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).
3. [Install cookiecutter](https://cookiecutter.readthedocs.io/en/latest/installation.html#install-cookiecutter).
4. (Recommended): [Initialize snakemake cluster profile](https://github.com/brcopeland/pbs-torque).
5. Copy this repository with `cookiecutter gh:brcopeland/wxs_pipeline` or `cookiecutter git+ssh://git@github.com/brcopeland/wxs_pipeline.git`
---

# Configuration:
Now within the main directory, the pipeline requires the directories `logs`, `output`, and `scratch` to exist.  These will be created elsewhere and a symlink generated as these should not be in your home directory or `/projects` for I/O and storage reasons.  A default value is provided which you can adjust (`$USER` will be replaced by your user name).

The pipeline will look to process every sample that it finds in `input/bams` and `input/fastqs`.  These directories already exist within the repository but each BAM/FASTQ should be a symlink (again for storage reasons) to where the file is actually located.  The pipeline expects each BAM to be named `{sample}.bam` in the `input/bams` directory, where `{sample}` matches the regular expression `[a-zA-Z\d_-]+`.  Likewise, the pipeline expects each FASTQ to be named `{sample}_{rg}.R{read}.fastq.gz` in the `input/fastqs` directory, where `{rg}` is a non-negative integer, and `{read}` is `1` or `2` (i.e., read one or read two).  Individual samples must have either one BAM or multiple FASTQs; some samples may have a BAM while others have FASTQs (but do not provide a BAM + FASTQs for the sample sample).

If a BAM is given for a sample, the read groups are discovered from the `@RG` tags in the header.  For FASTQs, each `{rg}` referred to above should be unique for the sample to be handled properly.

Here is an example of the commands necessary to set up the input links for sample `9988423398`: this sample has two read groups.
1. `ln -s /projects/ps-gleesonlab3/RCIGM_WGS_NTD_samples/202105_NTD_WGS_del22q11_subjects/210419_A01342_0041_AHYHG7DRXX/R21AB247-L1/UCSD-WGS-CSFFC2-IN74E9-R21AB247-L1-WGS-0-V1_S1_L001_R1_001.fastq.gz input/fastqs/9988423398_0.R1.fastq.gz`
2. `ln -s /projects/ps-gleesonlab3/RCIGM_WGS_NTD_samples/202105_NTD_WGS_del22q11_subjects/210419_A01342_0041_AHYHG7DRXX/R21AB247-L1/UCSD-WGS-CSFFC2-IN74E9-R21AB247-L1-WGS-0-V1_S1_L001_R1_001.fastq.gz input/fastqs/9988423398_0.R2.fastq.gz`
3. `ln -s /projects/ps-gleesonlab3/RCIGM_WGS_NTD_samples/202105_NTD_WGS_del22q11_subjects/210419_A01342_0041_AHYHG7DRXX/R21AB247-L1/UCSD-WGS-CSFFC2-IN74E9-R21AB247-L1-WGS-0-V1_S1_L002_R1_001.fastq.gz input/fastqs/9988423398_1.R1.fastq.gz`
4. `ln -s /projects/ps-gleesonlab3/RCIGM_WGS_NTD_samples/202105_NTD_WGS_del22q11_subjects/210419_A01342_0041_AHYHG7DRXX/R21AB247-L1/UCSD-WGS-CSFFC2-IN74E9-R21AB247-L1-WGS-0-V1_S1_L002_R1_001.fastq.gz input/fastqs/9988423398_1.R2.fastq.gz`

I assigned (arbitrarily) read group `0` to the first pair of FASTQs, and read group `1` to the second pair of FASTQs.

Now actually running the pipeline should be fairly simple.  If you have set up the `pbs-torque` cluster profile, all that is necessary is to run the command `snakemake --profile pbs-torque --use-conda` from within the `workflow` directory.  Log files will be saved to the `logs` directory and output files to the `output` directory.  You will need to archive these files to long-term storage after successfully running the pipeline.
---

# Notes:
1.  When installing the pipeline you will be prompted to ask which reference genome to use.  This is stored in `config/project.yaml` and can be changed between `hg38` and `grch37` if needed later.
2.  When installing the pipeline you will be prompted to ask whether the data is WES or WGS.  The only change that is made is that with WGS the DP (read depth) annotation is used for VQSR, while it is not for WES.  Please see [Which training sets arguments should I use for running VQSR?](https://gatk.broadinstitute.org/hc/en-us/articles/4402736812443-Which-training-sets-arguments-should-I-use-for-running-VQSR-) for the rationale.  In principle you can align and call variants on a mixture of exomes and genomes but there are biases introduced (and VQSR cannot be run properly). If you need to change this later you can either add or remove `DP` from the `annotation` section of `snp_recalibration` and/or `indel_recalibration` in `config/config.yaml`.
