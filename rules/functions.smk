##### Helper functions #####

import re
def getpath(str):
	if str in ['', '.', './']:
		return ''
	if str.startswith('./'):
		regex = re.compile('^\./?')
		str = regex.sub('', str)
	if not str.endswith('/'):
		str += '/'
	return str

###########################

def get_fastq(wildcards):
    """Get fastq files of given sample."""
    fastqs = samples.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}

def is_single_end(sample):
    """Return True if sample is single end."""
    return pd.isnull(samples.loc[(sample), "fq2"])

def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample."""
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand(outputdir + "trimmed/{sample}_{group}_fastq.gz",
                      group=[1, 2], **wildcards)
    # single end sample
    return outputdir + "trimmed/{sample}.fastq.gz".format(**wildcards)

def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=wildcards.sample,
        platform=samples.loc[(wildcards.sample), "platform"])

###########################

def norm_fastq(wildcards):
    fastqs = normals.loc[(wildcards.id), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}

def is_single_norm(id):
    """Return True if sample is single end."""
    return pd.isnull(normals.loc[(id), "fq2"])
        
def norm_trimmed_reads(wildcards):
    if not is_single_norm(**wildcards):
        return expand(normaldir + "trimmed/{id}_{group}_fastq.gz",
                      group=[1, 2], **wildcards)
    # single end sample
    return normaldir + "trimmed/{id}.fastq.gz".format(**wildcards)                  

def norm_read_group(wildcards):
    return r"-R '@RG\tID:{id}\tSM:{id}\tPL:{platform}'".format(
        id=wildcards.id,
        platform=normals.loc[(wildcards.id), "platform"])
