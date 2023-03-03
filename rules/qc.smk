rule fastqc:
    input:
        outputdir + "trimmed/{sample}_{fq}_fastq.gz"
    output:
        html=outputdir + "qc/fastqc/{sample}_{fq}.html",
        zip=outputdir + "qc/fastqc/{sample}_{fq}_fastqc.zip"
    threads: 8
    log:
        outputdir + "logs/fastqc/{sample}_{fq}.log"
    benchmark:
        outputdir + "benchmarks/fastqc/{sample}_{fq}.fastqc.benchmark.txt"
    wrapper:
        "v1.23.3/bio/fastqc"

rule samtools_stats:
    input:
        bam=outputdir + "dedup/{sample}.bam"
    output:
        outputdir + "qc/samtools-stats/{sample}.txt"
    log:
        outputdir + "logs/samtools-stats/{sample}.log"
    benchmark:
        outputdir + "benchmarks/samtools-stats/{sample}.samtools-stats.benchmark.txt"
    wrapper:
        "v1.23.3/bio/samtools/stats"

rule mosdepth:
    input:
        bam=outputdir + "dedup/{sample}.bam",
        bed=config["filtering"]["restrict-regions"],
    output:
        outputdir + "qc/mosdepth/{sample}.mosdepth.global.dist.txt",
        outputdir + "qc/mosdepth/{sample}.mosdepth.region.dist.txt",
        outputdir + "qc/mosdepth/{sample}.regions.bed.gz",
        summary= outputdir + "qc/mosdepth/{sample}.mosdepth.summary.txt"
    log:
        outputdir + "logs/mosdepth/{sample}.log"
    benchmark:
        outputdir + "benchmarks/mosdepth/{sample}.mosdepth.benchmark.txt"
    threads: config["ncores"]
    wrapper:
        "v1.23.3/bio/mosdepth"
        
rule multiqc:
    input:
        expand([outputdir + "qc/fastqc/{sample}_{fq}_fastqc.zip",
                outputdir + "qc/samtools-stats/{sample}.txt",
                outputdir + "qc/dedup/{sample}.metrics.txt",
                outputdir + "qc/snpeff/{sample}.snpEff.csv",
                outputdir + "qc/mosdepth/{sample}.mosdepth.region.dist.txt"],
                sample=samples[samples.columns[0]],fq=["1", "2"])
    output:
        report(outputdir + "results/Report.html", caption="multiqc.rst", category="Quality control")
    log:
        outputdir + "logs/multiqc.log"
    benchmark:
        outputdir + "benchmarks/multiqc/multiqc.benchmark.txt"
    wrapper:
        "v1.23.3/bio/multiqc"
