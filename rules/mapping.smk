rule trim_reads_pe:
    input:
        unpack(get_fastq)
    output:
        r1=temp(outputdir + "trimmed/{sample}_1_fastq.gz"),
        r2=temp(outputdir + "trimmed/{sample}_2_fastq.gz"),
        r1_unpaired=temp(outputdir + "trimmed/{sample}_1.unpaired.fastq.gz"),
        r2_unpaired=temp(outputdir + "trimmed/{sample}_2.unpaired.fastq.gz"),      
    params:
        **config["params"]["trimmomatic"]["pe"]
    threads:
        config["ncores"]    
    log:
        outputdir + "logs/trimmomatic/{sample}.log"
    benchmark:
        outputdir + "benchmarks/trimmomatic/{sample}.trim.benchmark.txt"
    wrapper:
        "0.78.0/bio/trimmomatic/pe"
        
rule map_reads:
    input:
        reads=get_trimmed_reads,
        idx=rules.bwa_index.output,
        fai=rules.genome_faidx.output
    output:
        temp(outputdir + "mapped/{sample}.sorted.bam")
    log:
        outputdir + "logs/bwa_mem/{sample}.log"
    benchmark:
        outputdir + "benchmarks/bwa_mem/{sample}.bwa.benchmark.txt"
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=get_read_group,
        sorting="samtools",
        sort_order="coordinate"
    threads: config["ncores"]
    wrapper:
        "0.78.0/bio/bwa/mem"

rule mark_duplicates:
    input:
        outputdir + "mapped/{sample}.sorted.bam"
    output:
        bam=temp(outputdir + "dedup/{sample}.bam"),
        metrics=outputdir + "qc/dedup/{sample}.metrics.txt"
    log:
        outputdir + "logs/picard/dedup/{sample}.log"
    benchmark:
        outputdir + "benchmarks/picard/{sample}.picard.benchmark.txt"
    params:
        "REMOVE_DUPLICATES=true OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT"
    resources:
        mem_mb=4096    
    wrapper:
        "0.78.0/bio/picard/markduplicates"
        
rule gatk_baserecalibrator:
    input:
        bam=outputdir + "dedup/{sample}.bam",
        ref=expand("resources/reference_genome/{ref}/{species}.fasta",ref=config["ref"]["build"],species=config["ref"]["species"]),
        dict=expand("resources/reference_genome/{ref}/{species}.dict",ref=config["ref"]["build"],species=config["ref"]["species"]),
        known=expand("resources/database/{ref}/variation.noiupac.vcf.gz",ref=config["ref"]["build"]),
        tbi=expand("resources/database/{ref}/variation.noiupac.vcf.gz.tbi",ref=config["ref"]["build"])
    output:
        recal_table=outputdir + "recal/{sample}.table"
    log:
        outputdir + "logs/gatk/bqsr/{sample}.log"
    benchmark:
        outputdir + "benchmarks/gatk/bqsr/{sample}.gatk.bqsr.benchmark.txt"
    wrapper:
        "0.78.0/bio/gatk/baserecalibrator"    

rule gatk_applybqsr:
    input:
        bam=outputdir + "dedup/{sample}.bam",
        ref=expand("resources/reference_genome/{ref}/{species}.fasta",ref=config["ref"]["build"],species=config["ref"]["species"]),
        dict=expand("resources/reference_genome/{ref}/{species}.dict",ref=config["ref"]["build"],species=config["ref"]["species"]),
        recal_table=outputdir + "recal/{sample}.table"
    output:
        bam=outputdir + "recal/{sample}.bam"
    log:
        outputdir + "logs/gatk/gatk_applybqsr/{sample}.log"
    benchmark:
        outputdir + "benchmarks/gatk/gatk_applybqsr/{sample}.gatk.applybqsr.benchmark.txt"
    wrapper:
        "0.78.0/bio/gatk/applybqsr"
