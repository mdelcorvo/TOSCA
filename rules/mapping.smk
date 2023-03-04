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
        "v1.23.3/bio/trimmomatic/pe"
        
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
        "v1.23.3/bio/bwa/mem"

rule mark_duplicates:
    input:
        bams=outputdir + "mapped/{sample}.sorted.bam"
    output:
        bam=outputdir + "dedup/{sample}.bam",
        metrics=outputdir + "qc/dedup/{sample}.metrics.txt"
    log:
        outputdir + "logs/picard/dedup/{sample}.log"
    benchmark:
        outputdir + "benchmarks/picard/{sample}.picard.benchmark.txt"
    params:
        extra="--REMOVE_DUPLICATES true --OPTICAL_DUPLICATE_PIXEL_DISTANCE 100 --CREATE_INDEX true --VALIDATION_STRINGENCY LENIENT"
    resources:
        mem_mb=4096    
    wrapper:
        "v1.23.3/bio/picard/markduplicates"
