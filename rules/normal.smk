rule norm_trim_reads_pe:
    input:
        unpack(norm_fastq)
    output:
        r1=temp(normaldir + "trimmed/{id}_1_fastq.gz"),
        r2=temp(normaldir + "trimmed/{id}_2_fastq.gz"),
        r1_unpaired=temp(normaldir + "trimmed/{id}_1.unpaired.fastq.gz"),
        r2_unpaired=temp(normaldir + "trimmed/{id}_2.unpaired.fastq.gz")     
    params:
        **config["params"]["trimmomatic"]["pe"]
    threads:
        config["ncores"]    
    wrapper:
        "0.78.0/bio/trimmomatic/pe"
        
rule norm_map_reads:
    input:
        reads=norm_trimmed_reads,
        idx=rules.bwa_index.output,
        fai=rules.genome_faidx.output
    output:
        temp(normaldir + "mapped/{id}.sorted.bam")
    log:
        normaldir + "logs/bwa_mem/{id}.log"
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=norm_read_group,
        sorting="samtools",
        sort_order="coordinate"
    threads: config["ncores"]
    wrapper:
        "0.78.0/bio/bwa/mem"

rule norm_mark_duplicates:
    input:
        normaldir + "mapped/{id}.sorted.bam"
    output:
        bam=temp(normaldir + "dedup/{id}.bam"),
        metrics=normaldir + "qc/dedup/{id}.metrics.txt"
    log:
        outputdir + "logs/picard/dedup/{id}.log"
    params:
        "REMOVE_DUPLICATES=true OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT"
    resources:
        mem_mb=4096
    wrapper:
        "0.78.0/bio/picard/markduplicates"
        
rule norm_gatk_baserecalibrator:
    input:
        bam=normaldir + "dedup/{id}.bam",
        ref=expand("resources/reference_genome/{ref}/{species}.fasta",ref=config["ref"]["build"],species=config["ref"]["species"]),
        dict=expand("resources/reference_genome/{ref}/{species}.dict",ref=config["ref"]["build"],species=config["ref"]["species"]),
        known=expand("resources/database/{ref}/variation.noiupac.vcf.gz",ref=config["ref"]["build"]),
        tbi=expand("resources/database/{ref}/variation.noiupac.vcf.gz.tbi",ref=config["ref"]["build"])
    output:
        recal_table=normaldir + "recal/{id}.table"
    log:
        normaldir + "logs/gatk/bqsr/{id}.log"
    wrapper:
        "0.78.0/bio/gatk/baserecalibrator"      

rule norm_gatk_applybqsr:
    input:
        bam=normaldir + "dedup/{id}.bam",
        ref=expand("resources/reference_genome/{ref}/{species}.fasta",ref=config["ref"]["build"],species=config["ref"]["species"]),
        dict=expand("resources/reference_genome/{ref}/{species}.dict",ref=config["ref"]["build"],species=config["ref"]["species"]),
        recal_table=normaldir + "recal/{id}.table"
    output:
        bam=temp(normaldir + "recal/{id}.bam")
    log:
        normaldir + "logs/gatk/gatk_applybqsr/{id}.log"
    wrapper:
        "0.78.0/bio/gatk/applybqsr"

rule norm_mutect2:
    input:
        ref=expand("resources/reference_genome/{ref}/{species}.fasta",ref=config["ref"]["build"],species=config["ref"]["species"]),
        bam = normaldir + "recal/{id}.bam",
        exac= config["database_url"]["GRCh38"]["germline"]["ExAC"] if config["ref"]["build"]=='GRCh38' else config["database_url"]["GRCh37"]["germline"]["ExAC"],
        pon= config["database_url"]["GRCh38"]["germline"]["PON"] if config["ref"]["build"]=='GRCh38' else config["database_url"]["GRCh37"]["germline"]["PON"]
    output:
        vcf_raw= normaldir + "variant/{id}.vcf"
    message:
        "create a panel of normals"
    conda:
        "../envs/gatk4.yaml"    
    shell:
        "gatk Mutect2 -R {input.ref} -I {input.bam} --germline-resource {input.exac} --panel-of-normals {input.pon} -O {output.vcf_raw}; "

rule norm_bgzip:
    input:
        vcf= normaldir + "variant/{id}.vcf"
    output:
        vcf= normaldir + "variant/{id}.vcf.gz"
    conda:
        "../envs/r4.yaml"    
    shell:
        "bgzip -c {input.vcf} > {output.vcf};"
        "tabix -p vcf {output.vcf}"

rule bcftools_merge:
    input:
        vcf= expand(normaldir + "variant/{id}.vcf.gz",id=normals['id'])
    output:
        list= normaldir + "variant/normal_vcf_list.txt",
        norm= normaldir + "variant/normals.merged.min5.vcf"
    conda:
        "../envs/bcftools.yaml"    
    shell:
        "ls {input.vcf} > {output.list};"
        "bcftools merge --file-list {output.list} | bcftools view --min-ac=5 > {output.norm};"       
        
rule norm_bgzip_merged:
    input:
        vcf= normaldir + "variant/normals.merged.min5.vcf"
    output:
        vcf= normaldir + "variant/normals.merged.min5.vcf.gz"
    conda:
        "../envs/r4.yaml"    
    shell:
        "bgzip -c {input.vcf} > {output.vcf};"
        "tabix -p vcf {output.vcf}"
