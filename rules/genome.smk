rule get_genome:
    output:
        expand("resources/reference_genome/{ref}/{species}.fasta",ref=config["ref"]["build"],species=config["ref"]["species"])
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"]
    log:
        outputdir + "logs/ensembl/get_genome.log"    
    cache: True
    wrapper:
        "0.78.0/bio/reference/ensembl-sequence"

rule get_annotation:
    output:
        expand("resources/reference_genome/{ref}/{species}_annotation.gtf",ref=config["ref"]["build"],species=config["ref"]["species"])
    params:
        species=config["ref"]["species"],
        release=config["ref"]["release"] if config["ref"]["release"]=='GRCh38' else 87,
        build=config["ref"]["build"],
        fmt="gtf"
    cache: True  # save space and time with between workflow caching (see docs)
    wrapper:
        "0.78.0/bio/reference/ensembl-annotation"

rule gtf:
    input:
        gtf=expand("resources/reference_genome/{ref}/{species}_annotation.gtf",ref=config["ref"]["build"],species=config["ref"]["species"]),
        script = "scripts/exon_ranking.R"
    params:
        build=config["ref"]["build"]    
    output:
        rdata=expand("resources/reference_genome/{ref}/exon_ranking.RData",ref=config["ref"]["build"])
    conda:
        "../envs/r4.yaml"
    shell:
        "Rscript --vanilla {input.script} {input.gtf} {output.rdata} "

rule genome_dict:
    input:
         expand("resources/reference_genome/{ref}/{species}.fasta",ref=config["ref"]["build"],species=config["ref"]["species"])
    output:
         expand("resources/reference_genome/{ref}/{species}.dict",ref=config["ref"]["build"],species=config["ref"]["species"])
    conda:
        "../envs/samtools.yaml"     
    cache: True
    shell:
        "samtools dict {input} > {output} "
        
rule bwa_index:
    input:
        expand("resources/reference_genome/{ref}/{species}.fasta",ref=config["ref"]["build"],species=config["ref"]["species"])
    output:
        expand("resources/reference_genome/{ref}/{species}.fasta.{bwa}",ref=config["ref"]["build"],species=config["ref"]["species"],bwa=["amb", "ann" ,"bwt", "pac", "sa"])
    params:
        algorithm="bwtsw"
    log:
        outputdir + "logs/bwa_index/bwa_index.log"          
    resources:
        mem_mb=369000
    cache: True
    wrapper:
        "0.78.0/bio/bwa/index"        

rule genome_faidx:
    input:
         expand("resources/reference_genome/{ref}/{species}.fasta",ref=config["ref"]["build"],species=config["ref"]["species"]),
         expand("resources/reference_genome/{ref}/{species}.fasta.{bwa}",ref=config["ref"]["build"],species=config["ref"]["species"],bwa=["amb", "ann" ,"bwt", "pac", "sa"])
    output:
         expand("resources/reference_genome/{ref}/{species}.fasta.fai",ref=config["ref"]["build"],species=config["ref"]["species"])
    cache: True
    wrapper:
        "0.78.0/bio/samtools/faidx"

rule get_known_variation:
    input:
        # use fai to annotate contig lengths for GATK BQSR
        fai=expand("resources/reference_genome/{ref}/{species}.fasta.fai",ref=config["ref"]["build"],species=config["ref"]["species"])
    output:
        vcf=temp(expand("resources/database/{ref}/variation.vcf.gz",ref=config["ref"]["build"]))
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        type="all"
    cache: True
    wrapper:
        "0.78.0/bio/reference/ensembl-variation"

rule remove_iupac_codes:
    input:
        expand("resources/database/{ref}/variation.vcf.gz",ref=config["ref"]["build"])
    output:
        expand("resources/database/{ref}/variation.noiupac.vcf.gz",ref=config["ref"]["build"])
    conda:
        "../envs/rbt.yaml"
    cache: True
    shell:
        "rbt vcf-fix-iupac-alleles < {input} | bcftools view -Oz > {output}"
              
rule tabix_known_variants:
    input:
        expand("resources/database/{ref}/variation.noiupac.vcf.gz",ref=config["ref"]["build"])
    output:
        expand("resources/database/{ref}/variation.noiupac.vcf.gz.tbi",ref=config["ref"]["build"])
    params:
        "-p vcf"
    cache: True
    wrapper:
        "0.78.0/bio/tabix"   

rule snpeff_download:
    output:
        directory(expand("resources/database/snpEff/{snpeff_dir}",snpeff_dir=config["ref"]["snpeff_db"]))
    params:
         db=config["ref"]["snpeff_db"]   
    conda:
        "../envs/snp_eff.yaml"    
    shell:
        "snpEff download -v {params.db} -dataDir $(pwd)/resources/database/snpEff; "
