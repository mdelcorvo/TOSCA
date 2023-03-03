rule get_genome:
    output:
        expand("resources/reference_genome/{ref}/homo_sapiens.fasta",ref=config["ref"]["build"])
    params:
        species="homo_sapiens",
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"]
    log:
        outputdir + "logs/ensembl/get_genome.log"    
    cache: "omit-software"
    wrapper:
        "v1.23.3/bio/reference/ensembl-sequence"

rule get_annotation:
    output:
       expand("resources/reference_genome/{ref}/homo_sapiens.gtf",ref=config["ref"]["build"])
    params:
        species="homo_sapiens",
        release=config["ref"]["release"] if config["ref"]["release"]=='GRCh38' else 87,
        build=config["ref"]["build"]
    cache: "omit-software" 
    wrapper:
        "v1.23.3/bio/reference/ensembl-annotation"

rule gtf:
    input:
        gtf=expand("resources/reference_genome/{ref}/homo_sapiens.gtf",ref=config["ref"]["build"]),
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
         get_reference
    output:
         expand("resources/reference_genome/{ref}/homo_sapiens.dict",ref=config["ref"]["build"])
    conda:
        "../envs/samtools.yaml"     
    shell:
        "samtools dict {input} > {output} "
        
rule bwa_index:
    input:
        get_reference
    output:
        idx=expand("resources/reference_genome/{ref}/homo_sapiens.fasta.{bwa}",ref=config["ref"]["build"],bwa=["amb", "ann" ,"bwt", "pac", "sa"])
    params:
        algorithm="bwtsw"
    log:
        outputdir + "logs/bwa_index/bwa_index.log"          
    wrapper:
        "v1.23.3/bio/bwa/index"        

rule genome_faidx:
    input:
         get_reference,
         expand("resources/reference_genome/{ref}/homo_sapiens.fasta.{bwa}",ref=config["ref"]["build"],bwa=["amb", "ann" ,"bwt", "pac", "sa"])
    output:
         expand("resources/reference_genome/{ref}/homo_sapiens.fasta.fai",ref=config["ref"]["build"])
    wrapper:
        "v1.23.3/bio/samtools/faidx"

rule get_known_variation:
    input:
        # use fai to annotate contig lengths for GATK BQSR
        fai=expand("resources/reference_genome/{ref}/homo_sapiens.fasta.fai",ref=config["ref"]["build"])
    output:
        vcf=expand("resources/database/{ref}/variation.vcf.gz",ref=config["ref"]["build"])
    params:
        species="homo_sapiens",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        type="all"
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v1.23.3/bio/reference/ensembl-variation"
              
rule tabix_known_variants:
    input:
        expand("resources/database/{ref}/variation.vcf.gz",ref=config["ref"]["build"])
    output:
        expand("resources/database/{ref}/variation.vcf.gz.tbi",ref=config["ref"]["build"])
    params:
        "-p vcf"
    wrapper:
        "v1.23.3/bio/tabix/index"   

rule snpeff_download:
    output:
        directory(expand("resources/database/snpEff/{snpeff_dir}",snpeff_dir=config["ref"]["snpeff_db"]))
    params:
         db=config["ref"]["snpeff_db"]   
    conda:
        "../envs/snp_eff.yaml"    
    shell:
        "snpEff download -v {params.db} -dataDir $(pwd)/resources/database/snpEff; "
