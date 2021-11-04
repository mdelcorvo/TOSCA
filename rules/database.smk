##### Download database #####

rule download_database:
    output:
        expand(["resources/database/{ref}/temp_ESP6500SI.vcf.tar.gz",
        "resources/database/{ref}/1000G-phase_3.vcf.gz",
        "resources/database/{ref}/ClinVar.vcf.gz",
        "resources/database/{ref}/COSMIC.vcf.gz"],
        ref=config["ref"]["build"])
    params:
        build=config["ref"]["build"],
        Cosmic = config["database_url"]["GRCh38"]["somatic"]["Cosmic"] if config["ref"]["build"]=='GRCh38' else config["database_url"]["GRCh37"]["somatic"]["Cosmic"],
	Gen1K = config["database_url"]["GRCh38"]["germline"]["1000G"] if config["ref"]["build"]=='GRCh38' else config["database_url"]["GRCh37"]["germline"]["1000G"],
	ESP = config["database_url"]["GRCh38"]["germline"]["ESP"] if config["ref"]["build"]=='GRCh38' else config["database_url"]["GRCh37"]["germline"]["ESP"],
	Clinvar = config["database_url"]["GRCh38"]["germline"]["ClinVar"] if config["ref"]["build"]=='GRCh38' else config["database_url"]["GRCh37"]["germline"]["ClinVar"]
    cache: True
    shell:
        "curl -k -L  '{params.Cosmic}' > resources/database/{params.build}/COSMIC.vcf.gz; "
        "curl -k -L  '{params.Gen1K}' > resources/database/{params.build}/1000G-phase_3.vcf.gz; "
        "curl -k -L  '{params.ESP}' > resources/database/{params.build}/temp_ESP6500SI.vcf.tar.gz; "
        "curl -k -L  '{params.Clinvar}' > resources/database/{params.build}/ClinVar.vcf.gz; "
        "curl -k -L  '{params.Clinvar}.tbi' > resources/database/{params.build}/ClinVar.vcf.gz.tbi; "

rule download_mapping:
    output:
        expand("resources/database/{ref}/raw_mappability_{len}.bw",ref=config["ref"]["build"],len=config["kmer"])
    params:
        build=config["ref"]["build"],
      	mappability = config["mappability"]["GRCh38"] if config["ref"]["build"]=='GRCh38' else config["mappability"]["GRCh37"],
        kmer = config["kmer"]
    cache: True
    run:
        if config["ref"]["build"]=='GRCh38': 
          shell("curl -L  '{params.mappability}_{params.kmer}.bw' > resources/database/{params.build}/raw_mappability_{params.kmer}.bw; ")
        else:
          shell("curl -L  '{params.mappability}{params.kmer}mer.bigWig' > resources/database/{params.build}/raw_mappability_{params.kmer}.bw; ")
        	
rule ESP6500SI:
    input:
         expand("resources/database/{ref}/temp_ESP6500SI.vcf.tar.gz",ref=config["ref"]["build"]),
         script = "scripts/ESP6500SI.R"
    output:
        expand("resources/database/{ref}/ESP6500SI.vcf.gz",ref=config["ref"]["build"])
    params:
        build=config["ref"]["build"]      
    cache: True
    conda:
        "../envs/r4.yaml"
    shell:
        "tar -xzvf resources/database/{params.build}/temp_ESP6500SI.vcf.tar.gz -C resources/database/{params.build}; "
        "Rscript --vanilla {input.script} {params.build}; "
        "vcf-concat resources/database/{params.build}/*.vcf | sort -k1,1V -k2,2n | bgzip -c > resources/database/{params.build}/ESP6500SI.vcf.gz; "
        "rm resources/database/{params.build}/temp_ESP6500SI.vcf.tar.gz; "
        "rm resources/database/{params.build}/*.vcf; "
	
rule tabix_Cosmic:
    input:
        expand("resources/database/{ref}/COSMIC.vcf.gz",ref=config["ref"]["build"])
    output:
        expand("resources/database/{ref}/COSMIC.vcf.gz.tbi",ref=config["ref"]["build"])
    params:
        "-p vcf"
    cache: True
    wrapper:
        "0.78.0/bio/tabix" 

rule tabix_1000G:
    input:
        expand("resources/database/{ref}/1000G-phase_3.vcf.gz",ref=config["ref"]["build"])
    output:
       expand("resources/database/{ref}/1000G-phase_3.vcf.gz.tbi",ref=config["ref"]["build"])
    params:
        "-p vcf"
    cache: True
    wrapper:
        "0.78.0/bio/tabix" 

rule tabix_ESP:
    input:
         expand("resources/database/{ref}/ESP6500SI.vcf.gz",ref=config["ref"]["build"])
    output:
         expand("resources/database/{ref}/ESP6500SI.vcf.gz.tbi",ref=config["ref"]["build"])
    params:
        "-p vcf"
    cache: True
    wrapper:
        "0.78.0/bio/tabix"     

rule edit_mappability:
    input:
        raw_map=expand("resources/database/{ref}/raw_mappability_{len}.bw",ref=config["ref"]["build"],len=config["kmer"]),
        ref=expand("resources/reference_genome/{ref}/{species}.fasta",ref=config["ref"]["build"],species=config["ref"]["species"]),
        script = "scripts/edit_mappability.R"
    output:
        raw_bed= temp(expand("resources/database/{ref}/raw_mappability_{len}.bed",ref=config["ref"]["build"],len=config["kmer"])),
        bed= temp(expand("resources/database/{ref}/mappability_{len}.bed",ref=config["ref"]["build"],len=config["kmer"])),
        chromsizes= temp(expand("resources/reference_genome/{ref}/sizes.genome",ref=config["ref"]["build"])),
        map= expand("resources/database/{ref}/mappability_{len}.bw",ref=config["ref"]["build"],len=config["kmer"])
    conda:
        "../envs/r4.yaml"
    shell:
        "bigWigToBedGraph {input.raw_map} {output.raw_bed}; "
        "Rscript {input.script} {output.raw_bed} {output.bed}; "
        "faidx {input.ref} -i chromsizes > {output.chromsizes}; "
        "bedGraphToBigWig {output.bed} {output.chromsizes} {output.map}; "
                
##########################
