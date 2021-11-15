rule interval_file:
    input:
        baits=config["baits"],
        ref=expand("resources/reference_genome/{ref}/{species}.fasta",ref=config["ref"]["build"],species=config["ref"]["species"]),
        mappability= expand("resources/database/{ref}/mappability_{len}.bw",ref=config["ref"]["build"],len=config["kmer"]),
        gtf=expand("resources/reference_genome/{ref}/{species}_annotation.gtf",ref=config["ref"]["build"],species=config["ref"]["species"]),
        script = "scripts/PureCN/Intervalfile.R"
    params:
        min_target=100,
        min_off_target=20000  
    output:
        intervals= expand(normaldir + "intervals/baits_{build}_intervals.txt",build='hg38' if config["ref"]["build"]=='GRCh38' else 'hg19')
    conda:
        "../envs/r4.yaml"    
    shell:
        "Rscript {input.script} {input.baits} {input.ref} {input.mappability} {input.gtf} {params.min_target} {params.min_off_target} {output.intervals} "
        
rule norm_cov:
    input:
        intervals= expand(normaldir + "intervals/baits_{build}_intervals.txt",build='hg38' if config["ref"]["build"]=='GRCh38' else 'hg19'),
        bam = normaldir + "recal/{id}.bam",
        script = "scripts/PureCN/Coverage.R" 
    output:
        coverage_raw= normaldir + "cov_normal/{id}_coverage.txt",
        coverage_loess= normaldir + "cov_normal/{id}_coverage_loess.txt"
    conda:
        "../envs/r4.yaml"    
    shell:
        "Rscript {input.script} {input.intervals} {input.bam} {output.coverage_raw} {output.coverage_loess} " 

rule tum_cov:
    input:
        intervals= expand(normaldir + "intervals/baits_{build}_intervals.txt",build='hg38' if config["ref"]["build"]=='GRCh38' else 'hg19'),
        bam = outputdir + "recal/{sample}.bam",
        script = "scripts/PureCN/Coverage.R" 
    output:
        coverage_raw= outputdir + "cov_tumor/{sample}_coverage.txt",
        coverage_loess= outputdir + "cov_tumor/{sample}_coverage_loess.txt"
    conda:
        "../envs/r4.yaml"    
    shell:
        "Rscript {input.script} {input.intervals} {input.bam} {output.coverage_raw} {output.coverage_loess} " 
        
rule norm_db:
    input:
        cov_norm= expand(normaldir + "cov_normal/{id}_coverage_loess.txt",id=normals['id']),
        vcf= normaldir + "variant/normals.merged.min5.vcf.gz",
        script = "scripts/PureCN/NormalDB.R" 
    output:
        list= normaldir + "cov_normal/normal_coverage.list",
        rds= expand(normaldir + "normalDB/normalDB_{build}.rds",build='hg38' if config["ref"]["build"]=='GRCh38' else 'hg19'),
        rds_bias= expand(normaldir + "normalDB/mapping_bias_{build}.rds",build='hg38' if config["ref"]["build"]=='GRCh38' else 'hg19')
    params:
        ref='hg38' if config["ref"]["build"]=='GRCh38' else 'hg19'
    conda:
        "../envs/r4.yaml"    
    shell:
        "ls {input.cov_norm} > {output.list};"
        "Rscript {input.script} {input.vcf} {output.list} {output.rds} {output.rds_bias} {params.ref};"               

rule purecn:
    input:
        call=outputdir + "variant/{sample}.txt",
        rds= expand(normaldir + "normalDB/normalDB_{build}.rds",build='hg38' if config["ref"]["build"]=='GRCh38' else 'hg19'),
        cov_tum= outputdir + "cov_tumor/{sample}_coverage_loess.txt",
        vcf_tum = outputdir + "variant/{sample}.vcf",
        intervals= expand(normaldir + "intervals/baits_{build}_intervals.txt",build='hg38' if config["ref"]["build"]=='GRCh38' else 'hg19'),
        rds_bias= expand(normaldir + "normalDB/mapping_bias_{build}.rds",build='hg38' if config["ref"]["build"]=='GRCh38' else 'hg19'),
        script = "scripts/PureCN/PureCN.R"
    params:
        ref='hg38' if config["ref"]["build"]=='GRCh38' else 'hg19',
        sample = "{sample}"
    conda:
        "../envs/r4.yaml"
    output:
        result_rds= outputdir + "PureCN/{sample}.rds",
        ps= outputdir + "PureCN/{sample}.predictSomatic.txt",
        plot_output=outputdir + "PureCN/{sample}.pdf" 
    benchmark:
        outputdir + "benchmarks/purecn/{sample}.purecn.benchmark.txt"    
    shell:
        "Rscript {input.script} {input.rds} {input.cov_tum} {input.vcf_tum} {input.intervals} {input.rds_bias} {params.ref} {params.sample} {output.result_rds} {output.ps} {output.plot_output} "   
