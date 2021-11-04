rule snpeff:
    input:
        expand("resources/database/snpEff/{snpeff_dir}",snpeff_dir=config["ref"]["snpeff_db"]),
        vcf=outputdir + "variant/{sample}.vcf"
    output:
        call=outputdir + "annotation/{sample}.snpEff.vcf",
        csv=outputdir + "qc/snpeff/{sample}.snpEff.csv",
        html=outputdir + "qc/snpeff/{sample}.snpEff.html"
    params:
         ref=config["ref"]["snpeff_db"],
         mem_gb="4" 
    conda:
        "../envs/snp_eff.yaml"  
    benchmark:
        outputdir + "benchmarks/snpeff/{sample}.snpeff.benchmark.txt"      
    shell:
        "snpEff -Xmx{params.mem_gb}g -dataDir $(pwd)/resources/database/snpEff -v -csvStats {output.csv} -s {output.html} {params.ref} {input.vcf} > {output.call}; "
        
rule snpsift_annotate_1000G:
    input:
        call=outputdir + "variant/{sample}.vcf",
        database=expand("resources/database/{ref}/1000G-phase_3.vcf.gz",ref=config["ref"]["build"]),
        tabix=expand("resources/database/{ref}/1000G-phase_3.vcf.gz.tbi",ref=config["ref"]["build"])
    output:
        call=outputdir + "annotation/{sample}.1000G.vcf"
    resources:
        mem_mb=1024
    benchmark:
        outputdir + "benchmarks/1000G/{sample}.1000G.benchmark.txt"    
    wrapper:
        "0.78.0/bio/snpsift/annotate"

rule snpsift_annotate_ESP:
    input:
        call=outputdir + "variant/{sample}.vcf",
        database=expand("resources/database/{ref}/ESP6500SI.vcf.gz",ref=config["ref"]["build"]),
        tabix=expand("resources/database/{ref}/ESP6500SI.vcf.gz.tbi",ref=config["ref"]["build"])
    output:
        call=outputdir + "annotation/{sample}.ESP.vcf"
    resources:
        mem_mb=1024
    benchmark:
        outputdir + "benchmarks/ESP/{sample}.ESP.benchmark.txt"     
    wrapper:
        "0.78.0/bio/snpsift/annotate"

rule snpsift_annotate_ExAC:
    input:
        call=outputdir + "variant/{sample}.vcf",
        database=config["database_url"]["GRCh38"]["germline"]["ExAC"]
    output:
        call=outputdir + "annotation/{sample}.ExAC.vcf"
    resources:
        mem_mb=1024
    benchmark:
        outputdir + "benchmarks/ExAC/{sample}.ExAC.benchmark.txt"    
    wrapper:
        "0.78.0/bio/snpsift/annotate"

rule snpsift_annotate_dbSNP:
    input:
        call=outputdir + "variant/{sample}.vcf",
        database=expand("resources/database/{ref}/variation.noiupac.vcf.gz",ref=config["ref"]["build"])
    output:
        call=outputdir + "annotation/{sample}.dbSNP.vcf"
    resources:
        mem_mb=1024
    benchmark:
        outputdir + "benchmarks/dbSNP/{sample}.dbSNP.benchmark.txt"     
    wrapper:
        "0.78.0/bio/snpsift/annotate"

rule Cosmic_annotate:
    input:
        call=outputdir + "variant/{sample}.vcf",
        database=expand("resources/database/{ref}/COSMIC.vcf.gz",ref=config["ref"]["build"]),
        tabix=expand("resources/database/{ref}/COSMIC.vcf.gz.tbi",ref=config["ref"]["build"]),
        script = "scripts/COSMIC.R"
    output:
        call=outputdir + "annotation/{sample}.Cosmic.vcf"
    conda:
        "../envs/r4.yaml"
    benchmark:
        outputdir + "benchmarks/Cosmic/{sample}.Cosmic.benchmark.txt"    
    shell:
        "Rscript --vanilla {input.script} {input.call} {input.database} {output.call} "

rule snpsift_annotate_ClinVar:
    input:
        call=outputdir + "variant/{sample}.vcf",
        database=expand("resources/database/{ref}/ClinVar.vcf.gz",ref=config["ref"]["build"])
    output:
        call=outputdir + "annotation/{sample}.ClinVar.vcf"
    resources:
        mem_mb=1024
    benchmark:
        outputdir + "benchmarks/ClinVar/{sample}.ClinVar.benchmark.txt"      
    wrapper:
        "0.78.0/bio/snpsift/annotate"

rule vcf_coverage:
    input:
         call=outputdir + "variant/{sample}.vcf",
         bam=outputdir + "recal/{sample}.bam"
    params:
         config["bamstats"]     
    output:      
         cov=outputdir + "annotation/{sample}.cov"
    shell:   
        "java -jar {params} -B {input.call} {input.bam} > {output.cov}"
        
rule annotation_assembler:
    input:
        vcf = outputdir + "variant/{sample}.vcf",
        snpeff=outputdir + "annotation/{sample}.snpEff.vcf",
        gen1k=outputdir + "annotation/{sample}.1000G.vcf",
        esp=outputdir + "annotation/{sample}.ESP.vcf",
        exac=outputdir + "annotation/{sample}.ExAC.vcf",
        dbsnp=outputdir + "annotation/{sample}.dbSNP.vcf",
        cosmic=outputdir + "annotation/{sample}.Cosmic.vcf",
        clinvar=outputdir + "annotation/{sample}.ClinVar.vcf",
        cov=outputdir + "annotation/{sample}.cov",
        af_table= outputdir + "annotation/{sample}.af.table",
        gtf= expand("resources/reference_genome/{ref}/exon_ranking.RData",ref=config["ref"]["build"]),
        script = "scripts/merge_annotation.R"
    params:
         ref=config["ref"]["build"],
         depth=config["filtering"]["min_depth"],
         vaf=config["filtering"]["vaf"] 
    output:
        call=outputdir + "variant/{sample}.txt"
    conda:
        "../envs/r4.yaml"
    benchmark:
        outputdir + "benchmarks/SomCalling/{sample}.SomCalling.benchmark.txt"    
    shell:
        "Rscript --vanilla {input.script} {input.vcf} {input.snpeff} {input.gen1k} {input.esp} {input.exac} {input.dbsnp} {input.cosmic} {input.clinvar} {input.cov} {input.af_table} {input.gtf} {params.ref} {params.depth} {params.vaf} {output.call} "

rule merge_results:
    input:
        call=expand(outputdir + "variant/{sample}.txt",sample=samples['sample']),
        ps=expand(outputdir + "PureCN/{sample}.predictSomatic.txt",sample=samples['sample']) if "meta_N" in config else [],
        script = "scripts/merge_results.R"
    output:
        list= temp(outputdir + "results/results.list"),
        list_ps= temp(outputdir + "results/results_ps.list") if "meta_N" in config else [],
        res= outputdir + "results/Somatic_Prediction.txt"
    conda:
        "../envs/r4.yaml"    
    shell:
        "ls {input.call} > {output.list};ls {input.ps} > {output.list_ps};Rscript {input.script} {output.list} {output.res} {output.list_ps} " if "meta_N" in config else "ls {input.call} > {output.list};Rscript {input.script} {output.list} {output.res} {output.list_ps} "
