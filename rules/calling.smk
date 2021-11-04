rule callable:
    input:
         ref=expand("resources/reference_genome/{ref}/{species}.fasta",ref=config["ref"]["build"],species=config["ref"]["species"]),
         bam=outputdir + "recal/{sample}.bam",
         bed=config["filtering"]["restrict-regions"],
         script = "scripts/callable_filt.R"
    output:
         out_bed=outputdir + "qc/coverage-stats/{sample}.bed",
         out_filt=outputdir + "qc/coverage-stats/{sample}.filt.bed",
         out_sum=outputdir + "qc/coverage-stats/{sample}.summary.txt",
         out_final=outputdir + "qc/coverage-stats/{sample}.final.bed"
    params:
         minD=config["filtering"]["min_depth"],
         gatk3=config["gatk3"]
    conda:
         "../envs/r4.yaml"
    benchmark:
        outputdir + "benchmarks/callable/{sample}.callable.benchmark.txt"
    shell:
        "java -jar {params.gatk3} -T CallableLoci -R {input.ref} -I {input.bam} --minDepth {params.minD} -summary {output.out_sum} -o {output.out_bed} -U ALLOW_SEQ_DICT_INCOMPATIBILITY -L {input.bed};" 
        "grep 'CALLABLE' {output.out_bed} > {output.out_filt};"      
        "Rscript --vanilla {input.script} {output.out_filt} {output.out_final}"
        
rule mutect2:
    input:
        ref=expand("resources/reference_genome/{ref}/{species}.fasta",ref=config["ref"]["build"],species=config["ref"]["species"]),
        bam = outputdir + "recal/{sample}.bam",
        pon= config["database_url"]["GRCh38"]["germline"]["PON"] if config["ref"]["build"]=='GRCh38' else config["database_url"]["GRCh37"]["germline"]["PON"],
        exac= config["database_url"]["GRCh38"]["germline"]["ExAC"] if config["ref"]["build"]=='GRCh38' else config["database_url"]["GRCh37"]["germline"]["ExAC"],
        bed=outputdir + "qc/coverage-stats/{sample}.final.bed"
    output:       
        vcf_raw = outputdir + "variant/{sample}.vcf",
        af_table= outputdir + "annotation/{sample}.af.table"
    message:
        "Implementing Mutect2"
    log:
        outputdir + "logs/mutect2/{sample}.log"
    benchmark:
        outputdir + "benchmarks/mutect2/{sample}.mutect2.benchmark.txt"
    conda:
        "../envs/gatk4.yaml"    
    shell:
        "gatk Mutect2 -R {input.ref} -I {input.bam} --germline-resource {input.exac} --panel-of-normals {input.pon} -L {input.bed} -O {output.vcf_raw}; "
        "gatk VariantsToTable -V {output.vcf_raw} -F CHROM -F POS -F TYPE -GF AF --show-filtered true -O {output.af_table}"
