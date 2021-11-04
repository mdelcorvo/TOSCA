# Copyright 2021 Marcello Del Corvo.
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed
# except according to those terms.


include: "rules/functions.smk"
# configfile: "./config.yaml"
outputdir = getpath(config["output"])

import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version
samples = pd.read_table(config["meta"], dtype=str).set_index(["sample"], drop=False)
if "meta_N" in config:
    normals = pd.read_table(config["meta_N"], dtype=str).set_index(["id"], drop=False)
    normaldir = getpath(config["normDB"])
    
rule all:
    input:
       outputdir + "results/Report.html",
       outputdir + "results/Somatic_Prediction.txt",
       expand(outputdir + "PureCN/{sample}.rds",sample=samples["sample"])  if "meta_N" in config else []    
           
           
##### Modules #####

include: "rules/genome.smk"
include: "rules/database.smk"
include: "rules/mapping.smk"
include: "rules/calling.smk"
include: "rules/filtering.smk"
include: "rules/qc.smk"
if "meta_N" in config:
    include: "rules/normal.smk"
    include: "rules/purecn.smk"
    
