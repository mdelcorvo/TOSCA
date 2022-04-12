# TOSCA workflow <img src="img/tosca_image.png" width="200" align="right" />


**TOSCA** (**T**umor **O**nly **S**omatic **CA**lling) is a [Snakemake workflow](https://snakemake.readthedocs.io/en/stable/index.html), aimed at performing a somatic variant calling (without matched normal samples) analysis in a reproducible, automated, and open source workflow.

TOSCA consists of a `Snakefile`, a set of [`conda`](https://conda.io/docs/) environment files (`envs/*.yaml`) a configuration file (`config/config.yaml`) and a set of `R` scripts. It is able to perform an end-to-end analysis, from raw read files, via quality checks, alignment and variant calling, to functional annotation, databases filtering, tumor purity and ploidy estimation and finally variant classification in whole exome / target sequencing data.  

By default, the pipeline performs only mandatory steps shown in the [diagram](img/dag_tosca.png) below. However, you can turn on optional rules in the `config/config.yaml` file. 

*Advanced use*: If you wish to custom your analysis by adding or removing some steps (e.g. VarScan over Mutect2 or Bowtie2 over BWA), you can use the software of your preference provided you have your own script(s), and change some lines within the `Snakefile`. If you think your "custom rule" might be of interest to a broader audience, let us know by opening an issue.


## Using the TOSCA workflow

We assume that you already have conda and Snakemake installed, otherwise you can easily install them with the following commands:

To install conda: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

To install Snakemake via conda: conda install -c conda-forge -c bioconda snakemake=6.10.0 snakemake-wrapper-utils mamba
```
To use TOSCA:

git clone https://github.com/mdelcorvo/TOSCA.git
#edit config and meta file
cd TOSCA && snakemake --use-conda --configfile config/config.yaml
```

More details about TOSCA workflow can be found in the [wiki](https://github.com/mdelcorvo/TOSCA/wiki).

## Workflow graph
![DAG](img/dag_tosca.png) 

Simplified directed acyclic graph (DAG) of the TOSCA workflow. Oval shapes are referred to data (input and intermediate files) while rectangles to processes. Raw input/output files are depicted as green while files independently downloaded by TOSCA are depicted as orange ovals. Red blocks are mandatory rules, while dashed lines and blue blocks are optional rules, controlled in the `config.yaml`. By default only mandatory rules are executed.

## Contributors
Current contributors include:

- [Marcello Del Corvo](https://github.com/mdelcorvo)
