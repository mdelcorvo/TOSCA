# TOSCA workflow


**TOSCA** (**T**umor **O**nly **S**omatic **CA**lling) is a [Snakemake workflow](https://snakemake.readthedocs.io/en/stable/index.html), aimed at performing a somatic variant calling (without matched normal samples) workflow in a reproducible, automated, and partially contained manner. It is implemented such that alternative or similar analysis can be added or removed. 

TOSCA consists of a `Snakefile`, a [`conda`](https://conda.io/docs/) environment file (`envs/environment.yaml`) a configuration file (`config.yaml`) and a set of `R` scripts, to perform quality control, preprocessing, variant calling, functional annotation and databases filtering of whole exome / target sequencing.  

By default, the pipeline performs all the steps shown in the [diagram](img/dag_nice3.png) below. However, you can turn off any combination of the light-colored steps in the `config.yaml` file. 

*Advanced use*: If you prefer other software to run one of the outlined steps (e.g. `DESeq2` over `edgeR`, or `kallisto` over `Salmon`), you can use the software of your preference provided you have your own script(s), and change some lines within the `Snakefile`. If you think your "custom rule" might be of use to a broader audience, let us know by opening an issue.


## Using the TOSCA workflow

Assuming that snakemake and conda are installed (and your system has the necessary libraries to compile R packages), you can use the following commands on a test dataset:

```
git clone https://github.com/csoneson/TOSCA.git
cd TOSCA && snakemake --use-conda
```

To use the TOSCA workflow on your own data, follow the steps outlined in the [wiki](https://github.com/mdelcorvo/TOSCA/wiki).


## Contributors
Current contributors include:

- [Marcello Del Corvo](https://github.com/mdelcorvo)

