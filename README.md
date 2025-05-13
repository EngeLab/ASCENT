## ASCENT
### Absolute Single-cell Copynumber Evaluation and subcloNe Topology

ASCENT is an end-to-end method for direct-tagmentation based single cell WGS (and scRNA-seq).

The input to ASCENT are demultiplexed (paired end) fastq files, and the output are absolute single cell copynumbers, clones, allele specific copynumbers per clone, and clonal topology. 

ASCENT is implemented within Snakemake and runs within conda environments to ensure reproducibility. Config files set parameters for each run, and is run with: 
```bash
snakemake --configfile=config.yaml --use-conda
```


## Inputs
A configfile denoting parameters and paths for each run (see config/HCT.yaml)

A seedfile with the columns cell, fq1, fq2 with cell_id, and full paths to fq1 and fq2. 

In the resource file there are a lot of the files needed to run the pipeline, however larger files are available through figshare, alternatively any other reference can be used. 

Specifically about normal cells in resources/figshare: 
One of the normalisations performed in ASCENT depends on having normal cells from the same sequencer, with the same read length as your experiment was performed on. This is not necessary, but if those are available it is recommended, as it decreases noise and increases both segmentation accuracy and scaling accuracy. 
In the resources folder we include bincounts that can be used for experiments run on Nextseq 550 (2x37bp) and Novaseq 6000 (1x150bp), with resolutions 40kb-500kb. For higher resolution (as is recommended) we share the files through figshare. 

Where we got the reference genomes: 
We aligned against Grch38, d1, vd1 obtained from here https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files. 

We used an exclusion list from encode from here https://www.encodeproject.org/files/ENCFF356LFX/, for genome build Grch38. 

To perform phasing and haplotyping we downloaded data from here gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2 

Paths need to be adjusted in the rule phasing.smk 


