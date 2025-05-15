## ASCENT
### Absolute Single-cell Copynumber Evaluation and subcloNe Topology

ASCENT is an end-to-end method for direct-tagmentation based single-cell WGS (and scRNA-seq).

The input to ASCENT are demultiplexed (paired end) fastq files, and the output are absolute single cell copynumbers, clones, allele specific copynumbers per clone, and clonal topology. 

ASCENT is implemented within Snakemake and runs within conda environments to ensure reproducibility. 

To run ASCENT create a dntr conda environment 

```bash
conda env create -f workflow/envs/dntr.yaml
```

Then edit a config file to set parameters for each run, and run with:

```bash
snakemake --configfile=config.yaml --use-conda
```




## Input
The input to ASCENT are single-cell WGS, demultiplexed fastq files. If paired scRNA-seq has been performed (like in DNTR-seq) then demultiplexed scRNA-seq fastq files can also be input. RNA-seq data is aligned and processed to identify cells in S or G2/M phase. If S/G2/M cells are identified they are removed before segmentation. 

To run ASCENT the user needs to update a config file denoting parameters and paths (see config/HCT.yaml) and a seedfile listing cell ids and paths to fq1 and fq2, and update the paths to haplotype and SNP information in the rule phasing.smk  

The resources needed to run ASCENT are in one of three places: The resources folder, figshare or created by the user by running Create_Filtered_Snps or following the instructions in Download_Phased_Haplotypes or CreateStarIndex

These resources are reference genomes, snps, haplotype information, exclusion lists and files needed for normalisation (gc, map and normal cells). 

One of the normalisations performed in ASCENT depends on having normal cells from the same sequencer, with the same read length as your experiment was performed on. This is not necessary, but if those are available it is recommended, as it decreases noise and increases both segmentation accuracy and scaling accuracy. 
In the resources folder we include bincounts that can be used for experiments run on Nextseq 550 (2x37bp) and Novaseq 6000 (1x150bp), with resolutions 40kb-500kb. For higher resolution (as is recommended) we share the files through figshare. 

Reference files: 

We aligned against Grch38, d1, vd1 obtained from here https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files. 

We used an exclusion list from encode from here https://www.encodeproject.org/files/ENCFF356LFX/, for genome build Grch38. 

To perform phasing and haplotyping we downloaded data from here gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2 

## Output

The output of ASCENT is organised in a result folder. All QC metrics are within the file results/patient/qc, all clone metrics are within the file results/patieng/clones, allele specific output is within the results/patient/baf folder and medicc output is in results/patient/medicc 

The clone folder includes multiple intermediate stages, at which the data should be scrutinized. We recommend running ultra-low coverage (60-500k reads per cell) data at 10-40kb resolution with gamma values for mpcf between 0.5-10. Both of these parameters are tunable within the config file. The gamma parameter is multiplied with number of samples within the mpcf function, so if there are many cells in a sample, the gamma value will be effectively higher, than in a sample with few cells. Therefore it is recommended to run 3-5 different gammas, and do a visual inspection of CNV heatmaps (result/patient/clones/patient-heatmap_revised-umap-g*-*.png) before continuing to refinement. Clonal refinement runs automatically within ASCENT, with a gamma value set in the config file. However it is often a good idea to run that code in a more interactive manner outside the pipeline, until a desired outcome is reached. By using the function plot_clone_detail(clone, region) the user can get a good idea if the gamma parameter is correct or not, or if some clones are a mixture of clones. 

Medicc2 runs automatically within ASCENT, but to root the tree correctly in the "normal" clone from the current sample it can be run outside the pipeline to indicate which clone is normal.  



