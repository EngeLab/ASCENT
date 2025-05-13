## ASCENT
### Absolute Single-cell Copynumber Evaluation and subcloNe Topology

ASCENT is an end-to-end method for direct-tagmentation based single cell WGS (and scRNA-seq).

The input to ASCENT are demultiplexed (paired end) fastq files, and the output are absolute single cell copynumbers, clones, allele specific copynumbers per clone, and clonal topology. 

ASCENT is implemented within Snakemake and runs within conda environments to ensure reproducibility. Config files set parameters for each run, and is run with: 
```bash
snakemake --configfile=config.yaml --use-conda
```



