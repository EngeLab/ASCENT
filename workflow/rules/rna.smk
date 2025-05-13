# workflow/rules/rna.smk

# RNA preprocessing
rule rna_adapter_trim:
    input:
        fq1=get_rna_fq1,
        fq2=get_rna_fq2
    output:
        fq1trim=temp(rna_dir + "/{cell}/pre/{cell}_val_1.fq.gz"),
        fq2trim=temp(rna_dir + "/{cell}/pre/{cell}_val_2.fq.gz"),
        report1=rna_dir + "/{cell}/pre/{cell}.1.trimming_report.txt",
        report2=rna_dir + "/{cell}/pre/{cell}.2.trimming_report.txt"
    params:
        merged1=rna_dir + "/{cell}/pre/{cell}.1.merged.fq.gz",
        merged2=rna_dir + "/{cell}/pre/{cell}.2.merged.fq.gz"
    log: 
        rna_dir + "/{cell}/log/{cell}.trimgalore.log"
    shell:
        '''
        outdir=$(dirname {output[0]})
        mkdir -vp $outdir
        fq1ar=({input.fq1})
        fq2ar=({input.fq2})
        [[ ${{#fq1ar[@]}} -gt 1 ]] && cat ${{fq1ar[@]}} > {params.merged1} && fq1={params.merged1} || fq1={input.fq1}
        [[ ${{#fq2ar[@]}} -gt 1 ]] && cat ${{fq2ar[@]}} > {params.merged2} && fq2={params.merged2} || fq2={input.fq2}
        trim_galore --paired --nextera --nextseq 20 -e 0.1 --basename {wildcards.cell} -o $outdir $fq1 $fq2 2> {log}
        mv -v $outdir/$(basename $fq1)_trimming_report.txt {output.report1}
        mv -v $outdir/$(basename $fq2)_trimming_report.txt {output.report2}
        rm -fv {params.merged1} {params.merged2}
        '''

rule rna_align:
    input:
        fq1=rna_dir + "/{cell}/pre/{cell}_val_1.fq.gz",
        fq2=rna_dir + "/{cell}/pre/{cell}_val_2.fq.gz"
    output:
        bam=temp(rna_dir + "/{cell}/align/Aligned.sortedByCoord.out.bam"),
        counts=rna_dir + "/{cell}/align/ReadsPerGene.out.tab"
    params:
        ref_dir=config["ref"]["rna_ref_dir"]
    threads: 2
    shell:
        '''
        outdir=$(dirname {output.bam})
        mkdir -vp $outdir
        STAR \
        --genomeDir {params.ref_dir} \
        --readFilesIn {input.fq1} {input.fq2} \
        --readFilesCommand zcat \
        --outFileNamePrefix "$outdir"/ \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outSAMstrandField intronMotif \
        --runThreadN {threads} \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI NM MD RG \
        --outSAMattrRGline ID:{wildcards.cell} LB:{wildcards.cell} SM:{wildcards.cell} \
        --genomeLoad LoadAndKeep \
        --limitBAMsortRAM 12000000000 \
        --quantMode GeneCounts
        '''

# Deduplication and indexing
rule rna_dedup:
    input: 
        rna_dir + "/{cell}/align/Aligned.sortedByCoord.out.bam"
    output:
        bam=rna_dir + "/{cell}/{cell}.sorted.bam",
        qc=rna_dir + "/{cell}/{cell}-picard-mark_duplicates.txt"
    params:
        tmp=config["path"]["temp"]
    threads: 2
    log: 
        rna_dir + "/{cell}/log/{cell}.picard.log"
    shell:
        '''
        picard MarkDuplicates \
        MAX_RECORDS_IN_RAM=5000000 \
        TMP_DIR={params.tmp} \
        INPUT={input} \
        OUTPUT={output.bam} \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=LENIENT \
        METRICS_FILE={output.qc} 2> {log}
        '''

rule rna_index_bam:
    input: 
        "{cell}.sorted.bam"
    output: 
        "{cell}.sorted.bam.bai"
    shell: 
        "samtools index {input}"

# Read counting
rule rna_count:
    input: 
        "{cell}.sorted.bam"
    output: 
        "{cell}.counts.txt"
    params:
        gtf=config["ref"]["rna_gtf"]
    log: 
        "{cell}.counts.log"
    shell:
        '''
        featureCounts -p -Q 30 -O --countReadPairs \
        -t exon -g gene_id -a {params.gtf} -o {output} {input} 2> {log}
        '''

# Collect count matrices

rule rna_collect:
    input: 
        lambda wildcards: expand(
            os.path.join(rna_dir,"{cell}/{cell}.counts.txt"), 
            cell=get_cells_rna(wildcards.patient_id),
        )
    output: out + "/{patient_id}/{patient_id}-rna_counts.tsv.gz"
    run:
        import gzip
        import os
        with gzip.open(output[0], 'wt') as writer:
            print("gene_id\t" + "\t".join(get_cells_rna(wildcards.patient_id)), file=writer)
            readers = [open(filename) for filename in input]
            for j, lines in enumerate(zip(*readers)):
                if j > 1:
                    for i, line in enumerate(lines):
                        line_list = line.split("\t")
                        if i == 0:
                            counts = [line_list[0], line_list[-1].strip()]
                        else:
                            counts.append(line_list[-1].strip())
                    print('\t'.join(counts), file=writer)

# rule rna_collect_star:
#     input: 
#         expand(rna_dir + "/{cell}/align/ReadsPerGene.out.tab",
#                cell=list(dict.fromkeys(cells_rna.cell)))
#     output: 
#         out + "/rna_counts.star.tsv.gz"
#     run:
#         import gzip
#         import os
#         with gzip.open(output[0], 'wt') as writer:
#             print("gene_id\t" + "\t".join(list(dict.fromkeys(cells_rna.cell))), file=writer)
#             readers = [open(filename) for filename in input]
#             for j, lines in enumerate(zip(*readers)):
#                 if j > 3:
#                     for i, line in enumerate(lines):
#                         line_list = line.split("\t")
#                         if i == 0:
#                             counts = [line_list[0], line_list[1].strip()]
#                         else:
#                             counts.append(line_list[1].strip())
#                     print('\t'.join(counts), file=writer)

# Quality Control
rule rna_rseqc:
    input: 
        "{cell}.sorted.bam"
    output: 
        "{cell}.rseqc_readDistribution.txt"
    params:
        genes=config["ref"]["genes"]
    conda: 
        "../envs/rseqc.yaml"
    shell: 
        "read_distribution.py -i {input} -r {params.genes} > {output}"

# RNA analysis
rule rna_load:
    input:
        counts=out + "/{patient_id}/{patient_id}-rna_counts.tsv.gz",
        meta_wsort=out + "/{patient_id}/{patient_id}-metadata_wsort.tsv",
        meta_wfacs=out + "/{patient_id}/{patient_id}-metadata_wfacs.tsv"
    output:
        seurat_obj=out + "/{patient_id}/{patient_id}-rna_seurat.Rds",
        dimplot=out + "/{patient_id}/{patient_id}-rna_dimplot.pdf",
        plot_qc=out + "/{patient_id}/{patient_id}-rna_qc.pdf",
        rna_phase=out + "/{patient_id}/{patient_id}-rna_phases.txt"
    params:
        min_count=config["rna"]["min_count"],
        min_genes=config["rna"]["min_genes"],
        quantile_cut=config["rna"]["quantile_cut"],
        quantile_gene=config["rna"]["quantile_gene"],
        rna_perplexity=30,
        cc_expr_cutoff=4.5,
        plot_qc=True,
        mt_genes=config["ref"]["mt_genes"],
        gene_info=config["ref"]["gene_info"]
    script: 
        "../scripts/process_rna.R"