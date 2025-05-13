# DNA QC rules
rule dna_chr_counts: 
    input: "{cell}.dedup.bam"
    output: "{cell}-chr_counts.txt"
    shell: 
        '''
        samtools view -q 20 {input} | \
        awk '$3~/^chr([0-9]+|[X,Y,M])$/{{sum[$3]++}}; END{{OFS="\t"; for (i in sum) print i,sum[i]}}' | \
        sort -k1,1V > {output}
        '''

rule dna_inserts: 
    input: "{cell}.dedup.bam"
    output: 
        pdf="{cell}-picard-inserts.pdf",
        txt="{cell}-picard-inserts.txt"
    shell: 
        '''
        picard CollectInsertSizeMetrics LEVEL="SAMPLE" I={input} O={output.txt} H={output.pdf} M=0.5
        if [[ $? == 0 && ! -f {output.txt} ]]; then
        touch {output.txt}
        touch {output.pdf}
        fi
        '''

rule qc_chrcounts:
    input: 
        lambda wildcards: expand(
            dna_dir + "/{cell}/{cell}-chr_counts.txt", 
            cell=get_cells_dna(wildcards.patient_id)
        )
    output: 
        out + "/{patient_id}/qc/{patient_id}-qc_dna_chrcounts.tsv"
    shell:
        '''
        awk '{{OFS="\t"; n=split(FILENAME, a, "/"); print a[n],$0}}' {input} > {output}
        '''

rule qc_dupstats:
    input: 
        lambda wildcards: expand(
            dna_dir + "/{cell}/{cell}-picard-mark_duplicates.txt", 
            cell=get_cells_dna(wildcards.patient_id)
        )
    output: 
        out + "/{patient_id}/qc/{patient_id}-qc_dna_dupstats.tsv"
    shell:
        '''
        awk 'FNR==7 && NR==7{{print}}FNR==8{{print}}' {input} > {output}
        '''

rule qc_inserts:
    input: 
        lambda wildcards: expand(
            dna_dir + "/{cell}/{cell}-picard-inserts.txt", 
            cell=get_cells_dna(wildcards.patient_id)
        )
    output: 
        out + "/{patient_id}/qc/{patient_id}-qc_dna_inserts.tsv"
    shell:
        '''
        awk '/^insert_size/{{f=1;next}}/^##/{{f=0}}f{{OFS="\t"; n=split(FILENAME, a, "/"); print a[n],$1,$2}}' {input} > {output}
        '''

checkpoint qc_dna:
    '''Collect DNA QC stats. Made into checkpoint to filter input to phasing etc'''
    input:
        seedfile=config["cells_dna"] if is_dna_analysis() else [],
        dupstats=out + "/{patient_id}/qc/{patient_id}-qc_dna_dupstats.tsv",
        chrfiles=out + "/{patient_id}/qc/{patient_id}-qc_dna_chrcounts.tsv",
        insfiles=out + "/{patient_id}/qc/{patient_id}-qc_dna_inserts.tsv",
        trim1=lambda wildcards: expand(
            os.path.join(dna_dir, "{cell}/trimmed/{cell}.1.trimming_report.txt"), 
            cell=get_cells_dna(wildcards.patient_id)
        ),
        trim2=lambda wildcards: expand(
            os.path.join(dna_dir, "{cell}/trimmed/{cell}.2.trimming_report.txt"), 
            cell=get_cells_dna(wildcards.patient_id)
        ),
    output:
        qc_tsv=out + "/{patient_id}/qc/{patient_id}-qc_dna.tsv",
        fq_plate=out + "/{patient_id}/qc/{patient_id}-qc_dna_fastq.pdf",
        map_plate=out + "/{patient_id}/qc/{patient_id}-qc_dna_maprate.pdf",
        dup_plate=out + "/{patient_id}/qc/{patient_id}-qc_dna_duprate.pdf",
        ins_plate=out + "/{patient_id}/qc/{patient_id}-qc_dna_inserts.pdf",
        ins_hist=out + "/{patient_id}/qc/{patient_id}-qc_dna_inserts_hist.pdf",
        fq_boxplot=out + "/{patient_id}/qc/{patient_id}-qc_dna_fastq_boxplot.pdf",
        dup_boxplot=out + "/{patient_id}/qc/{patient_id}-qc_dna_duprate_boxplot.pdf",
        ins_boxplot=out + "/{patient_id}/qc/{patient_id}-qc_dna_inserts_boxplot.pdf",
        counts_boxplot=out + "/{patient_id}/qc/{patient_id}-qc_dna_countsByChr.pdf",
        counts_boxplot_chrM=out + "/{patient_id}/qc/{patient_id}-qc_dna_chrM_boxplot.pdf",
        counts_plate_chrM=out + "/{patient_id}/qc/{patient_id}-qc_dna_chrM_plate.pdf",
        reads_vs_dups=out + "/{patient_id}/qc/{patient_id}-qc_dna_reads_vs_dups.pdf",
    params:
        dnaout=dna_dir,
        cells=lambda wildcards: get_cells_dna(wildcards.patient_id)
    script:
        "../scripts/qc_dna.R"

# RNA QC rules
rule rna_qc_dupstats:
    input: 
        lambda wildcards: expand(
            rna_dir + "/{cell}/{cell}-picard-mark_duplicates.txt", 
            cell=get_cells_rna(wildcards.patient_id)
        )
    output: 
        out + "/{patient_id}/qc/{patient_id}-qc_rna_dupstats.tsv",
    run:
        with open(output[0], "w") as out:
            for idx, file in enumerate(input):
                with open(file) as fp:
                    for i, line in enumerate(fp):
                        if idx == 0:
                            if i == 6:
                                out.write(line)
                        if i == 7:
                            pp = line.split("\t")
                            newline = os.path.basename(file) + "\t" + "\t".join(pp[1:])
                            out.write(newline)
                        elif i > 7:
                            break

rule qc_rna:
    input:
        counts=out + "/{patient_id}/{patient_id}-rna_counts.tsv.gz",
        dupstats=out + "/{patient_id}/qc/{patient_id}-qc_rna_dupstats.tsv",
        trim1=lambda wildcards: expand(
            rna_dir + "/{cell}/pre/{cell}.1.trimming_report.txt", 
            cell=get_cells_rna(wildcards.patient_id)
        ),
        trim2=lambda wildcards: expand(
            rna_dir + "/{cell}/pre/{cell}.2.trimming_report.txt", 
            cell=get_cells_rna(wildcards.patient_id)
        ),
        distr=lambda wildcards: expand(
            rna_dir + "/{cell}/{cell}.rseqc_readDistribution.txt", 
            cell=get_cells_rna(wildcards.patient_id)
        ) if config["rna"].get("detailed_qc", False) else []
    output:
        qc=out + "/{patient_id}/qc/{patient_id}-qc_rna.tsv",
        counts_plate=out + "/{patient_id}/qc/{patient_id}-qc_rna_counts.pdf",
        genes_plate=out + "/{patient_id}/qc/{patient_id}-qc_rna_genes.pdf",
        ercc_plate=out + "/{patient_id}/qc/{patient_id}-qc_rna_ercc.pdf",
        reads_vs_dups=out + "/{patient_id}/qc/{patient_id}-qc_rna_reads_vs_dups.pdf",
        boxplots=out + "/{patient_id}/qc/{patient_id}-qc_rna_boxplots.pdf",
        rseqc_distr=out + "/{patient_id}/qc/{patient_id}-qc_rna_distributions.pdf" if config["rna"].get("detailed_qc", False) else [],
    params: 
        mt_genes=config["ref"]["mt_genes"],
        cells=lambda wildcards: get_cells_rna(wildcards.patient_id)
    script:
        "../scripts/qc_rna.R"

# Shared metadata collection
rule metadata_collect:
    input:
        metafile=config["metadata"] if os.path.exists(config["metadata"]) else [],
        dnafile=config["cells_dna"] if is_dna_analysis() else [],
        rnafile=config["cells_rna"] if is_rna_analysis() else []
    output:
        meta_per_cell=out + "/{patient_id}/{patient_id}-metadata_long.tsv",
        meta_wsort=out + "/{patient_id}/{patient_id}-metadata_wsort.tsv",
        meta_wfacs=out + "/{patient_id}/{patient_id}-metadata_wfacs.tsv"
    params: 
        cells_rna=lambda wildcards: get_cells_rna(wildcards.patient_id),
        cells_dna=lambda wildcards: get_cells_dna(wildcards.patient_id)
    script:
        "../scripts/collect_metadata.R"