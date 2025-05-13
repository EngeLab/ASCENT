# Bin creation and preparation
rule create_fixed_bins:
    input: 
    output: 
        "resources/fixed-{binsize}_arms.bed"
    params:
        genome=config["ref"]["chr_list"]
    shell:
        '''
        bedtools makewindows -g {params.genome} -w {wildcards.binsize} > {output}
        '''

rule fix_fixed_bins:
    input:
        bed="resources/fixed-{binsize}_arms.bed",
        arms="grch38-chrarms_start_end.tsv"
    output:
        bed="resources/fixed-{binsize}.bed"
    script: "../scripts/combine_chr_arms.R"


rule create_bin_tracks:
    '''Calculate fractional mappability and GC over fixed bins'''
    input:
        bins="resources/fixed-{binsize}.bed"
    output:
        map="resources/fixed-{binsize}.map.txt",
        gc="resources/fixed-{binsize}.gc.txt"
    params:
        ref=config["ref"]["fasta"],
        maptrack=config["ref"]["maptrack"]
    shell:
        '''
        bedtools coverage -a {input.bins} -b {params.maptrack} | cut -f7 > {output.map}
        bedtools nuc -fi {params.ref} -bed {input.bins} | cut -f5 | tail -n+2 > {output.gc}
        '''

# DNA preprocessing rules
rule dna_adapter_trim:
    '''
    Trim adapter sequences and poor quality bases
    Nextera flag == CTGTCTCTTATA
    If multiple fastqs available, will merge them here by concatenation of R1 and R2 respectively
    (by making fq1 and fq2 into arrays, and counting > 1 then merge them and change input name)
    '''
    input:
        fq1=get_dna_fq1,
        fq2=get_dna_fq2
    output:
        fq1trim=temp(dna_dir + "/{cell}/trimmed/{cell}_val_1.fq.gz"),
        fq2trim=temp(dna_dir + "/{cell}/trimmed/{cell}_val_2.fq.gz"),
        report1=dna_dir + "/{cell}/trimmed/{cell}.1.trimming_report.txt",
        report2=dna_dir + "/{cell}/trimmed/{cell}.2.trimming_report.txt"
    params:
        merged1=dna_dir + "/{cell}/trimmed/{cell}.1.merged.fq.gz",
        merged2=dna_dir + "/{cell}/trimmed/{cell}.2.merged.fq.gz"
    shell:
        '''
        outdir=$(dirname {output[0]})
        mkdir -vp $outdir
        fq1ar=({input.fq1})
        fq2ar=({input.fq2})
        [[ ${{#fq1ar[@]}} -gt 1 ]] && cat ${{fq1ar[@]}} > {params.merged1} && fq1={params.merged1} || fq1={input.fq1}
        [[ ${{#fq2ar[@]}} -gt 1 ]] && cat ${{fq2ar[@]}} > {params.merged2} && fq2={params.merged2} || fq2={input.fq2}
        trim_galore --paired --nextera --nextseq 20 -e 0.1 --basename {wildcards.cell} -o $outdir $fq1 $fq2
        mv -v $outdir/$(basename $fq1)_trimming_report.txt {output.report1}
        mv -v $outdir/$(basename $fq2)_trimming_report.txt {output.report2}
        rm -fv {params.merged1} {params.merged2}
        '''

rule align_bwa:
    '''Align paired DNA reads with bwa'''
    input:
        fq1=os.path.join(dna_dir,"{cell}/trimmed/{cell}_val_1.fq.gz"),
        fq2=os.path.join(dna_dir,"{cell}/trimmed/{cell}_val_2.fq.gz")
    params:
        rg="'@RG\\tID:{cell}\\tPL:ILLUMINA\\tPU:{cell}\\tLB:{cell}\\tSM:{cell}'",
        ref=config["ref"]["fasta"]
    output:
        temp(os.path.join(dna_dir,"{cell}/{cell}.aligned.bam"))
    threads: 1
    shell:
        '''
        bwa mem -t {threads} -M -k 25 -R {params.rg} {params.ref} {input.fq1} {input.fq2} | \
        samtools sort -@{threads} -O BAM -o {output} -
        '''

rule dna_dedup:
    input: "{cell}.aligned.bam"
    output:
        bam="{cell}.dedup.bam",
        qc="{cell}-picard-mark_duplicates.txt"
    params:
        tmp=config["path"]["temp"]
    threads: 2
    shell:
        '''
        picard MarkDuplicates \
        MAX_RECORDS_IN_RAM=5000000 \
        TMP_DIR={params.tmp} \
        INPUT={input} \
        OUTPUT={output.bam} \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=LENIENT \
        METRICS_FILE={output.qc}
        '''

rule index_bam:
    input: "{cell}.dedup.bam"
    output: "{cell}.dedup.bam.bai"
    shell: "samtools index {input}"

rule index_bam_compl:
    input: 
        lambda wildcards: expand(
            os.path.join(dna_dir,"{cell}/{cell}.dedup.bam.bai"), 
            cell=get_cells_dna(wildcards.patient_id)
            )
    output: out + "/{patient_id}/dna_bam_index_done"
    shell: "touch {output}"

# Fragment creation and binning
rule create_bed:
    input: "{cell}.dedup.bam"
    output: "{cell}.bed.gz"
    params:
        blacklist=config["ref"]["blacklist"],
        min_mapq=config["dna"]["min_mapq"],
        max_insert=config["dna"]["max_insert"],
        min_insert=config["dna"]["min_insert"],
        temp=config["path"]["temp"]
    shell:
        '''
        samtools sort -T {params.temp} -n -m 4G {input} | \
        samtools view -q {params.min_mapq} -b -f 0x2 -F 0x100 - | \
        bedtools bamtobed -bedpe -i stdin | \
        awk -F'\t' '$1~/^chr([0-9]+|[XY])$/ && $1==$4 {{ if($9=="-") {{OFS="\t"; print $1, $5, $3, $7, $8 }} else {{OFS="\t"; print $1, $2, $6, $7, $8 }}}}' | \
        awk -F'\t' -v ins_max={params.max_insert} -v ins_min={params.min_insert} ' ($3-$2)>=ins_min  && ($3-$2)<=ins_max' | \
        bedtools intersect -a stdin -b {params.blacklist} -v | \
        sort -S 2G -k1,1V -k2,2n | gzip -c > {output}
        '''

rule create_scp_bed:
    input: "{cell}.dedup.bam"
    output: "{cell}.scp.bed.gz"
    params:
        blacklist=config["ref"]["blacklist"],
        min_mapq=config["dna"]["min_mapq"],
        max_insert=config["dna"]["max_insert"],
        min_insert=config["dna"]["min_insert"],
        temp=config["path"]["temp"],
        dedup_script="workflow/scripts/fuzzydupsRADICAL.pl"
    shell:
        '''
        samtools sort -T {params.temp} -n -m 4G {input} | \
        samtools view -q {params.min_mapq} -b -f 0x2 -F 0x100 - | \
        bedtools bamtobed -bedpe -i stdin | \
        awk -F'\t' '$1~/^chr([0-9]+|[XY])$/ && $1==$4 {{ if($9=="-") {{OFS="\t"; print $1, $5, $3, $7, $8 }} else {{OFS="\t"; print $1, $2, $6, $7, $8 }}}}' | \
        awk -F'\t' -v ins_max={params.max_insert} -v ins_min={params.min_insert} ' ($3-$2)>=ins_min  && ($3-$2)<=ins_max' | \
        bedtools intersect -a stdin -b {params.blacklist} -v | \
        sort -S 2G -k1,1V -k2,2n | perl {params.dedup_script} | gzip -c > {output}
        ''' 

# Count fragments over bins
rule bincount_fixed:
    input:
        cell="{cell}.bed.gz",
        bins="resources/fixed-{binsize}.bed"
    output: "{cell}-bincounts-{binsize}.tsv"
    params:
        cellid=lambda wildcards: os.path.basename(wildcards.cell),
        genome=config["ref"]["chr_list_intersect"]
    shell:
        '''
        echo {params.cellid} > {output}
        bedtools intersect -g {params.genome} -sorted -F 10E-9 -c \
        -a {input.bins} -b {input.cell} | cut -f4 >> {output}
        '''

# Collect bincounts per patient
rule collect_bincounts_patient:
    input: 
        lambda wildcards: expand(
            os.path.join(dna_dir,"{cell}/{cell}-bincounts-{binsize}.tsv"), 
            cell=get_cells_dna(wildcards.patient_id), 
            binsize=wildcards.binsize
        )
    output: out + "/{patient_id}/{patient_id}-bincounts-{binsize}.tsv.gz"
    run:
        import gzip
        with gzip.open(output[0], 'wt') as writer:
            readers = [open(filename) for filename in input]
            for lines in zip(*readers):
                counts = [line.strip() for i,line in enumerate(lines)]
                print('\t'.join(counts), file=writer)


# Normal panel processing: use input base name to output file, as identifier
rule collect_bincounts_normals_bins:
    input: 
        lambda wildcards: expand(
            os.path.join(dna_dir,"{cell}/{cell}-bincounts-{binsize}.tsv"),
            cell=normals_bins["cell"], 
            binsize=wildcards.binsize
        )
    output: "resources/normals_bins-" + str(normals_bins_id) + "-{binsize}.tsv.gz"
    run:
        import gzip
        with gzip.open(output[0], 'wt') as writer:
            readers = [open(filename) for filename in input]
            for lines in zip(*readers):
                counts = [line.strip() for i,line in enumerate(lines)]
                print('\t'.join(counts), file=writer)

rule collect_bincounts_normals_scaling:
    input: 
        lambda wildcards: expand(
            os.path.join(dna_dir,"{cell}/{cell}-bincounts-{binsize}.tsv"),
            cell=normals_scaling["cell"], 
            binsize=wildcards.binsize
        )
    output: "resources/normals_scaling-" + str(normals_scaling_id) + "-{binsize}.tsv.gz"
    run:
        import gzip
        with gzip.open(output[0], 'wt') as writer:
            readers = [open(filename) for filename in input]
            for lines in zip(*readers):
                counts = [line.strip() for i,line in enumerate(lines)]
                print('\t'.join(counts), file=writer)

rule get_goodbins:
    input:
        counts="resources/normals_bins-" + normals_bins_id + "-{binsize}.tsv.gz",
        map="resources/fixed-{binsize}.map.txt",
        gc="resources/fixed-{binsize}.gc.txt",
        bins="resources/fixed-{binsize}.bed"
    output:
        pdf="resources/badbins-{binsize}.pdf",
        badbins="resources/badbins-{binsize}.bed",
        goodbins="resources/goodbins-{binsize}.bed"
    params:
        gc_min=config["dna"]["bin_min_gc"],
        map_min=config["dna"]["bin_min_map"],
        centromeres=config["ref"].get("centromeres", None),
        genome_gaps=config["ref"].get("gaps", None)
    script: "../scripts/call_badbins.R"
