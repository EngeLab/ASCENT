def get_filtered_dna(wildcards):
    '''
    Returns list of cells from patient that pass QC thresholds
    rule qc_dna becomes a checkpoint, to force re-evaluation after running
    '''
    qc_file = checkpoints.qc_dna.get(**wildcards).output.qc_tsv
    qc_data = pd.read_table(qc_file)
    good_cells = qc_data[
        (qc_data['bam_read_pairs'] > config["dna"]["min_count"]) & 
        (qc_data['bam_dup_frac'] < config["dna"]["max_dup_frac"]) &
        (qc_data['dna_library_id'].isin(cells_dna['cell'][cells_dna.patient_id==wildcards.patient_id]))
    ]['dna_library_id']
    return list(good_cells)

rule merge_dna_bam:
    ''' 
    Make pseudobulk bam from filtered list of cells
    Note: takes dependency on all bam-files from patient and selects before merger
    '''
    input: 
        bams=lambda wildcards: expand(os.path.join(dna_dir,"{cell}/{cell}.dedup.bam"), cell=get_filtered_dna(wildcards)),
    output: 
        bam=out + "/{patient_id}/{patient_id}.pseudobulk.dna.bam",
        idx=out + "/{patient_id}/{patient_id}.pseudobulk.dna.bam.bai"
    params:
        genome_bed=config["ref"]["genome_bed"]
    threads: 12
    shell:
        '''
        samtools merge -@ {threads} -L {params.genome_bed} -o {output.bam} {input.bams}
        samtools index -o {output.idx} {output.bam}
        '''

rule call_germline_hets: 
    '''
    Use 1kg callset on all cells to find likely heterozygous germline variants
    MAF 0.01 means at least 1/100 reads
    AD >= 2 in practice means a SNP is detected in at least 2 different cells
    '''
    input: 
        bam=out + "/{patient_id}/{patient_id}.pseudobulk.dna.bam",
        idx=out + "/{patient_id}/{patient_id}.pseudobulk.dna.bam.bai",
        vcf_ref="/wrk/resources/genomes/1kGP_highcoverage_hg38/1kGP_highcoverage.{chr}.snv_filtered.vcf.gz"
    output: 
        cellsnp_dir=temp(directory(out + "/{patient_id}/{patient_id}_{chr}")),
        vcf=temp(out + "/{patient_id}/{patient_id}.het.{chr}.vcf.gz"),
        tbi=temp(out + "/{patient_id}/{patient_id}.het.{chr}.vcf.gz.tbi")
    params:
        fasta=config["ref"]["fasta"],
        min_AD=2
    threads: 4
    conda: "../envs/snv.yaml"
    shell: 
        '''
        cellsnp-lite \
        -s {input.bam} \
        -I {wildcards.patient_id} \
        -O {output.cellsnp_dir} \
        --cellTAG None \
        --genotype \
        --gzip \
        -p {threads} \
        -f {params.fasta} \
        -R {input.vcf_ref} \
        --UMItag None \
        --minMAF 0.01 \
        --minCOUNT 5 \
        --minMAPQ 20 \
        --maxDEPTH 100 \
        --countORPHAN \
        --exclFLAG UNMAP,SECONDARY,QCFAIL,DUP \
        --gzip

        bcftools filter -i 'GT="het" & INFO/AD[0] >= {params.min_AD}' {output.cellsnp_dir}/cellSNP.cells.vcf.gz | bgzip > {output.vcf}
        tabix {output.vcf}
        '''

rule phase_variants: 
    '''
    Phasing panel from HGDP+1kg. Info: https://gnomad.broadinstitute.org/downloads#v3-hgdp-1kg-tutorials
    Downloaded from gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2
    NOTE! Not using chrX currently. Renamed the non-PAR chrX, available
    '''
    input: 
        vcf=out + "/{patient_id}/{patient_id}.het.{chr}.vcf.gz",
        idx=out + "/{patient_id}/{patient_id}.het.{chr}.vcf.gz.tbi",
        phase_ref="/wrk/resources/genomes/hgdp_1kg/phased_haplotypes_v2/hgdp1kgp_{chr}.filtered.SNV_INDEL.phased.shapeit5.bcf",
    output: temp(out + "/{patient_id}/{patient_id}.phased.{chr}.vcf.gz")
    params: 
        genetic_map="/wrk/resources/genomes/hgdp_1kg/genetic_map_hg38_withX.txt.gz",
        outprefix=out + "/{patient_id}/{patient_id}.phased.{chr}"
    conda: "../envs/snv.yaml"
    shell:
        '''
    eagle \
    --vcfTarget={input.vcf} \
    --vcfRef={input.phase_ref} \
    --geneticMapFile={params.genetic_map} \
    --vcfOutFormat=z \
    --outPrefix={params.outprefix}
        '''

rule collect_phased_vcfs:
    input: expand(out + "/{{patient_id}}/{{patient_id}}.phased.{chr}.vcf.gz", chr=['chr' + str(i) for i in range(1,23)])
    output: out + "/{patient_id}/{patient_id}.phased.vcf.gz"
    conda: "../envs/snv.yaml"
    shell: 
        ''' 
        bcftools concat -o {output} {input}
        '''

def get_patient_vcf(wildcards):
    cells_with_patients = cells_dna[cells_dna['patient_id'].notna()]
    # Get (1) patient_id for this cell
    patient_id = cells_with_patients.loc[cells_with_patients['cell'] == wildcards.cell, 'patient_id'].iloc[0]
    return f"{out}/{patient_id}/{patient_id}.phased.vcf.gz"

rule genotype_cells:
    """
    Get allele-specific read counts at phased heterozygous sites
    1. Only allow alt genotypes that match phased VCF (input.vcf)
    """
    input:
        vcf=get_patient_vcf,
        bam=dna_dir + "/{cell}/{cell}.dedup.bam"
    output:
        baf=dna_dir + "/{cell}/{cell}.baf.txt"
    params:
        fasta=config["ref"]["fasta"],
        min_mapq=20,
        min_baseq=1,
    run:
        import subprocess
        
        # Get phase information
        phase_info = {}
        gt_info = {}
        vcf_ref = {}
        vcf_alt = {}
        cmd = f"bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' {input.vcf}"
        proc = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        
        for line in proc.stdout.splitlines():
            chrom, pos, ref, alt, gt = line.strip().split('\t')
            phase_info[f"{chrom}_{pos}"] = gt.split('|')[0] == '1'
            gt_info[f"{chrom}_{pos}"] = gt
            vcf_ref[f"{chrom}_{pos}"] = ref
            vcf_alt[f"{chrom}_{pos}"] = alt
        
        print(f"\n[{wildcards.cell}] Loaded {len(phase_info)} phased positions")

        # Run mpileup to get read counts 
        cmd = (f"bcftools mpileup -d 50 -f {params.fasta} -a INFO/AD,FORMAT/DP,FORMAT/SP,FORMAT/QS "
               f"-q {params.min_mapq} -Q {params.min_baseq} "
               f"-T {input.vcf} {input.bam} | "
               "bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%AD{0}\t%AD{1}\t[%QS{0}]\t[%QS{1}]\n'")

        #print(f"\nRunning mpileup command: {cmd}")
        proc = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)

        n_written = 0
        n_mismatched = 0
        with open(output.baf, 'w') as outf:
            # header = ['chrom', 'pos', 'a_count', 'b_count',
            #           'ref', 'alt', 'vcf_ref', 'vcf_alt', 
            #           'ref_count', 'alt_count', 'ref_quality', 'alt_quality', 'original_genotype']
            # outf.write('\t'.join(header) + '\n')

            for line in proc.stdout.splitlines():
                chrom, pos, ref, alt, ref_count, alt_count, ref_q, alt_q = line.strip().split('\t')

                # Check for alt mismatch ALT allele and VCF ALT
                if alt != "<*>" and alt != vcf_alt[f"{chrom}_{pos}"]:
                    n_mismatched += 1
                    continue
                
                # Flip counts if phase is 1|0
                is_flipped = phase_info[f"{chrom}_{pos}"]
                o_gt = gt_info[f"{chrom}_{pos}"]
                if is_flipped:
                    a_count, b_count = alt_count, ref_count
                else:
                    a_count, b_count = ref_count, alt_count
                
                # fields = [chrom, pos, a_count, b_count]
                fields = [chrom, pos, wildcards.cell, a_count, b_count,
                          ref, alt, vcf_ref[f"{chrom}_{pos}"], vcf_alt[f"{chrom}_{pos}"], 
                          ref_count, alt_count, ref_q, alt_q, o_gt]

                outf.write('\t'.join(map(str, fields)) + '\n')
                n_written += 1

        print(f"[{wildcards.cell}] Excluded {n_mismatched} mismatched ALT positions")
        print(f"[{wildcards.cell}] Wrote {n_written} positions to {output.baf}")
        

rule collect_genotyped_cells:
    """Concatenate A/B counts of cells that pass QC and sort"""
    input: lambda wildcards: expand(dna_dir + "/{cell}/{cell}.baf.txt", cell=get_filtered_dna(wildcards))
    output: out + "/{patient_id}/{patient_id}-baf.txt.gz"
    threads: 4  
    shell:
        '''
        cat {input} | sort -k1,1V -k2,2n --parallel={threads} | gzip > {output}
        '''
