rule segment_joint: 
    input:
        counts=out + "/{patient_id}/{patient_id}-bincounts-{binsize}.tsv.gz",
        goodbins="resources/goodbins-{binsize}.bed",
        map="resources/fixed-{binsize}.map.txt",
        gc="resources/fixed-{binsize}.gc.txt",
        normal_bincounts="resources/normals_scaling-" + str(normals_scaling_id) + "-{binsize}.tsv.gz",
        rna_phase=lambda wildcards: out + "/{patient_id}/{patient_id}-rna_phases.txt" if has_rna_data(wildcards.patient_id) and config["dna"]["exclude_sphase"] else [],
    output: 
        mpcf=out + "/{patient_id}/clones/{patient_id}-mpcf-g{gamma}-{binsize}.txt.gz",
        cn=out + "/{patient_id}/clones/{patient_id}-copynumber-g{gamma}-{binsize}.txt.gz",
        scaleplot=out + "/{patient_id}/clones/{patient_id}-scaleplot-g{gamma}-{binsize}.pdf",
        scalefactor=out + "/{patient_id}/clones/{patient_id}-scalefactors-g{gamma}-{binsize}.txt",
        sc_segs=out + "/{patient_id}/clones/{patient_id}-sc_segments-g{gamma}-{binsize}.txt",
    params: 
        cells=lambda wildcards: get_cells_dna(wildcards.patient_id),
        gamma=lambda wildcards: float(wildcards.gamma),
        counts_min=config["dna"]["min_count"],
        mpcf_exact=True,
        min_scale=config["dna"]["min_scale_factor"],
        max_scale=config["dna"]["max_scale_factor"],
        cytoband_file=config["ref"]["cytobands"],
    threads: 12
    script: "../scripts/segment_cells_multipcf.R"

rule cnv_calc_scp_scCN:
    '''
    Calculate overlap statistics based on median copynumber per clone
    '''
    input:
        frag_files=lambda wildcards: expand(
            os.path.join(dna_dir,"{cell}/{cell}.scp.bed.gz"),
            cell=get_cells_dna(wildcards.patient_id)
        ),
        cn=out + "/{patient_id}/clones/{patient_id}-copynumber-g{gamma}-{binsize}.txt.gz",
        sc_segments=out + "/{patient_id}/clones/{patient_id}-sc_segments-g{gamma}-{binsize}.txt",
        bins="resources/goodbins-{binsize}.bed",
        all_bins="resources/fixed-{binsize}.bed"
    output:
        scps=out + "/{patient_id}/clones/{patient_id}-scps-scCN-{binsize}-g{gamma}.txt"
    params:
        cells=lambda wildcards: get_cells_dna(wildcards.patient_id),
        min_binsize=100
    threads: 12
    script: "../scripts/scp_calc_cnv_scCN.R"

rule log_odds_scp:
    '''
    Log odds and plotting of copy number changes
    '''
    input:
        scps_sc=out + "/{patient_id}/clones/{patient_id}-scps-scCN-{binsize}-g{gamma}.txt",
        bins="resources/goodbins-{binsize}.bed"
    output:
        log_odds_sc=out + "/{patient_id}/clones/{patient_id}-log_odds-scCN-{binsize}-g{gamma}.png",
        log_odds_df_sc=out + "/{patient_id}/clones/{patient_id}-log_odds_df-scCN-{binsize}-g{gamma}.tsv"
    script: "../scripts/scp_log_odds.R"


rule find_clones:
    '''
    Note: dna_qc and rna_phase are still cohort-dependent now
    '''
    input:
        mpcf=out + "/{patient_id}/clones/{patient_id}-mpcf-g{gamma}-{binsize}.txt.gz",
        cn=out + "/{patient_id}/clones/{patient_id}-copynumber-g{gamma}-{binsize}.txt.gz",
        bins="resources/goodbins-{binsize}.bed",
        all_bins="resources/fixed-{binsize}.bed",
        dna_qc=out + "/{patient_id}/qc/{patient_id}-qc_dna.tsv",
        metadata_patient=out + "/{patient_id}/{patient_id}-metadata_long.tsv",
        logodds=out + "/{patient_id}/clones/{patient_id}-log_odds_df-scCN-{binsize}-g{gamma}.tsv",
        rna_phase=lambda wildcards: out + "/{patient_id}/{patient_id}-rna_phases.txt" if has_rna_data(wildcards.patient_id) and config["dna"]["exclude_sphase"] else []
    output:
        clones=out + "/{patient_id}/clones/{patient_id}-clones-{method}-g{gamma}-{binsize}.txt",
        heatmap=out + "/{patient_id}/clones/{patient_id}-heatmap-{method}-g{gamma}-{binsize}.png",
        heatmap_revised=out + "/{patient_id}/clones/{patient_id}-heatmap_revised-{method}-g{gamma}-{binsize}.png",
        qc_corplot=out + "/{patient_id}/clones/{patient_id}-qc_cor-{method}-g{gamma}-{binsize}.pdf",
        qc_kmplot=out + "/{patient_id}/clones/{patient_id}-qc_km-{method}-g{gamma}-{binsize}.pdf",
        qc_dimred=out + "/{patient_id}/clones/{patient_id}-qc_dimred-{method}-g{gamma}-{binsize}.pdf"
    params:
        kmeans_threshold=lambda wildcards: get_patient_param(wildcards.patient_id, 'kmeans_threshold'),
        min_cells=lambda wildcards: get_patient_param(wildcards.patient_id, 'min_cells_per_cluster'),
        bin_filter=lambda wildcards: get_patient_param(wildcards.patient_id, 'min_bins_per_segment'),
        kmeans_max=config["dna"]["kmeans_max"],
        umap_metric="manhattan",
        umap_unique_threshold=35,
        umap_cluster_dims=30, 
        umap_neighbors=20, 
        umap_min_dist=0.1, 
        umap_spread=1.2, 
        umap_seed=42,
        src_general="workflow/scripts/general.R",
        bin_diff_max=0.2,
        bin_diff_zlim=3,
        clone_frac_diff=0.5
    threads: 6
    script: "../scripts/evaluate_clones.R"

rule collect_clone_stats:
    """Summarize clone statistics per patient across method, gamma, binsize."""
    input:
        lambda wildcards: [
            out + f"/{wildcards.patient_id}/clones/{wildcards.patient_id}-clones-{method}-g{g}-{binsize}.txt" 
            for g in get_patient_param(wildcards.patient_id, 'gamma')
            for method in config["dna"]["dim_reduction"]
            for binsize in config["dna"]["binsize"]
        ]
    output:
        out + "/{patient_id}/{patient_id}-clone_stats.txt"
    run:
        # Collect stats for each gamma and method
        results = []
        for clone_file in input:
            # Extract params from filename
            parts = clone_file.strip('.txt').split('-')
            method = parts[-3]
            gamma = parts[-2].lstrip('g')
            binsize = parts[-1]
            
            df = pd.read_csv(clone_file, sep='\t')
            clone_counts = df['clone_final'].value_counts()
            total_cells = len(df)
            unassigned = clone_counts.get('0', 0)
            n_clones = len(clone_counts[clone_counts.index != '0'])
            
            results.append({
                'binsize': binsize,
                'gamma': gamma,
                'method': method,
                'n_cells': total_cells,
                'n_assigned': total_cells - unassigned,
                'n_clones': n_clones
            })
            
        # Convert to table and save
        summary_df = pd.DataFrame(results)
        summary_df = summary_df.sort_values(['method', 'gamma'])
        summary_df.to_csv(output[0], sep='\t', index=False)


rule refine_clones_automatic:
    '''
    Refine clone automatically with parameters set in config file
    '''
    input:
        clones=out + "/{patient_id}/clones/{patient_id}-clones-umap-g{gamma}-{binsize}.txt",
        counts=out +"/{patient_id}/{patient_id}-bincounts-{binsize_refine}.tsv.gz",
        normal_cells="resources/normals_scaling-" + str(normals_scaling_id) + "-{binsize_refine}.tsv.gz",
        sf=out + "/{patient_id}/clones/{patient_id}-scalefactors-g{gamma}-{binsize}.txt",
        logodds=out + "/{patient_id}/clones/{patient_id}-log_odds_df-scCN-{binsize}-g{gamma}.tsv",
        bins="resources/fixed-{binsize_refine}.bed",
        map="resources/fixed-{binsize_refine}.map.txt",
        gc="resources/fixed-{binsize_refine}.gc.txt",
        cytoband="/wrk/resources/genomes/hg38-iGenome-ERCC/cytoBand.txt.gz",
        good_bins="resources/goodbins-{binsize_refine}.bed",
        meta=out + "/{patient_id}/{patient_id}-metadata_long.tsv",
        qc_dna=out + "/{patient_id}/qc/{patient_id}-qc_dna.tsv"
    params:
        clone_gamma=0.5
    output:
        chr_heatmap=out+ "/{patient_id}/clones/{patient_id}-final-clones-refined-g{gamma}-b{binsize}-br{binsize_refine}.pdf",
        sc_heatmap=out + "/{patient_id}/clones/{patient_id}-final-refined-clones-heatmap-g{gamma}-b{binsize}-br{binsize_refine}.png",
        final_clones=out + "/{patient_id}/clones/{patient_id}-final-refined-clones-g{gamma}-b{binsize}-br{binsize_refine}.txt",
        final_clone_object=out + "/{patient_id}/clones/{patient_id}-final_clone_object-g{gamma}-b{binsize}-br{binsize_refine}.Rds"
    script:
        "../scripts/refine_clones_automatic.R"


rule haplotype_refined_clones:
    '''
    Haplotype refined clones
    '''
    input:
       cn_obj=out + "/{patient_id}/clones/{patient_id}-final_clone_object-g{gamma}-b{binsize}-br{binsize_refine}.Rds",
       phasing=out + "/{patient_id}/{patient_id}-baf.txt.gz"
    output:
       heatmap=out + "/{patient_id}/clones/{patient_id}-allele_specific-heatmap-g{gamma}-b{binsize}-br{binsize_refine}.pdf",
       cmBAF=out + "/{patient_id}/clones/{patient_id}-cmBAFs-g{gamma}-b{binsize}-br{binsize_refine}.pdf",
       medicc=out + "/{patient_id}/clones/{patient_id}-medicc_input-g{gamma}-b{binsize}-br{binsize_refine}.txt",
       cn_obj_out=out + "/{patient_id}/clones/{patient_id}-final_clone_object-haplotyped-g{gamma}-b{binsize}-br{binsize_refine}.Rds"
    script:
       "../scripts/haplotype_clones.R"


