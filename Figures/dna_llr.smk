# DNA segmentation by cell
#For downsampling and segmentation with llr (hypothesis free segmentation) exchange this rule for 
#Current dna_llr.smk in pipeline 
 
rule downsample_dna:
    input:
        cell="{cell}-bincounts-{binsize}.tsv",
        good_bins="resources/goodbins-{binsize}.bed"
    output:
        cell_out="{cell}-bincounts-{binsize}-{downsample}.tsv"
    params:
        downsample="{downsample}"
    script: "../scripts/Downsample.R"

rule segment_downsampled_dna:
    input:
        cell="{cell}-bincounts-{binsize}-{downsample}.tsv",
        map="resources/fixed-{binsize}.map.txt",
        gc="resources/fixed-{binsize}.gc.txt",
        bins="resources/fixed-{binsize}.bed",
        badbins="resources/badbins-{binsize}.bed" if config["dna"]["normalize_to_panel"] else []
    output:
        seg="{cell}-segments-{binsize}-{downsample}.txt",
        seg_cn="{cell}-segments_cn-{binsize}-{downsample}.txt",
        cn="{cell}-copynumber-{binsize}-{downsample}.txt",
        pdf_qc="{cell}-llr-qc-{binsize}-{downsample}.pdf",
        cnv_png="{cell}-cnv-{binsize}-{downsample}.png",
        cnv2_png="{cell}-cnv2-{binsize}-{downsample}.png",
        normcounts="{cell}-normcounts-{binsize}-{downsample}.txt",
        scalefactors="{cell}-scalefactors-{binsize}-{downsample}.txt"
    params:
        gc_min=config["dna"]["bin_min_gc"],
        map_min=config["dna"]["bin_min_map"],
        reads_per_window=config["dna"]["reads_per_window"] if "reads_per_window" in config["dna"].keys() else 300,
        max_w=10000,
        min_w=100,
        min_scale=config["dna"]["min_scale_factor"],
        max_scale=config["dna"]["max_scale_factor"],
        cn_trunc=config["dna"]["plot_max_cn"],
        src_general="workflow/scripts/general.R",
        src_cpp="workflow/scripts/calcLLR.cpp"
    script: "../scripts/call_breakpoints.R"


rule segment_dna:
    input:
        cell="{cell}-bincounts-{binsize}.tsv",
        map="resources/fixed-{binsize}.map.txt",
        gc="resources/fixed-{binsize}.gc.txt",
        bins="resources/fixed-{binsize}.bed",
        badbins="resources/badbins-{binsize}.bed" if config["dna"]["normalize_to_panel"] else []
    output:
        seg="{cell}-segments-{binsize}.txt",
        seg_cn="{cell}-segments_cn-{binsize}.txt",
        cn="{cell}-copynumber-{binsize}.txt",
        pdf_qc="{cell}-llr-qc-{binsize}.pdf",
        cnv_png="{cell}-cnv-{binsize}.png",
        cnv2_png="{cell}-cnv2-{binsize}.png",
        normcounts="{cell}-normcounts-{binsize}.txt",
        scalefactors="{cell}-scalefactors-{binsize}.txt"
    params:
        gc_min=config["dna"]["bin_min_gc"],
        map_min=config["dna"]["bin_min_map"],
        reads_per_window=config["dna"]["reads_per_window"] if "reads_per_window" in config["dna"].keys() else 300,
        max_w=10000,
        min_w=100,
        min_scale=config["dna"]["min_scale_factor"],
        max_scale=config["dna"]["max_scale_factor"],
        cn_trunc=config["dna"]["plot_max_cn"],
        src_general="workflow/scripts/general.R",
        src_cpp="workflow/scripts/calcLLR.cpp"
    script: "../scripts/call_breakpoints.R"

# Result aggregation
rule concat_copynumber_llr:
    input: 
        lambda wildcards: expand(
            os.path.join(dna_dir,"{cell}/{cell}-copynumber-{binsize}-{downsample}.txt"),
            cell=get_cells_dna(wildcards.patient_id), 
            binsize=wildcards.binsize,
            downsample=wildcards.downsample	
        )
    output: out + "/{patient_id}/{patient_id}-copynumber_llr-{binsize}-{downsample}.tsv.gz"
    run:
        import gzip
        with gzip.open(output[0], 'wt') as writer:
            readers = [open(filename) for filename in input]
            for lines in zip(*readers):
                counts = [line.strip() for i,line in enumerate(lines)]
                print('\t'.join(counts), file=writer)

rule concat_normcounts:
    input:
        lambda wildcards: expand(
            os.path.join(dna_dir,"{cell}/{cell}-normcounts-{binsize}.txt"),
            cell=get_cells_dna(wildcards.patient_id), 
            binsize=wildcards.binsize
        )
    output: out + "/{patient_id}/{patient_id}-normcounts-{binsize}.tsv.gz"
    run:
        import gzip
        cell_ids = get_cells_dna(wildcards.patient_id)
        with gzip.open(output[0], 'wt') as writer:
            print("\t".join(cell_ids), file=writer)
            readers = [open(filename) for filename in input]
            for lines in zip(*readers):
                counts = [line.strip() for i,line in enumerate(lines)]
                print('\t'.join(counts), file=writer)

# Rules for individual cell llr clustering
rule cluster_llr:
    input:
        copynumber=out + "/{patient_id}/{patient_id}-copynumber_llr-{binsize}.tsv.gz",
        counts=out + "/{patient_id}/{patient_id}-bincounts-{binsize}.tsv.gz",
        map="resources/fixed-{binsize}.map.txt",
        gc="resources/fixed-{binsize}.gc.txt",
        bins="resources/fixed-{binsize}.bed",
        badbins="resources/badbins-{binsize}.bed" if config["dna"]["normalize_to_panel"] else [],
        normal_bincounts="resources/normals_scaling-" + normals_scaling_id + "-{binsize}.tsv.gz" if config["dna"]["normalize_to_panel"] else [],
        meta_per_cell=out + "/{patient_id}/{patient_id}-metadata_long.tsv",
        dna_qc=out + "/{patient_id}/qc/{patient_id}-qc_dna.tsv",
        rna_phase=lambda wildcards: out + "/{patient_id}/{patient_id}-rna_phases.txt" if has_rna_data(wildcards.patient_id) and config["dna"]["exclude_sphase"] else []
    output:
        rda=out + "/{patient_id}/clones_llr/{patient_id}-clones_llr-{binsize}.Rda",
        clones=out + "/{patient_id}/clones_llr/{patient_id}-clones_llr-{binsize}.txt",
        qc_clones=out + "/{patient_id}/clones_llr/{patient_id}-qc_clones_llr-{binsize}.pdf",
        heatmap=out + "/{patient_id}/clones_llr/{patient_id}-heatmap_llr-{binsize}.png",
        heatmap_plate=out + "/{patient_id}/clones_llr/{patient_id}-heatmap_llr-byplate-{binsize}.png",
        heatmap_raw=out + "/{patient_id}/clones_llr/{patient_id}-heatmap_llr_all_cells-{binsize}.png",
        consensus_clone_profiles=out + "/{patient_id}/clones_llr/{patient_id}-clone_consensus_llr-{binsize}.pdf",
        segments=out + "/{patient_id}/clones_llr/{patient_id}-segments_llr-{binsize}.tsv"
    params:
        cells=lambda wildcards: get_cells_dna(wildcards.patient_id),
        normals=config["dna"]["normals_llr"],
        min_cells=config["dna"]["min_cells_per_cluster"],
        umap_n_components=30,
        umap_log=config["dna"].get("umap_log", True),
        umap_metric="manhattan",
        umap_neighbors=20,
        umap_mindist=0.01,
        umap_spread=1,
        max_clone_diff_frac=0.15,
        max_clone_excl_frac=0.3,
        tree_dist_max=1500,
        tree_method=config["dna"].get("tree_method", "nj"),
        gc_min=config["dna"]["bin_min_gc"],
        map_min=config["dna"]["bin_min_map"],
        counts_min=config["dna"]["min_count"],
        heatmap_smooth_factor=10,
        gamma=config["dna"].get("multipcf_gamma", 10),
        prune1=config["dna"].get("prune1", 5),
        prune2=config["dna"].get("prune2", 10),
        min_scale=config["dna"]["min_scale_factor"],
        max_scale=config["dna"]["max_scale_factor"],
        cytoband_file=config["ref"]["cytobands"],
        exclude_y=config["dna"].get("exclude_chrY", True),
        src_general="workflow/scripts/general.R"
    threads: 8
    script: "../scripts/cluster_cells_byPatient.R"

rule collect_clones_llr:
    input: 
        expand(out + "/{patient_id}/clones_llr/{patient_id}-clones_llr-{{binsize}}.txt",
               patient_id=get_patients_dna()) if is_dna_analysis() else []
    output: 
        out + "/clones_llr-{binsize}.txt"
    shell: "cat {input} > {output}"
