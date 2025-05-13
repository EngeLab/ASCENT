rule rna_fusions_seedfile:
    input: 
        config["cells_rna"]
    output: 
        fusion_dir + "/{patient_id}.seedfile.txt"
    run:
        patient_df = cells_rna[cells_rna['patient_id'] == wildcards.patient_id].sort_index()
        with open(output[0], 'w') as f:
            for _, row in patient_df.iterrows():
                f.write(f"{row['cell']}\t{row['fq1']}\t{row['fq2']}\n")

rule rna_fusions:
    """
    Run STAR-Fusion as pseudobulk per patient.
    Samplefile format: [cell_id]<tab>[fq1]<tab>[fq2]
    Note! Creates massive temporary folder, which is automatically cleaned up by Snakemake's temp() directive.
    """
    input: 
        samplefile = fusion_dir + "/{patient_id}.seedfile.txt"
    output: 
        run = temp(directory(fusion_dir + "/{patient_id}/run")),
        predictions = fusion_dir + "/{patient_id}/star-fusion.fusion_predictions.tsv",        
        predictions_abridged = fusion_dir + "/{patient_id}/star-fusion.fusion_predictions.abridged.tsv",        
        deconvolved = fusion_dir + "/{patient_id}/star-fusion.fusion_predictions.tsv.samples_deconvolved.tsv",        
        deconvolved_abridged = fusion_dir + "/{patient_id}/star-fusion.fusion_predictions.tsv.samples_deconvolved.abridged.tsv"
    params:
        refdir = config["ref"]["starfusion_ref"],
        result_prefix = "star-fusion.fusion_predictions"
    conda: "../envs/starfusion.yaml"
    threads: 15
    log: fusion_dir + "/{patient_id}/star-fusion.log"
    shell:
        """
        STAR-Fusion --version
        STAR-Fusion --samples_file {input.samplefile} \
            --genome_lib_dir {params.refdir} \
            -O {output.run} \
            --CPU {threads}
        
        # Explicitly copy 
        for out in {output.predictions} {output.predictions_abridged} {output.deconvolved} {output.deconvolved_abridged}; do
            cp -v {output.run}/$(basename $out) $out
        done
        """

rule collect_rna_fusions:
    input: 
        expand(fusion_dir + "/{patient_id}/star-fusion.fusion_predictions.tsv",
               patient_id=dict.fromkeys(cells_rna['patient_id']) if is_rna_analysis() else [])
    output: 
        fusion_dir + "/rna_fusions.txt"
    shell: 
        "cat {input} > {output}"