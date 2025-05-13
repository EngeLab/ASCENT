rule run_medicc:
    '''
    Run medicc 
    It is good to also run medicc outside the pipeline and give -n "diploid clone" so that it roots the tree in the diploid cells
    '''
    input:
        medicc_input=out + "/{patient_id}/clones/{patient_id}-medicc_input-g{gamma}-{binsize}.txt"
    output:
        medicc_dir = directory(out + "/{patient_id}/medicc/g{gamma}-{binsize}"),
        marker = out + "/{patient_id}/medicc/g{gamma}-{binsize}/.medicc_done"
    conda: "../envs/medicc_env.yaml"
    shell:
        '''
        medicc2 {input.medicc_input} {output.medicc_dir}
        touch {output.marker}
        '''



