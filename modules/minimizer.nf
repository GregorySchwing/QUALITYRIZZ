#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process minimize_ligand {
    container "${params.container__biobb_amber}"
    publishDir "${params.output_folder}/${params.database}/minimizations/${molecule}_${model}_${temperature}", mode: 'copy', overwrite: false

    debug true
    input:
    tuple val(molecule), path(prm), path(crd), val(model), val(temperature), path(xvv) 
    output:
    path("sander.*"), emit: paths

    script:
    """
    #!/usr/bin/env python
    print("Hello from ${molecule} ${prm} ${crd} ${model} ${temperature} ${xvv}")
    # Import module
    from biobb_amber.sander.sander_mdrun import sander_mdrun

    # Create prop dict and inputs/outputs
    output_n_min_traj_file = 'sander.n_min.x'
    output_n_min_rst_file = 'sander.n_min.rst7'
    output_n_min_log_file = 'sander.n_min.log'

    # Minimization script from PC_PLUS
    prop = {
        'simulation_type' : "minimization",
        "mdin" : { 
            'imin' : 1, # perform minimization
            'maxcyc' : 500, # The maximum number of cycles of minimization
            'drms' : 1e-3, # RMS force
            'ntmin' : 3, # xmin algorithm
            'ntb' : 0, # no periodic boundary
            'cut' : 999, # non-bonded cutoff
            'ntpr' : 5, # printing frequency
            'ntxo' : 1, # asci formatted rst7
        },
        "dev" : "-xvv {xvv_file}".format(xvv_file="${xvv}"),
    }

    # Create and launch bb
    sander_mdrun(input_top_path="${prm}",
                input_crd_path="${crd}",
                output_traj_path=output_n_min_traj_file,
                output_rst_path=output_n_min_rst_file,
                output_log_path=output_n_min_log_file,
                properties=prop)

    """
}

workflow minimize_ligands {
    take:
    systems
    solvents
    main:
    // Process each JSON file asynchronously
    all=systems.combine(solvents)
    minimize_ligand(all)
}