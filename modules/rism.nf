#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process solvation {
    container "${params.container__biobb_amber}"
    publishDir "${params.output_folder}/${params.database}/3DRISM/${molecule}_${partial_charge_method}_${model}_${temperature}", mode: 'copy', overwrite: false

    debug false
    input:
    tuple val(molecule), path(prm), path(crd), val(partial_charge_method), val(model), val(temperature), path(xvv), path(pdb), path(rst)
    output:
    path("*"), emit: all
    path("rism.log"), emit: log
    path("${molecule}_${partial_charge_method}_${model}_${temperature}.json"), emit: json
    maxRetries 20

    script:
    """
    #!/usr/bin/env python
    print("Hello from ${molecule} ${prm} ${crd} ${model} ${temperature} ${xvv} ${pdb} ${rst}")
    def setup_calculation(output_file_prefix="3DRISM",
                        closure='pse3', write_g=False, write_h=False,
                        write_c=False,
                        write_u=False, write_asymp=False,
                        noasympcorr=False,
                        buffer_distance=25.0,
                        solvbox=False,
                        grdspc=(0.5, 0.5, 0.5),
                        tolerance=[1e-5], polar_decomp=False,
                        verbose=2, maxstep=500,
                        rism3d_path='rism3d.snglpnt'):

        grdspc = ','.join(map(str, grdspc))
        if solvbox:
            solvbox = ','.join(map(str, solvbox))
        run_flags_list = [rism3d_path,
            '--pdb', "${pdb}",
            '--prmtop', "${prm}",
            '--rst', "${rst}",
            '--xvv', "${xvv}",
            '--grdspc', grdspc,]
        run_flags_list.extend(['--tolerance'] + tolerance)
        run_flags_list.extend(['--closure', closure])
        if solvbox:
            run_flags_list.extend(['--solvbox', solvbox])
        else:
            run_flags_list.extend(['--buffer', str(buffer_distance)])
        if write_g:
            run_flags_list.extend(['--guv',
                                        'g_{}'.format(output_file_prefix)])
        if write_h:
            run_flags_list.extend(['--huv',
                                        'h_{}'.format(output_file_prefix)])
        if write_c:
            run_flags_list.extend(['--cuv',
                                        'c_{}'.format(output_file_prefix)])
        if write_u:
            run_flags_list.extend(['--uuv',
                                        'u_{}'.format(output_file_prefix)])
        if write_asymp:
            run_flags_list.extend(['--asymp',
                                        'a_{}'.format(output_file_prefix)])
        if noasympcorr:
            run_flags_list.extend(['--noasympcorr'])
        if polar_decomp:
            run_flags_list.extend(['--polarDecomp'])
        if verbose:
            run_flags_list.extend(['--verbose'])
            run_flags_list.extend(['{}'.format(verbose)])
        if maxstep:
            run_flags_list.extend(['--maxstep'])
            run_flags_list.extend(['{}'.format(maxstep)])
        run_flags_list.extend(['--pc+'])
        return run_flags_list

    run_flags_list = setup_calculation()

    #print(run_flags_list)
    # Convert non-string elements to strings and create a space-separated string
    space_separated_string = ' '.join(map(str, run_flags_list))
    #print(space_separated_string)
    import subprocess
    extracted_number = 0.0
    foundNumber = False
    # Command to run
    # Run the command and capture the output
    #output = subprocess.check_output(space_separated_string, shell=True, text=True)
    process = subprocess.run(['bash', '-c', space_separated_string], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    lines = process.stdout.splitlines()
    for line in lines:
        # Split the line into parts based on whitespace
        parts = line.split()
        # Check if the line starts with the desired prefix
        if len(parts) >= 2 and parts[0] == "rism_excessChemicalPotentialPCPLUS":
            # Extract the number from the second part
            extracted_number = float(parts[1])
            print(extracted_number)
            foundNumber = True

    if foundNumber:
        # Extract the value from the match
        print("PC+dG*(solv)",extracted_number,"kcal/mol")
    else:
        print("ERROR!")

    # Open the file in write mode and write the string
    with open('rism.log', 'w') as file:
        file.write(process.stdout)

    import json
    results = {"molecule":"${molecule}", "partial_charge_method":"${partial_charge_method}", "solvent" : "${model}", "temperature" : ${temperature}, "PC+dG*(solv)(kcal/mol)":extracted_number}
    with open("${molecule}_${partial_charge_method}_${model}_${temperature}.json", 'w') as file:
        json.dump(results, file)
    """
}

process write_solvation_script {
    container "${params.container__biobb_amber}"

    debug false
    input:
    tuple val(molecule), path(prm), path(crd), val(partial_charge_method), val(model), val(temperature), path(xvv), path(pdb), path(rst)
    output:
    tuple val(molecule), path(prm), path(crd), val(partial_charge_method), val(model), val(temperature), path(xvv), path(pdb), path(rst), emit: system
    path("rism.sh"), emit: bash_script

    script:
    """
    #!/usr/bin/env python
    print("Hello from ${molecule} ${prm} ${crd} ${model} ${temperature} ${xvv} ${pdb} ${rst}")
    def setup_calculation(output_file_prefix="3DRISM",
                        closure='pse3', write_g=False, write_h=False,
                        write_c=False,
                        write_u=False, write_asymp=False,
                        noasympcorr=False,
                        buffer_distance=25.0,
                        solvbox=False,
                        grdspc=(0.5, 0.5, 0.5),
                        tolerance=[1e-5], polar_decomp=False,
                        verbose=2, maxstep=500,
                        rism3d_path='rism3d.snglpnt'):

        grdspc = ','.join(map(str, grdspc))
        if solvbox:
            solvbox = ','.join(map(str, solvbox))
        run_flags_list = [rism3d_path,
            '--pdb', "${pdb}",
            '--prmtop', "${prm}",
            '--rst', "${rst}",
            '--xvv', "${xvv}",
            '--grdspc', grdspc,]
        run_flags_list.extend(['--tolerance'] + tolerance)
        run_flags_list.extend(['--closure', closure])
        if solvbox:
            run_flags_list.extend(['--solvbox', solvbox])
        else:
            run_flags_list.extend(['--buffer', str(buffer_distance)])
        if write_g:
            run_flags_list.extend(['--guv',
                                        'g_{}'.format(output_file_prefix)])
        if write_h:
            run_flags_list.extend(['--huv',
                                        'h_{}'.format(output_file_prefix)])
        if write_c:
            run_flags_list.extend(['--cuv',
                                        'c_{}'.format(output_file_prefix)])
        if write_u:
            run_flags_list.extend(['--uuv',
                                        'u_{}'.format(output_file_prefix)])
        if write_asymp:
            run_flags_list.extend(['--asymp',
                                        'a_{}'.format(output_file_prefix)])
        if noasympcorr:
            run_flags_list.extend(['--noasympcorr'])
        if polar_decomp:
            run_flags_list.extend(['--polarDecomp'])
        if verbose:
            run_flags_list.extend(['--verbose'])
            run_flags_list.extend(['{}'.format(verbose)])
        if maxstep:
            run_flags_list.extend(['--maxstep'])
            run_flags_list.extend(['{}'.format(maxstep)])
        run_flags_list.extend(['--pc+'])
        return run_flags_list

    run_flags_list = setup_calculation()

    #print(run_flags_list)
    # Convert non-string elements to strings and create a space-separated string
    space_separated_string = ' '.join(map(str, run_flags_list))

    RISM_SCRPT = '''#!/bin/bash
    {run_command} > rism.out
    '''
    my_string=RISM_SCRPT.format(run_command=space_separated_string)
    # Open the file in write mode ('w' or 'wt' for text mode)
    with open("rism.sh", 'w') as file:
        # Write the string to the file
        file.write(my_string)
    """
}


process run_solvation_script {
    container "${params.container__biobb_amber}"
    publishDir "${params.output_folder}/${params.database}/3DRISM/${molecule}_${partial_charge_method}_${model}_${temperature}", mode: 'copy', overwrite: false

    debug false
    errorStrategy 'ignore'
    input:
    tuple val(molecule), path(prm), path(crd), val(partial_charge_method), val(model), val(temperature), path(xvv), path(pdb), path(rst)
    path(bash_script)
    output:
    tuple val(molecule), val(partial_charge_method), val(model), val(temperature), env(chemical_potential), emit: system
    tuple path(bash_script), path("rism.out"), emit: files
    shell:
    """
    # Turn off exit-on-error temporarily
    set +e
    source !{bash_script}
    # Revert to the default behavior (exit on error)
    set -e
    cat rism.out
    chemical_potential=\$(awk '/^rism_excessChemicalPotentialPCPLUS/ {print \$2}' rism.out)
    if [ -z "\$chemical_potential" ]; then
        chemical_potential="NA"
    fi
    """
}

workflow rism_solvation {
    take:
    minimized_systems
    main:
    write_solvation_script(minimized_systems) |
    run_solvation_script

    emit:
    system = run_solvation_script.out.system
}