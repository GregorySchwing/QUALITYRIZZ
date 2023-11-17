#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process build_ligand {
    container "${params.container__openff_toolkit}"
    publishDir "${params.output_folder}/${params.database}/systems/${pathToJson.baseName}", mode: 'copy', overwrite: true

    debug false
    input:
    path pathToJson
    output:
    val(pathToJson.baseName), emit: molecule
    path("ligand.prmtop"), emit: prm
    path("ligand.inpcrd"), emit: crd
    tuple val(pathToJson.baseName), path("ligand.prmtop"), path("ligand.inpcrd"), emit: system
    script:
    """
    #!/usr/bin/env python
    import json
    # Open and read the JSON file
    with open("${pathToJson}", "r") as file:
        print("reading ${pathToJson}")
        data = json.load(file)
    smiles=data["smiles"]
    # Imports from the toolkit
    from openff.toolkit import ForceField, Molecule, Topology
    from openff.units import Quantity, unit
    from openff.interchange import Interchange

    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem import rdmolfiles
    def embed(mol, seed=None):
        params = AllChem.ETKDGv2()
        if seed is not None:
            params.randomSeed = seed
        else:
            params.randomSeed = 123
        AllChem.EmbedMolecule(mol, params)
        return mol

    forcefield = ForceField("openff-2.1.0.offxml")
    openff_mol = Molecule.from_smiles(smiles,allow_undefined_stereo=True)
    rdmol = openff_mol.to_rdkit()
    rdmol3D = embed(rdmol,123)
    # Create a PDB file
    writer = rdmolfiles.MolToPDBFile(rdmol3D, "ligand.pdb")
    # Load the topology from a PDB file and `Molecule` objects
    topology = Topology.from_pdb(
        "ligand.pdb",
        unique_molecules=[openff_mol],
    )

    interchange = Interchange.from_smirnoff(
        force_field=forcefield,
        topology=topology,
    )

    interchange.to_prmtop("ligand.prmtop")
    interchange.to_inpcrd("ligand.inpcrd")
    
    """
}

process build_solvent {
    container "${params.container__biobb_amber}"
    publishDir "${params.output_folder}/${params.database}/1DRISM/${model}_${T}", mode: 'copy', overwrite: true

    debug true
    input:
    tuple val(model), val(T), path(mdl)
    output:
    path("${model}_${T}*.*"), emit: paths
    val(model), emit: model
    val(T), emit: temperature
    tuple val(model), val(T), path("${model}_${T}.xvv"), emit: solvent optional true
    shell:
    """
    #!/usr/bin/env python
    import json
    # Open and read the JSON file
    def water_dielectric_const(T):
        if not 253.15 <= T <= 383.15:
            raise ValueError("Temperature is outside of allowed range.")
        T_star = T/300.0
        coefs = [-43.7527, 299.504, -399.364, 221.327]
        exp_f = [-0.05, -1.47, -2.11, -2.31]
        e = 0
        for i in range(4):
            e += coefs[i]*T_star**(exp_f[i])
        return e


    def water_concentration(T):
        if not 253.15 <= T <= 383.15:
            raise ValueError("Temperature is outside of allowed range.")
        p0 = 10.0**5    # Pa
        R = 8.31464     # J/mol/K
        Tr = 10.0
        Ta = 593.0
        Tb = 232.0
        a = [1.93763157E-2,
            6.74458446E+3,
            -2.22521604E+5,
            1.00231247E+8,
            -1.63552118E+9,
            8.32299658E+9]
        b = [5.78545292E-3,
            -1.53195665E-2,
            3.11337859E-2,
            -4.23546241E-2,
            3.38713507E-2,
            -1.19946761E-2]
        n = [None, 4., 5., 7., 8., 9.]
        m = [1., 2., 3., 4., 5., 6.]
        def alpha(T):
            return Tr/(Ta - T)
        def beta(T):
            return Tr/(T - Tb)
        coef = a[0] + b[0]*beta(T)**m[0]
        for i in range(1, 6):
            coef += a[i]*alpha(T)**n[i] + b[i]*beta(T)**m[i]
        v0 = R*Tr/p0*coef  # m3/mol
        return 1/(v0*1000)    # mol/L
    dieps=round(water_dielectric_const($T), 3)
    conc = round(water_concentration($T), 3)

    df = {"model":"${model}", "T":$T, "dieps":dieps, "conc":conc}

    with open("solvParams.json", 'w') as file:
        json.dump(df, file)

    SOLV_SUCEPT_SCRPT = '''#!/bin/bash
    cat > ${model}_${T}_{closure}_{iter}.inp <<EOF
    &PARAMETERS
        THEORY='{rism1d}', CLOSURE='{closure}',           !Theory
        selftest=1, ! Verify inputs/outputs
        NR=16384, DR=0.025,                    !Grid
        OUTLIST='xCGT', rout=0,                !Output
        mdiis_nvec=20, mdiis_del=0.3, tolerance=1.e-12,       !MDIIS
        KSAVE={NUMSTEPS}, KOUT=1, maxstep={NUMSTEPS},       !Check pointing and iterations
        SMEAR=1, ADBCOR=0.5,                   !Electrostatics
        TEMPERATURE={temp}, DIEps={diel},      !bulk solvent properties
        NSP=1,
        entropicDecomp=0                       ! While debugging convergence
    /
        &SPECIES                               !SPC water
        DENSITY={conc}d0,
        MODEL="${mdl}"
    /
    EOF
    {renameCommand}
    rism1d "${model}_${T}_{closure}_{iter}"  > "${model}_${T}_{closure}_{iter}.out"
    '''
    
    T=${T}
    smodel="${model}"
    rism1d="DRISM"
    closure="PSE3"
    closure="KH"
    rism1d_name = '{smodel}_{temp}'.format(smodel=smodel,temp=T)
    diel = dieps
    conc = conc

    closure_list = ["KH", "PSE1", "PSE2", "PSE3"]
    RENAMECOMMAND = ""
    LASTRUN = ""
    for closure in closure_list:
        for iter in range(6):
            if (LASTRUN != ""):
                THISRUN = "${model}_${T}_{closure}_{iter}".format(closure=closure,iter=iter)
                RENAMECOMMAND = "mv {LASTRUN}.sav {THISRUN}.sav".format(LASTRUN=LASTRUN,THISRUN=THISRUN)
            input_string=SOLV_SUCEPT_SCRPT.format(temp=T, diel=diel, conc=conc,\
                                smodel=smodel, rism1d=rism1d,\
                                closure=closure.upper(),\
                                NUMSTEPS=10000,\
                                iter=iter,\
                                renameCommand=RENAMECOMMAND,\
                                name1d=rism1d_name)

            import subprocess
            # Run the script in the current Bash environment
            process = subprocess.run(['bash', '-c', input_string], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

            # Print the output
            print("Output:", process.stdout)

            # Print any errors
            if process.stderr:
                print("Errors:", process.stderr)

            LASTRUN = "${model}_${T}_{closure}_{iter}".format(closure=closure,iter=iter)
            if (process.stderr == ""):
                break
    
    """

}


process build_water_model {
    container "${params.container__biobb_amber}"
    publishDir "${params.output_folder}/${params.database}/solvents/c${model}", mode: 'copy', overwrite: true

    debug false
    input:
    tuple val(model), val(temperature)
    output:
    tuple val(model), val(temperature), path("c*.mdl")

    shell:
    """
    #!/usr/bin/env python
    import subprocess
    import parmed
    import numpy as np

    # This script uses tleap to build the single residue 
    # CRD/TOP, but it could be extended to take these
    # As arguments.  

    # If you are going to have ionic elements in the solvent
    # special care needs to be taken to
    # Ensure parameters for ions are sourced correctly!

    TLEAP_SCRPT = '''#!/bin/bash
    cat > "${model}".tleap <<EOF
    source leaprc.water.${model}

    # Create a system with one residue of OPC water
    mol = sequence { WAT }

    # Save the system to a parameter and coordinate file
    saveAmberParm mol ${model}.prmtop ${model}.inpcrd

    quit
    EOF

    tleap -f "${model}".tleap
    '''

    # Run the script in the current Bash environment
    process = subprocess.run(['bash', '-c', TLEAP_SCRPT], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    # Print the output
    print("Output:", process.stdout)

    # Print any errors
    if process.stderr:
        print("Errors:", process.stderr)

    parm = parmed.amber.LoadParm("${model}.prmtop", xyz="${model}.inpcrd")

    print(parm.bonds[1].atom1.type)
    print(parm.bonds[1].atom2.type)
    print(parm.bonds[1].type)
    print(parm.LJ_types)
    print(parm.LJ_radius)
    # Create a list of indices
    indices_list = list(range(len(parm.LJ_types)))
    print(parm.LJ_depth)

    for LJ_params in zip(parm.LJ_types,parm.LJ_radius,parm.LJ_depth):

        if (LJ_params[1]==0):
            print("embedding", LJ_params, " radius")
            print(LJ_params[0])
            my_type = LJ_params[0]

            matching_bond = next((bond for bond in parm.bonds \
            if (bond.atom1.type == my_type and bond.atom2.type != my_type) \
            or (bond.atom2.type == my_type and bond.atom1.type != my_type)), None)
            bondLength = matching_bond.type.req
            print("bondLength",bondLength)

            other_atom_type = matching_bond.atom2.type if \
            matching_bond.atom1.type == my_type else matching_bond.atom1.type \
            if matching_bond.atom2.type == my_type else None

            other_radius = parm.LJ_radius[parm.LJ_types[other_atom_type]-1]
            other_sigma = other_radius*(2 ** (-1/6))
            new_sigma = other_sigma - bondLength
            parm.LJ_radius[parm.LJ_types[my_type]-1]= new_sigma/(2 ** (-1/6))
        if (LJ_params[2]==0):
            print("embedding", LJ_params, " depth")
            my_type = LJ_params[0]

            matching_bond = next((bond for bond in parm.bonds \
            if (bond.atom1.type == my_type and bond.atom2.type != my_type) \
            or (bond.atom2.type == my_type and bond.atom1.type != my_type)), None)
            bondLength = matching_bond.type.req
            print("bondLength",bondLength)

            other_atom_type = matching_bond.atom2.type if \
            matching_bond.atom1.type == my_type else matching_bond.atom1.type \
            if matching_bond.atom2.type == my_type else None

            other_depth = parm.LJ_depth[parm.LJ_types[other_atom_type]-1]
            new_depth = other_depth * 0.1
            parm.LJ_depth[parm.LJ_types[my_type]-1]= new_depth
    print(parm.LJ_radius)
    print(parm.LJ_depth)

    parm.recalculate_LJ()
    parm.write_parm("c{}.prmtop".format("${model}".upper()))

    print(parm.pointers)
    print(parm.ptr)
    print(parm._AMBERPARM_ATTRS)
    print(dir(parm))

    get_name = lambda ltype: next(atom.name for atom in parm.residues[0].atoms if atom.type == ltype)
    get_mass = lambda ltype: next(atom.mass for atom in parm.residues[0].atoms if atom.type == ltype)
    get_charge = lambda ltype: next(atom.charge for atom in parm.residues[0].atoms if atom.type == ltype)
    get_multi = lambda ltype: sum(1 for atom in parm.residues[0].atoms if atom.type == ltype)
    #get_coords = lambda ltype: next(atom.charge for atom in parm.residues[0].atoms if atom.type == ltype)

    emptyMDL = parmed.amber.AmberFormat()
    emptyMDL.charge_flag="CHG"
    #emptyMDL.add_flag(flag_name='TITLE',flag_format=str(parm.formats['TITLE']),data=parm.parm_data['TITLE'])
    emptyMDL.add_flag(flag_name='TITLE',flag_format=str(parm.formats['TITLE']),data=["c{}".format("${model}".upper())])
    emptyMDL.add_flag(flag_name='POINTERS',flag_format=str(parm.formats['POINTERS']),data=[parm.pointers['NATOM'],parm.pointers['NTYPES']])
    emptyMDL.add_flag(flag_name='ATMTYP',flag_format=str(parm.formats['ATOM_TYPE_INDEX']),data=[value for value in parm.LJ_types.values()])
    emptyMDL.add_flag(flag_name='ATMNAME',flag_format=str(parm.formats['TITLE']),data=[get_name(ltype) for ltype in parm.LJ_types.keys()])
    emptyMDL.add_flag(flag_name='MASS',flag_format=str(parm.formats['RADII']),data=[get_mass(ltype) for ltype in parm.LJ_types.keys()])
    emptyMDL.add_flag(flag_name='CHG',flag_format=str(parm.formats['RADII']),data=[get_charge(ltype) for ltype in parm.LJ_types.keys()])
    emptyMDL.add_flag(flag_name='LJEPSILON',flag_format=str(parm.formats['RADII']),data=[parm.LJ_depth[value-1] for value in parm.LJ_types.values()])
    emptyMDL.add_flag(flag_name='LJSIGMA',flag_format=str(parm.formats['RADII']),data=[parm.LJ_radius[value-1] for value in parm.LJ_types.values()])
    emptyMDL.add_flag(flag_name='MULTI',flag_format=str(parm.formats['POINTERS']),data=[get_multi(ltype) for ltype in parm.LJ_types.keys()])
    emptyMDL.add_flag(flag_name='COORD',flag_format=str(parm.formats['RADII']),data=np.concatenate(list(np.array(item).flatten() for item in parm.get_coordinates())))

    emptyMDL.write_parm("c{}.mdl".format("${model}".upper()))
    quit()

    """

}

workflow build_ligands {
    take:
    extract_database_ch
    main:
    // Process each JSON file asynchronously
    build_ligand(extract_database_ch)
    emit:
    //molecule = build_ligand.out.molecule
    //prm = build_ligand.out.prm
    //crd = build_ligand.out.crd
    system = build_ligand.out.system
}

workflow build_solvents {
    take:
    solv_temp_pairs
    main:
    //build_solvent(solv_temp_pairs)
    build_water_model(solv_temp_pairs) | build_solvent
    emit:
    xvv = build_solvent.out.paths.flatten().filter { file -> file.name.endsWith("xvv") }
    model = build_solvent.out.model
    temperature = build_solvent.out.temperature
    solvent = build_solvent.out.solvent
}
