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
    path("ligand.crd"), emit: crd
    tuple val(pathToJson.baseName), path("ligand.prmtop"), path("ligand.crd"), emit: system
    script:
    """
    #!/usr/bin/env python
    import json
    # Open and read the JSON file
    with open("${pathToJson}", "r") as file:
        print("reading ${pathToJson}")
        data = json.load(file)
    key=next(iter(data))
    smiles=data[key]["SMILES"]
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
    interchange.to_inpcrd("ligand.crd")
    
    """
}


process build_ligand_list {
    container "${params.container__openff_toolkit}"
    publishDir "${params.output_folder}/${params.database}/systems/${molecule}", mode: 'copy', overwrite: true

    debug false
    input:
    tuple  val(molecule), val(SMILES), val(reference)
    output:
    val(molecule), emit: molecule
    path("ligand.prmtop"), emit: prm
    path("ligand.crd"), emit: crd
    tuple val(molecule), path("ligand.prmtop"), path("ligand.crd"), emit: system
    script:
    """
    #!/usr/bin/env python

    smiles="$SMILES"
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
    interchange.to_inpcrd("ligand.crd")
    
    """
}

process build_water_parameters {
    container "${params.container__openff_toolkit}"
    publishDir "${params.output_folder}/${params.database}/water_parameters/${model}", mode: 'copy', overwrite: true

    debug false
    input:
    tuple val(model), val(temperature)
    output:
    tuple val(model), val(temperature), path("*.mdl")

    script:
    """
    #!/usr/bin/env python
    # Imports from the toolkit
    from openff.toolkit import ForceField, Molecule, Topology
    from openff.units import Quantity, unit
    from openff.interchange import Interchange
    from openff.toolkit.typing.engines.smirnoff import ForceField
    from openff.toolkit.typing.engines import smirnoff

    # Many thanks: openff-toolkit/examples/virtual_sites/vsite_showcase.ipynb
    # And https://www.calculator.net/triangle-calculator.html
    def build_water(OH_bond_length=1, HH_bond_length=1.63298086):
        from math import cos, sin, acos, degrees
        calc_theta = lambda a, b: degrees(acos(1 - ((a * a) / (2 * b * b)))) if a != 0 else None    
        water_reference = Molecule.from_mapped_smiles("[O:1]([H:2])[H:3]")
        water_reference.atoms[0].name = "O"
        water_reference.atoms[1].name = "H1"
        water_reference.atoms[2].name = "H2"

        # Add ideal TIP5P geometry
        bond_length = Quantity(OH_bond_length, unit.angstrom)
        print("theta",calc_theta(HH_bond_length,OH_bond_length))
        print("rOH",bond_length)

        theta = Quantity(calc_theta(HH_bond_length,OH_bond_length), unit.degree).to(unit.radian)
        water_reference.add_conformer(
            bond_length
            * Quantity(
                [
                    [0.0, 0.0, 0.0],
                    [-sin(theta / 2), cos(theta / 2), 0.0],
                    [sin(theta / 2), cos(theta / 2), 0.0],
                ]
            )
        )

        return water_reference

    #print(smirnoff.get_available_force_fields())
    print("${model}")
    water_interchange = ForceField("${model}").create_interchange(Molecule.from_mapped_smiles(
        "[O:1]([H:2])[H:3]"
    ).to_topology())
    OH_bond, HH_bond = water_interchange.collections['Constraints'].potentials.values()

    # Build correct water molecule
    water_reference = build_water(OH_bond.parameters['distance'],HH_bond.parameters['distance'])

    # Scaffold interchange, with bond information needed to construct water
    interchange = ForceField("openff-2.0.0.offxml").create_interchange(water_reference.to_topology())


    # Assigns desired forcefield parameters to interchange
    interchange.collections['vdW'].potentials=water_interchange.collections['vdW'].potentials
    interchange.collections['Electrostatics'].potentials=water_interchange.collections['Electrostatics'].potentials

    # Tweak forcefield if you wish
    if (True):
        from openff.interchange.models import TopologyKey
        oxygen = TopologyKey(atom_indices=(0,))
        hydrogen = TopologyKey(atom_indices=(1,))
        oxygen_pot_key = interchange.collections['vdW'].key_map[oxygen]
        hydrogen_pot_key = interchange.collections['vdW'].key_map[hydrogen]
        interchange.collections['vdW'].potentials[hydrogen_pot_key].parameters['sigma']= interchange.collections['vdW'].potentials[oxygen_pot_key].parameters['sigma'] - 2*OH_bond.parameters['distance']
        interchange.collections['vdW'].potentials[hydrogen_pot_key].parameters['epsilon']= interchange.collections['vdW'].potentials[oxygen_pot_key].parameters['epsilon'] * 0.1
        print("sigma_Hy",interchange.collections['vdW'].potentials[hydrogen_pot_key].parameters['sigma'])
        print("epsilon_Hy",interchange.collections['vdW'].potentials[hydrogen_pot_key].parameters['epsilon'])
        print("sigma_O",interchange.collections['vdW'].potentials[oxygen_pot_key].parameters['sigma'])
        print("epsilon_Oxy",interchange.collections['vdW'].potentials[oxygen_pot_key].parameters['epsilon'])


    interchange.to_prmtop("water.prmtop")
    interchange.to_inpcrd("water.crd")
    from parmed.amber import LoadParm
    parm=LoadParm("water.prmtop","water.crd")
    parm.write_mdl("water.mdl")
    """
}


process build_solvent_parameters {
    container "${params.container__openff_toolkit}"
    publishDir "${params.output_folder}/${params.database}/solvent_parameters/${model}", mode: 'copy', overwrite: true

    debug true
    input:
    tuple val(model), val(SMILES), val(FF), val(temperature)
    output:
    tuple val(model), val(temperature), path("solvent.prmtop"), path("solvent.crd")

    script:
    """
    #!/usr/bin/env python
    # Imports from the toolkit
    from openff.toolkit import ForceField, Molecule, Topology
    from openff.units import Quantity, unit
    from openff.interchange import Interchange
    from openff.toolkit.typing.engines.smirnoff import ForceField
    from openff.toolkit.typing.engines import smirnoff

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


    # Many thanks: openff-toolkit/examples/virtual_sites/vsite_showcase.ipynb
    # And https://www.calculator.net/triangle-calculator.html
    def build_water(OH_bond_length=1, HH_bond_length=1.63298086):
        from math import cos, sin, acos, degrees
        calc_theta = lambda a, b: degrees(acos(1 - ((a * a) / (2 * b * b)))) if a != 0 else None    
        water_reference = Molecule.from_mapped_smiles("[O:1]([H:2])[H:3]")
        water_reference.atoms[0].name = "O"
        water_reference.atoms[1].name = "H1"
        water_reference.atoms[2].name = "H2"

        # Add ideal TIP5P geometry
        bond_length = Quantity(OH_bond_length, unit.angstrom)
        print("theta",calc_theta(HH_bond_length,OH_bond_length))
        print("rOH",bond_length)

        theta = Quantity(calc_theta(HH_bond_length,OH_bond_length), unit.degree).to(unit.radian)
        water_reference.add_conformer(
            bond_length
            * Quantity(
                [
                    [0.0, 0.0, 0.0],
                    [-sin(theta / 2), cos(theta / 2), 0.0],
                    [sin(theta / 2), cos(theta / 2), 0.0],
                ]
            )
        )

        return water_reference

    #print(smirnoff.get_available_force_fields())

    if ("$model"=="WATER"):
        #print(smirnoff.get_available_force_fields())
        print("${FF}")
        water_interchange = ForceField("${FF}").create_interchange(Molecule.from_mapped_smiles(
            "[O:1]([H:2])[H:3]"
        ).to_topology())
        OH_bond, HH_bond = water_interchange.collections['Constraints'].potentials.values()

        # Build correct water molecule
        water_reference = build_water(OH_bond.parameters['distance'],HH_bond.parameters['distance'])

        # Scaffold interchange, with bond information needed to construct water
        interchange = ForceField("openff-2.0.0.offxml").create_interchange(water_reference.to_topology())


        # Assigns desired forcefield parameters to interchange
        interchange.collections['vdW'].potentials=water_interchange.collections['vdW'].potentials
        interchange.collections['Electrostatics'].potentials=water_interchange.collections['Electrostatics'].potentials

        # Tweak forcefield if you wish
        if (True):
            from openff.interchange.models import TopologyKey
            oxygen = TopologyKey(atom_indices=(0,))
            hydrogen = TopologyKey(atom_indices=(1,))
            oxygen_pot_key = interchange.collections['vdW'].key_map[oxygen]
            hydrogen_pot_key = interchange.collections['vdW'].key_map[hydrogen]
            interchange.collections['vdW'].potentials[hydrogen_pot_key].parameters['sigma']= interchange.collections['vdW'].potentials[oxygen_pot_key].parameters['sigma'] - 2*OH_bond.parameters['distance']
            interchange.collections['vdW'].potentials[hydrogen_pot_key].parameters['epsilon']= interchange.collections['vdW'].potentials[oxygen_pot_key].parameters['epsilon'] * 0.1
            print("sigma_Hy",interchange.collections['vdW'].potentials[hydrogen_pot_key].parameters['sigma'])
            print("epsilon_Hy",interchange.collections['vdW'].potentials[hydrogen_pot_key].parameters['epsilon'])
            print("sigma_O",interchange.collections['vdW'].potentials[oxygen_pot_key].parameters['sigma'])
            print("epsilon_Oxy",interchange.collections['vdW'].potentials[oxygen_pot_key].parameters['epsilon'])

        interchange.to_prmtop("solvent.prmtop")
        interchange.to_inpcrd("solvent.crd")
        from parmed.amber import LoadParm
        parm=LoadParm("solvent.prmtop","solvent.crd")
        parm.write_mdl("solvent.mdl")
    else:
        forcefield = ForceField("${FF}")
        openff_mol = Molecule.from_smiles("${SMILES}",allow_undefined_stereo=True)
        rdmol = openff_mol.to_rdkit()
        rdmol3D = embed(rdmol,123)
        # Create a PDB file
        writer = rdmolfiles.MolToPDBFile(rdmol3D, "solvent.pdb")
        # Load the topology from a PDB file and `Molecule` objects
        topology = Topology.from_pdb(
            "solvent.pdb",
            unique_molecules=[openff_mol],
        )

        interchange = Interchange.from_smirnoff(
            force_field=forcefield,
            topology=topology,
        )
        interchange.minimize()
        interchange.to_prmtop("solvent.prmtop")
        interchange.to_inpcrd("solvent.crd")
        from parmed.amber import LoadParm
        parm=LoadParm("solvent.prmtop","solvent.crd")
        parm.write_mdl("solvent.mdl")

    """
}


process minimize_solvent {
    container "${params.container__biobb_amber}"
    publishDir "${params.output_folder}/${params.database}/solvent_minimized/${model}_${temperature}", mode: 'copy', overwrite: false

    debug false
    input:
    tuple val(model), val(temperature), path(prm), path(crd)
    output:
    path("sander.*"), emit: paths
    tuple val(model), val(temperature), path(prm), path("sander.n_min.rst7"), emit: minimized_solvent
    maxRetries 20
    script:
    """
    #!/usr/bin/env python
    print("Hello from ${model} ${temperature} ${prm} ${crd}")
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
    }

    # Create and launch bb
    sander_mdrun(input_top_path="${prm}",
                input_crd_path="${crd}",
                output_traj_path=output_n_min_traj_file,
                output_rst_path=output_n_min_rst_file,
                output_log_path=output_n_min_log_file,
                properties=prop)

    # Import module
    from biobb_amber.process.process_minout import process_minout

    # Create prop dict and inputs/outputs
    output_n_min_dat_file = 'sander.n_min.energy.dat'
    prop = {
        "terms" : ['ENERGY']
    }

    # Create and launch bb
    process_minout(input_log_path=output_n_min_log_file,
                output_dat_path=output_n_min_dat_file,
                properties=prop)

    from biobb_amber.ambpdb.amber_to_pdb import amber_to_pdb
    output_n_min_pdb_file = "sander.n_min.pdb"
    prop = {
        'remove_tmp': True,
        'check_extensions' : False,
    }
    amber_to_pdb(input_top_path="${prm}",
                input_crd_path=output_n_min_rst_file,
                output_pdb_path=output_n_min_pdb_file,
                properties=prop)
    """
}



process build_mdl {
    container "${params.container__openff_toolkit}"
    publishDir "${params.output_folder}/${params.database}/mdl/${model}_${temperature}", mode: 'copy', overwrite: true

    debug true
    input:
    tuple val(model), val(temperature), path(prm), path(rst)
    output:
    tuple val(model), val(temperature), path("solvent.mdl")

    script:
    """
    #!/usr/bin/env python
    # Imports from the toolkit
    from parmed.amber import LoadParm
    parm=LoadParm("${prm}","${rst}")
    parm.write_mdl("solvent.mdl")

    """
}


process build_solvent {
    container "${params.container__biobb_amber}"
    publishDir "${params.output_folder}/${params.database}/1DRISM/${model}_${T}", mode: 'copy', overwrite: true

    debug false
    input:
    tuple val(model), val(T), path(mdl)
    output:
    path("${model}_${T}*.*"), emit: paths
    val(model), emit: model
    val(T), emit: temperature
    path("solvParams.json"), emit: solvParams
    tuple val(model), val(T), path("${model}_${T}_HNC_*.xvv"), emit: solvent optional true
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
        KSAVE={NUMSTEPS}, KOUT=0, maxstep={NUMSTEPS},       !Check pointing and iterations
        SMEAR=1, ADBCOR=0.5,                   !Electrostatics
        TEMPERATURE={temp}, DIEps={diel},      !bulk solvent properties
        NSP=1,
        entropicDecomp=1                       ! While debugging convergence set 0
    /
        &SPECIES                               !SPC water
        DENSITY={conc}d0,
        MODEL="${mdl}"
        !MODEL="/usr/local/dat/rism1d/mdl/cSPCE.mdl"
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

    closure_list = ["KH", "PSE1", "PSE2", "PSE3","HNC"]
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

process build_solvent_no_boot {
    container "${params.container__biobb_amber}"
    publishDir "${params.output_folder}/${params.database}/1DRISM/${model}_${T}", mode: 'copy', overwrite: true

    debug false
    input:
    tuple val(model), val(T), path(mdl)
    output:
    path("${model}_${T}*.*"), emit: paths
    val(model), emit: model
    val(T), emit: temperature
    path("solvParams.json"), emit: solvParams
    tuple val(model), val(T), path("${model}_${T}.xvv"), emit: solvent
    shell:

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
    cat > {name1d}.inp <<EOF
    &PARAMETERS
        THEORY='{rism1d}', CLOSUR='{closure}',           !Theory
        NR=16384, DR=0.025,                    !Grid
        OUTLST='xCGT', rout=0,                 !Output
        NIS=20, DELVV=0.3, TOLVV=1.e-12,       !MDIIS
        KSAVE=-1, KSHOW=1, maxstep=10000,      !Check pointing and iterations
        SMEAR=1, ADBCOR=0.5,                   !Electrostatics
        TEMPER={temp}, DIEps={diel},           !bulk solvent properties
        NSP=1
    /
        &SPECIES                               !SPC water
        DENSITY={conc}d0,
        MODEL="/usr/local/dat/rism1d/mdl/cSPCE.mdl"
    /
    EOF

    rism1d "${model}_${T}" > "${model}_${T}.out"
    '''
    
    T=${T}
    smodel="${model}"
    rism1d="DRISM"
    closure="PSE3"
    rism1d_name = '{smodel}_{temp}'.format(smodel=smodel,temp=T)
    diel = dieps
    conc = conc
    input_string=SOLV_SUCEPT_SCRPT.format(temp=T, diel=diel, conc=conc,\
                        smodel=smodel, rism1d=rism1d,\
                        closure=closure.upper(),\
                        name1d=rism1d_name)

    import subprocess
    # Run the script in the current Bash environment
    process = subprocess.run(['bash', '-c', input_string], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    # Print the output
    print("Output:", process.stdout)

    # Print any errors
    if process.stderr:
        print("Errors:", process.stderr)
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
    //build_water_parameters(solv_temp_pairs)
    build_solvent_parameters(solv_temp_pairs)
    | minimize_solvent
    
    build_mdl(minimize_solvent.out.minimized_solvent)
    | build_solvent

    
    emit:
    xvv = build_solvent.out.paths.flatten().filter { file -> file.name.endsWith("xvv") }
    model = build_solvent.out.model
    temperature = build_solvent.out.temperature
    solvent = build_solvent.out.solvent
    
}
