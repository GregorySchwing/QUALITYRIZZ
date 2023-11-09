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
    path("ligand.prmtop"), emit: prm
    path("ligand.inpcrd"), emit: crd

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
    openff_mol = Molecule.from_smiles(smiles)
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

process get_conc_dieps {
    container "${params.container__openff_toolkit}"

    debug true
    input:
    tuple val(model), val(T)
    output:
    path('solvParams.json'), emit: json
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
    """
}

process build_solvent {
    container "${params.container__biobb_amber}"
    //publishDir "${params.output_folder}/${params.database}/solvents/${model}_${T}", mode: 'copy', overwrite: true

    debug true
    input:
    tuple val(model), val(T), path('solvParams.json')

    script:
    """
    #!/usr/bin/env python
    ## RISM-related globals ##

    import json
    # Open and read the JSON file
    with open("solvParams.json", "r") as file:
        print("reading 'solvParams.json'")
        data = json.load(file)

    T=${T}
    smodel="${model}"
    rism1d="DRISM"
    closure="PSE3"
    rism1d_name = '{smodel}_{temp}'.format(smodel=smodel,temp=T)
    diel = data['diel']
    conc = data['conc']

    print(diel,conc)

    """
    //template 'rism1d.sh'

}

workflow build_ligands {
    take:
    extract_database_ch
    main:
    // Process each JSON file asynchronously
    build_ligand(extract_database_ch)
    emit:
    prm = build_ligand.out.prm
    crd = build_ligand.out.crd

}

workflow build_solvents {
    take:
    solv_temp_pairs
    main:
    // Process each JSON file asynchronously
    get_conc_dieps(solv_temp_pairs)
    get_conc_dieps.out.json.view()

    //build_solvent(get_conc_dieps.out.json)
}