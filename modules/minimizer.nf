#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process minimize_ligand {
    container "${params.container__biobb_amber}"
    publishDir "${params.output_folder}/${params.database}/minimizations/${model}_${T}_${molecule}", mode: 'copy', overwrite: false

    debug true
    input:
    tuple val(model), val(T), val(molecule)
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

workflow minimize_ligand {
    take:
    extract_database_ch
    main:
    // Process each JSON file asynchronously
    minimize_ligand(extract_database_ch)


}