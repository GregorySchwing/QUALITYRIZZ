#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// All of the default parameters are being set in `nextflow.config`

// Import sub-workflows
include { extract_database } from './modules/database_reader'
include { extract_database_deepchem } from './modules/database_reader'
include { extract_database_channel } from './modules/database_reader'


include { build_ligands } from './modules/system_builder'
include { build_solvents } from './modules/system_builder'

include { minimize_ligands } from './modules/minimizer'
include { rism_solvation } from './modules/rism'
include { analyze_list } from './modules/analysis'


// Function which prints help message text
def helpMessage() {
    log.info"""
Usage:

nextflow run [-profile local/docker/singularity/slurm] . [--database /Path/To/CSV --reference_col expt]

# FreeSolv
nextflow run . -profile docker

# User defined database
nextflow run . -profile docker --database USER --database_path databases/PFAS/PFAS_MAP_INPUT.csv --id_col NAME --structure_col SMILES --reference_col PROPERTY

Required Arguments:

  Input Data:
  --database            Name of database to run (options: FreeSolv - Default, User).
  --database_path       /Path/To/CSV.
  --id_col              Column in user defined csv to use as unique keys (default compoundid).
  --structure_col       Column in user defined csv to build ligands (default SMILES).
  --reference_col       Column in user defined csv to compare with RISM-HFE (default experimentalvalue(kcal/mol)).

Optional Arguments:
    --output_folder       Folder for output files (default $projectDir/results).
    --solute_samples      Number of solutes in solute database to run (default -1).
    --solvent_samples     Number of solutes in solvent list to run (default -1).

    """.stripIndent()
}


// Main workflow
workflow {

def myString = """ 
.----------------.  .----------------.  .----------------.  .----------------. 
| .--------------. || .--------------. || .--------------. || .--------------. |
| |  _______     | || |     _____    | || |   ________   | || |   ________   | |
| | |_   __ \\    | || |    |_   _|   | || |  |  __   _|  | || |  |  __   _|  | |
| |   | |__) |   | || |      | |     | || |  |_/  / /    | || |  |_/  / /    | |
| |   |  __ /    | || |      | |     | || |     .'.' _   | || |     .'.' _   | |
| |  _| |  \\ \\_  | || |     _| |_    | || |   _/ /__/ |  | || |   _/ /__/ |  | |
| | |____| |___| | || |    |_____|   | || |  |________|  | || |  |________|  | |
| |              | || |              | || |              | || |              | |
| '--------------' || '--------------' || '--------------' || '--------------' |
 '----------------'  '----------------'  '----------------'  '----------------' """





log.info """\
         ${myString}
=============================================================================
         output_folder          : ${params.output_folder}
         database               : ${params.database}
         database_path          : ${params.database_path}
         solvent_path           : ${params.solvent_path}
         id_col                 : ${params.id_col}
         structure_col          : ${params.structure_col}
         reference_col          : ${params.reference_col}
         solute_samples         : ${params.solute_samples}
         solvent_samples        : ${params.solvent_samples}
         """
         .stripIndent()


    // Supported charge types : 'am1bcc', 'am1-mulliken', 'gasteiger', 'resp'
    // Show help message if the user specifies the --help flag at runtime
    // or if any required params are not provided
    if ( params.help || params.database == null){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 1
    }

    if ( params.database ){

        //waterModels = ["tip3p_fb-1.1.1.offxml"]
        waterModels = ["tip3p_fb-1.1.1.offxml","tip3p-1.0.1.offxml","opc3-1.0.1.offxml","spce-1.0.0.offxml"]
        temperatures = ["298.15"]
        partial_charge_method = ["gasteiger","am1bcc","am1-mulliken","RESP"]
        input = Channel.fromPath( params.database_path ).splitCsv(header: true,limit: 1).map { 
            row -> [row."${params.id_col}", row."${params.structure_col}", row."${params.reference_col}"]
        }

        input_dict2 = Channel.fromPath(params.database_path)
            .splitCsv(header: true, limit: params.solute_samples)
            .flatMap { row ->
                def tumor_reads =  '{' + "\"${params.reference_col}\"" + ': ' + "\"${row."${params.reference_col}"}\"" + ', ' +  
                                        "\"${params.structure_col}\"" + ': ' + "\"${row."${params.structure_col}"}\"" + '}'
                ['{' + "\"${row."${params.id_col}"}\"" + ': ' + tumor_reads  + '}']
            }
        //input_dict2.view()
 
        solventListChannel = Channel.fromPath( params.solvent_path ).splitCsv(header: true,limit: 4).map { 
            row -> [row.NAME, row.SMILES, row.FF, row.TEMP, row.DIEPS, row.DENSITY]
        }
        solventListChannel.view()
        build_solvents(solventListChannel) 
        database=extract_database_channel(input_dict2)
        molecules_and_charges = database.combine(partial_charge_method)
        results = build_ligands(molecules_and_charges)
            | combine(build_solvents.out.solvent) 
            | minimize_ligands
            | rism_solvation
            | collect

        analyze_list(results,database.collect())
    } else {
        helpMessage()
    }
}