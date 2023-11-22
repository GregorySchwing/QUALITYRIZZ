#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// All of the default parameters are being set in `nextflow.config`

// Import sub-workflows
include { extract_database } from './modules/database_reader'
include { extract_database_deepchem } from './modules/database_reader'


include { build_ligands } from './modules/system_builder'
include { build_solvents } from './modules/system_builder'

include { minimize_ligands } from './modules/minimizer'
include { rism_solvation } from './modules/rism'
include { analyze_solvation } from './modules/analysis'
include { analyze_mobley_solvation } from './modules/analysis'


// Function which prints help message text
def helpMessage() {
    log.info"""
Usage:

nextflow run [-profile local/docker/singularity/slurm] . [--database /Path/To/CSV --reference_col expt]

Required Arguments:

  Input Data:
  --database            Name of database to run (options: FreeSolv - Default, /Path/To/CSV (REQUIRED COLUMNS: NAME,SMILES,reference_col} )).
  --id_col              Column in user defined csv to use as unique keys (default compoundid).
  --structure_col       Column in user defined csv to build ligands (default SMILES).
  --reference_col       Column in user defined csv to compare with RISM-HFE (default experimentalvalue(kcal/mol)).

Optional Arguments:
    --output_folder       Folder for output files (default $projectDir/results).
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
         id_col                 : ${params.id_col}
         structure_col          : ${params.structure_col}
         reference_col          : ${params.reference_col}
         """
         .stripIndent()

    // Show help message if the user specifies the --help flag at runtime
    // or if any required params are not provided
    if ( params.help || params.database == null){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 1
    }

    if ( params.database ){

        pathToDataBase = ""
        id_col = ""
        structure_col = ""
        reference_col = ""

        if (params.database == "FreeSolv"){
            pathToDataBase = "$projectDir/databases/FreeSolv/database.txt"
            database = Channel.fromPath( pathToDataBase )
        } else {
            databaseName = params.database.split("\\.")[0].split(File.separator)[-1]
            pathToDataBase = params.database
            reference_col = params.reference_col
            c2 = Channel.fromPath( pathToDataBase ).splitCsv(header: true)
            c2.view { row -> "${row.NAME} - ${row.SMILES} - ${row."${reference_col}"}" }
        }
        return


        //waterModels = ["spce-1.0.0.offxml","tip3p_fb-1.1.1.offxml"]
        waterModels = ["tip3p_fb-1.1.1.offxml","tip3p-1.0.1.offxml","opc3-1.0.1.offxml","spce-1.0.0.offxml"]
        temperatures = ["298.15"]
        waterChannel = Channel.from( waterModels )
        temperatureChannel = Channel.from( temperatures )
        solventChannel = waterChannel.combine(temperatureChannel)
        solventChannel.view()
        build_solvents(solventChannel) 
        results = extract_database(database)
            | flatten
            | build_ligands 
            | combine(build_solvents.out.solvent) 
            | minimize_ligands
            | rism_solvation
            | collect
        analyze_mobley_solvation(results,database)
                //| analyze_solvation
    } else {
        helpMessage()
    }
}