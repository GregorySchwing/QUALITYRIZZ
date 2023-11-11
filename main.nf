#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// All of the default parameters are being set in `nextflow.config`

// Import sub-workflows
include { extract_database } from './modules/database_reader'
include { build_ligands } from './modules/system_builder'
include { build_solvents } from './modules/system_builder'
include { minimize_ligands } from './modules/minimizer'
include { rism_solvation } from './modules/rism'
include { analyze_solvation } from './modules/analysis'


// Function which prints help message text
def helpMessage() {
    log.info"""
Usage:

nextflow run -profile docker . --param_name nextflow.config

Required Arguments:

  Input Data:
  --database            Name of database to run (default FreeSolv).

Optional Arguments:
    --output_folder       Folder for output files (default $projectDir/results).
    """.stripIndent()
}


// Main workflow
workflow {

    // Show help message if the user specifies the --help flag at runtime
    // or if any required params are not provided
    if ( params.help){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 1
    }

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

         """
         .stripIndent()

    if ( params.database ){

        // Set up a channel from the pairs of files found with that pattern

        pathToDataBase = ""
        if (params.database == "FreeSolv"){
            pathToDataBase = "$projectDir/databases/FreeSolv/database.pickle"
        }

        database = Channel
            .fromPath( pathToDataBase )
        database.view()
        extract_database(
            database
        )
        nc = extract_database.out.json.flatten()
        build_ligands(nc)
        solvent = Channel.from( [["cSPCE","298.15"]] )
        build_solvents(solvent)
        minimize_ligands(build_ligands.out.system,build_solvents.out.solvent)
        rism_solvation(minimize_ligands.out.minimized_system)
        results = rism_solvation.out.json.collect()
        results.view()
        analyze_solvation(results,database)
    }
}