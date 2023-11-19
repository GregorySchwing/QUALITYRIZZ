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

        pathToDataBase = ""
        if (params.database == "FreeSolv"){
            pathToDataBase = "$projectDir/databases/FreeSolv/database.pickle"
        }
        database = Channel
            .fromPath( pathToDataBase )

        //waterModels = ["spce-1.0.0.offxml","tip3p_fb-1.1.1.offxml"]
        waterModels = ["tip3p_fb-1.1.1.offxml","tip3p-1.0.1.offxml","opc3-1.0.1.offxml","spce-1.0.0.offxml"]
        temperatures = ["298.15"]
        waterChannel = Channel.from( waterModels )
        temperatureChannel = Channel.from( temperatures )
        solventChannel = waterChannel.combine(temperatureChannel)
        //solventChannel.view()
        build_solvents(solventChannel) 
                //| combine(temperatureChannel)
        // Channel holds the indices to sample.
        // In future, will use ML to determine which samples along with parameters
        // In an iterative loop.
        /**
        Channel.of(0..2)
                | buffer(size: 642, remainder: true)
                | extract_database_deepchem 
        */
        results = extract_database(database)
                | flatten 
                | build_ligands 
                | combine(build_solvents.out.solvent) 
                | minimize_ligands
                | rism_solvation
                | collect
        analyze_mobley_solvation(results,database)
                //| analyze_solvation
    }
}