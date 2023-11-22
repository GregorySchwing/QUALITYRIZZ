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

        //waterModels = ["spce-1.0.0.offxml","tip3p_fb-1.1.1.offxml"]
        waterModels = ["tip3p_fb-1.1.1.offxml","tip3p-1.0.1.offxml","opc3-1.0.1.offxml","spce-1.0.0.offxml"]
        temperatures = ["298.15"]
        waterChannel = Channel.from( waterModels )
        temperatureChannel = Channel.from( temperatures )
        solventChannel = waterChannel.combine(temperatureChannel)
        solventChannel.view()
        build_solvents(solventChannel) 
        c2 = Channel.fromPath( params.database_path ).splitCsv(header: true,limit: 10).map { 
            row -> [row."${params.id_col}", row."${params.structure_col}"]
            //row -> [row."${params.id_col}", row."${params.structure_col}", row."${params.reference_col}"]
        }
            | build_ligands 
            | combine(build_solvents.out.solvent) 
            | minimize_ligands
            | rism_solvation
        //    | collect
        //analyze_mobley_solvation(results,database)
                //| analyze_solvation*/
    } else {
        helpMessage()
    }
}