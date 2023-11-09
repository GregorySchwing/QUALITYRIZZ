#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process build_ligand {
    container "${params.container__openff_toolkit}"
    publishDir "${params.output_folder}/${params.database}/systems/${pathToJson.baseName}", mode: 'copy', overwrite: true

    debug true
    input:
    path pathToJson
    output:
    path "${pathToJson.baseName}.txt"

    script:
    """
    echo JSON file name: ${pathToJson}
    touch ${pathToJson.baseName}.txt
    """
}

workflow build_ligands {
    take:
    extract_database_ch
    main:
    // Process each JSON file asynchronously
    build_ligand(extract_database_ch)
}