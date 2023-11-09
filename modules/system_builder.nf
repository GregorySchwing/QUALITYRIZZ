#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2


process echoJsonFileName {
    input:
    path json_file

    script:
    """
    echo "JSON file name: ${json_file}"
    """
}

workflow build_ligands {
    take:
    extract_database_ch
    main:
    // Process each JSON file asynchronously
    echoJsonFileName(extract_database_ch)
}