#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process depickle {
    container "${params.container__openff_toolkit}"
    publishDir "${params.output_folder}/${params.database}/json", mode: 'copy', overwrite: true

    debug false
    input:
    path pathToDatabase
    
    output:
    path('*.json'), emit: json

    script:
    """
    #!/usr/bin/env python
    import pickle
    import sys
    import json
    import pandas as pd
    from scipy.stats import pearsonr  # Import the pearsonr function
    print(f"Python Version: {sys.version}")
    # Use the pandas read_pickle function to read the DataFrame.
    df = pd.read_pickle("$pathToDatabase")
    #print(df['mobley_5708811'])

    # Write the dictionary to the file in JSON format
    counter = 0
    for key in df.keys():
        if counter > 2:
            break
        
        with open("{key}.json".format(key=key), 'w') as file:
            json.dump(df[key], file)
        
        counter = counter + 1

    """
}


process sample_deepchem {
    container "${params.container__deep_chem}"
    //publishDir "${params.output_folder}/${params.database}/json", mode: 'copy', overwrite: true

    debug false
    input:
    val indices
    output:
    path('*.json'), emit: json

    script:
    """
    #!/usr/bin/env python
    from deepchem.molnet import load_freesolv
    import json
    print($indices)
    tasks, datasets, transformers = load_freesolv(splitter=None)
    train = datasets[0]
    x,y,w,ids = train.X, train.y, train.w, train.ids
    for index in $indices:
        df = {"smiles": ids[index]}
        with open("{key}.json".format(key=index), 'w') as file:
            json.dump(df, file)
    """
}

workflow extract_database{

    take:
    database_pickle_ch
    main:
    depickle(database_pickle_ch) 
    emit:
    json = depickle.out.json

}

workflow extract_database_deepchem{

    take:
    database_pickle_ch
    main:
    sample_deepchem(database_pickle_ch) 
    emit:
    json = sample_deepchem.out.json

}