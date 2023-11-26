#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process analyze {
    container "${params.container__deep_chem}"
    publishDir "${params.output_folder}/${params.database}/results", mode: 'copy', overwrite: false

    debug true
    input:
        path(ids)
    output:
        path("*")
    script:
    """
        #!/usr/bin/env python
        import json
        from deepchem.molnet import load_freesolv

        import pandas as pd
        import matplotlib.pyplot as plt
        import numpy as np
        from scipy.stats import pearsonr  # Import the pearsonr function


        # Assuming the list of file paths is provided
        file_ids = "${ids}".split(" ")

        tasks, datasets, transformers = load_freesolv(splitter=None)
        train = datasets[0]
        x,y,w,ids = train.X, train.y, train.w, train.ids

        # Create an empty DataFrame to store the extracted data
        #"PC+dG*(solv)(kcal/mol)", 

        # Create an empty DataFrame to store the extracted data
        result_df = pd.DataFrame(columns=["expt (solv)(kcal/mol)"])
        result_df["expt (solv)(kcal/mol)"] = y.flatten()
        result_df.index.names = ['molecule']
        print(result_df)

        # Open and read the JSON files
        for file_id in file_ids:
            with open(file_id, "r") as file:
                data = json.load(file)

                # Extract relevant information
                molecule = data["molecule"]
                solvent = data["solvent"]
                pc_dg_solv = data["PC+dG*(solv)(kcal/mol)"]
                result_df.at[int(molecule), solvent] = pc_dg_solv

        result_df = result_df.dropna()
        print(result_df)

        # Save the result DataFrame to a CSV file
        result_df.to_csv("results.csv")



        # Create a scatter plot
        plt.figure(figsize=(10, 6))
        # Plotting
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(result_df['expt (solv)(kcal/mol)'], result_df['expt (solv)(kcal/mol)'], label="Reference")

        models = []
        r2_values = []
        mad_values = []
        for col in result_df.columns[1:]:
            # Calculate the Pearson correlation coefficient
            correlation_coefficient, _ = pearsonr(result_df[col], result_df["expt (solv)(kcal/mol)"])
            print(col,"R=",correlation_coefficient, _)
            # Calculate Mean Absolute Deviation
            mad = np.mean(np.abs(result_df[col]-result_df["expt (solv)(kcal/mol)"]))
            models.append(col)
            r2_values.append(correlation_coefficient)
            mad_values.append(mad)

            ax.scatter(result_df['expt (solv)(kcal/mol)'], result_df[col], label=col)
            # Calculate and plot a linear trendline
            trendline = np.polyfit(result_df["expt (solv)(kcal/mol)"], result_df[col], 1)
            plt.plot(result_df["expt (solv)(kcal/mol)"], np.polyval(trendline, result_df["expt (solv)(kcal/mol)"]), label=col)

        # Creating DataFrame
        data = {'R2': r2_values, 'MAD' : mad_values}
        stats = pd.DataFrame(data, index=models)

        # Displaying the DataFrame
        print(stats)
        stats.to_csv("stats.csv")


        plt.xlabel("expt (ΔGsolv)(kcal/mol)", fontsize=20)
        plt.ylabel("PC+dG*(ΔGsolv)(kcal/mol)", fontsize=20)
        ax.legend()
        plt.grid(True)

        # Save the plot as a PNG file
        plt.savefig("results.png")
        """
}


process analyze_mobley {
    container "${params.container__openff_toolkit}"
    publishDir "${params.output_folder}/${params.database}/results", mode: 'copy', overwrite: false

    debug true
    input:
        path(ids)
        path(pathToDatabase)
    output:
        path("*")
    script:
    """
    #!/usr/bin/env python
    import json
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.stats import pearsonr  # Import the pearsonr function

    # Assuming the path to the database is specified
    path_to_database = "${pathToDatabase}"
    json_data_list = []
    # Assuming the list of file paths is provided
    file_ids = "${pathToDatabase}".split(" ")
    for file_path in file_ids:
        # Read JSON file
        with open(file_path, 'r') as json_file:
            json_data = json.load(json_file)
            json_data_list.append(json_data)

    print(json_data_list)
    # Create an empty DataFrame
    database = pd.DataFrame()

    # Loop through each dictionary in the list
    for data_dict in json_data_list:
        # Extract the key (index) and value (column data)
        index, column_data = next(iter(data_dict.items()))

        # Create a DataFrame from the column data with the index as the first column
        temp_df = pd.DataFrame(column_data, index=[index])

        # Append the temporary DataFrame to the main DataFrame
        database = pd.concat([database, temp_df])

    database = database.rename_axis('molecule')
    # Drop the 'SMILES' column
    database = database.drop('SMILES', axis=1)

    print(database)
    result_df = database
    result_df2 = pd.DataFrame()
    # Assuming the list of file paths is provided
    file_ids = "${ids}".split(" ")

    # Open and read the JSON files
    for file_id in file_ids:
        with open(file_id, "r") as file:
            data = json.load(file)

        # Extract relevant information
        molecule = data["molecule"]
        solvent = data["solvent"]
        partial_charge_method = data["partial_charge_method"]
        if (solvent == 'tip3p_fb-1.1.1.offxml'):
            solvent = 'tip3p-fb-1.1.1.offxml'
        pc_dg_solv = data["PC+dG*(solv)(kcal/mol)"]
        expt_solv = result_df.at[molecule, '${params.reference_col}']
        # Append the data to the result DataFrame
        # result_df = pd.concat([result_df, pd.DataFrame({"molecule": [molecule], "PC+dG*(solv)(kcal/mol)": [pc_dg_solv], "expt (solv)(kcal/mol)": [expt_solv]})], ignore_index=True)
        #result_df.at[molecule, solvent] = pc_dg_solv
        result_df2 = pd.concat([result_df2, pd.DataFrame({"molecule": [molecule], "solvent":[solvent], "partial_charge_method":[partial_charge_method], "PC+dG*(solv)(kcal/mol)": [pc_dg_solv], "${params.reference_col}": [float(expt_solv)]})], ignore_index=True)

    print(result_df2.to_string())
    # Save the result DataFrame to a CSV file
    result_df2 = result_df2.dropna()
    # Convert all columns to numeric
    #result_df2 = result_df2.apply(pd.to_numeric, errors='coerce')

    normalize = False
    if (normalize):
        result_df2 = (result_df2 - result_df2.min()) / (result_df2.max() - result_df2.min())

    result_df2.to_csv("results.csv")


    # Create a scatter plot
    plt.figure(figsize=(10, 6))
    # Plotting
    print("${params.reference_col}")
    print(result_df["${params.reference_col}"])
    print(result_df2["${params.reference_col}"])
    plt.plot(result_df2["${params.reference_col}"], result_df2["${params.reference_col}"], label="Reference", color="black",linewidth=2.0)

    models = []
    r2_values = []
    mad_values = []
    mrd_values = []

    #markers = ['+', 'x', '.', '1']
    filtered_columns = [col for col in result_df.columns if "${params.reference_col}" not in col]
    for col in filtered_columns:
        # Calculate the Pearson correlation coefficient
        correlation_coefficient, _ = pearsonr(result_df[col], result_df["${params.reference_col}"])
        print(col,"R=",correlation_coefficient, _)
        # Calculate Mean Absolute Deviation
        mad = np.mean(np.abs(result_df[col]-result_df["${params.reference_col}"]))
        models.append(col)
        r2_values.append(correlation_coefficient)
        mad_values.append(mad)
        # Calculate Mean Absolute Deviation
        print(result_df["${params.reference_col}"], result_df[col])
        plt.scatter(result_df["${params.reference_col}"], result_df[col], label=col)
        # Calculate and plot a linear trendline
        trendline = np.polyfit(result_df["${params.reference_col}"], result_df[col], 1)
        plt.plot(result_df["${params.reference_col}"], np.polyval(trendline, result_df["${params.reference_col}"]), label=col)
    # Calculate the mean absolute deviation
    result_df2['absolute_deviation'] = abs(result_df2['PC+dG*(solv)(kcal/mol)'] - result_df2["${params.reference_col}"])

    # Filter by the desired 'partial_charge_method' values
    selected_charge_methods = ['am1bcc', 'gasteiger', 'RESP']
    filtered_df = result_df2[result_df2['partial_charge_method'].isin(selected_charge_methods)]

    # Calculate the mean absolute deviation for each 'partial_charge_method'
    mad_by_charge = filtered_df.groupby('partial_charge_method')['absolute_deviation'].mean()

    print(mad_by_charge)
    mad_by_charge.to_csv("stats.csv", header=True, index=True)


    plt.xlabel("${params.reference_col}", fontsize=20)
    plt.ylabel("PC+dG*(ΔGsolv)(kcal/mol)", fontsize=20)
    plt.legend()
    plt.grid(True)

    # Save the plot as a PNG file
    plt.savefig("results.png")
    """
}


process analyze_list_proc {
    container "${params.container__openff_toolkit}"
    publishDir "${params.output_folder}/${params.database}/results", mode: 'copy', overwrite: false

    debug true
    input:
        path(results)
        path(database)
    output:
        //path("*")
    script:
    """
    #!/usr/bin/env python
    import json
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.stats import pearsonr  # Import the pearsonr function

    result_dict = {}
    current_key = None
    
    print($database)
    return
    # Convert the JSON-formatted string to a Python dictionary
    result_dict = json.loads("$database")

    # Print the resulting dictionary
    print(result_dict)


    """
}


workflow analyze_solvation {
    take:
    solvation_results
    main:
    analyze(solvation_results)
}

workflow analyze_mobley_solvation {
    take:
    solvation_results
    database
    main:
    analyze_mobley(solvation_results,database)
}

workflow analyze_list {
    take:
    solvation_results
    database
    main:
    analyze_mobley(solvation_results,database)
}