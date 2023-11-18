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
        path("results.*")
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
    database = pd.read_pickle(path_to_database)

    # Extract 'expt' from inner dictionaries and create a DataFrame with keys as indices
    result_df = pd.DataFrame({'expt': [inner_dict['expt'] for inner_dict in database.values()]}, index=database.keys()).rename_axis('molecule')

    # Assuming the list of file paths is provided
    file_ids = "${ids}".split(" ")

    # Open and read the JSON files
    for file_id in file_ids:
        with open(file_id, "r") as file:
            data = json.load(file)

        # Extract relevant information
        molecule = data["molecule"]
        solvent = data["solvent"]
        pc_dg_solv = data["PC+dG*(solv)(kcal/mol)"]

        # Append the data to the result DataFrame
        # result_df = pd.concat([result_df, pd.DataFrame({"molecule": [molecule], "PC+dG*(solv)(kcal/mol)": [pc_dg_solv], "expt (solv)(kcal/mol)": [expt_solv]})], ignore_index=True)
        result_df.at[molecule, solvent] = pc_dg_solv

    # Save the result DataFrame to a CSV file
    result_df = result_df.dropna()

    #result_df.set_index("molecule", inplace=True)
    result_df.to_csv("results.csv")



    # Create a scatter plot
    plt.figure(figsize=(10, 6))
    # Plotting
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(result_df['expt'], result_df['expt'], label="Reference")

    models = []
    r2_values = []
    mad_values = []
    for col in result_df.columns[1:]:
        # Calculate the Pearson correlation coefficient
        correlation_coefficient, _ = pearsonr(result_df[col], result_df["expt"])
        print(col,"R=",correlation_coefficient, _)
        # Calculate Mean Absolute Deviation
        mad = np.mean(np.abs(result_df[col]-result_df["expt"]))
        models.append(col)
        r2_values.append(correlation_coefficient)
        mad_values.append(mad)

        ax.scatter(result_df['expt'], result_df[col], label=col)
        # Calculate and plot a linear trendline
        trendline = np.polyfit(result_df["expt"], result_df[col], 1)
        plt.plot(result_df["expt"], np.polyval(trendline, result_df["expt"]), label=col)

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