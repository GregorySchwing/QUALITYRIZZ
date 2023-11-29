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
    publishDir "${params.output_folder}/${params.database}/results", mode: 'copy', overwrite: true

    debug true
    input:
        path(results)
        path(database)
    output:
        path("merged_results.csv")
        path("stats.csv")
        path("results.png")
    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.stats import pearsonr  # Import the pearsonr function
    # Read the CSV file into a DataFrame
    df_database = pd.read_csv("${database}")
    print(df_database.head())
    # Read the CSV file into a DataFrame
    df_results = pd.read_csv("${results}")
    print(df_results.head())
    # Assuming your first dataframe is df1 and the second one is df2
    result_df2 = pd.merge(df_database, df_results, on='molecule')
    result_df2.to_csv("merged_results.csv")

    # Create a scatter plot
    plt.figure(figsize=(10, 6))
    # Plotting
    print("${params.reference_col}")
    print(result_df2["${params.reference_col}"])

    # Plotting
    fig, ax = plt.subplots(figsize=(10, 6))
    import seaborn as sns
    # Define markers for each charge group
    markers = {'am1bcc': 'o', 'gasteiger': 's', 'RESP': '^'}

    # Loop through each combination of 'partial_charge_method' and 'solvent'
    for (charge_method, solvent), group in result_df2.groupby(['partial_charge_method', 'solvent']):
        marker = markers.get(charge_method, 'o')  # Default to 'o' if charge_method is not found
        sns.regplot(
            x=group['${params.reference_col}'],
            y=group['PC+dG*(solv)(kcal/mol)'],
            label=f'{charge_method} - {solvent}',
            marker=marker,
            ax=ax
        )
    plt.plot(result_df2["${params.reference_col}"], result_df2["${params.reference_col}"], color='black', linewidth=2, linestyle='--', label='Reference')

    # Add labels and legend
    ax.set_xlabel('Experimental Value (kcal/mol)')
    ax.set_ylabel('PC+dG*(solv) (kcal/mol)')

    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    # Save the plot as a PNG file
    fig.savefig('results.png', bbox_inches='tight')
    # Calculate the mean absolute deviation
    result_df2['absolute_deviation'] = abs(result_df2['PC+dG*(solv)(kcal/mol)'] - result_df2["${params.reference_col}"])

    # Filter by the desired 'partial_charge_method' values
    selected_charge_methods = ['am1bcc', 'am1-mulliken', 'gasteiger', 'RESP']
    filtered_df = result_df2[result_df2['partial_charge_method'].isin(selected_charge_methods)]

    # Calculate the mean absolute deviation for each 'partial_charge_method'
    grouped_stats = filtered_df.groupby(['partial_charge_method', 'solvent'])['absolute_deviation'].mean()
    grouped_stats_std = filtered_df.groupby(['partial_charge_method', 'solvent'])['absolute_deviation'].std()

    # Calculate correlation coefficients
    correlation_coefficients = (
        result_df2.groupby(['partial_charge_method', 'solvent'])
        .apply(lambda group: group['${params.reference_col}'].corr(group['PC+dG*(solv)(kcal/mol)']))
    )

    # Combine both into a single DataFrame
    stats = pd.DataFrame({'MAD': grouped_stats, 'MAD_STDEV':grouped_stats_std, 'Correlation_Coefficients': correlation_coefficients})
    print(stats)

    # Write to a single CSV file
    stats.to_csv('stats.csv', header=True, index=True, mode='w', sep=',')

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