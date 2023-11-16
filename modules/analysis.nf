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
        path("results.*")
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
    result_df = pd.DataFrame(columns=["molecule", "PC+dG*(solv)(kcal/mol)", "expt (solv)(kcal/mol)"])

    # Open and read the JSON files
    for file_id in file_ids:
        with open(file_id, "r") as file:
            data = json.load(file)

        # Extract relevant information
        molecule = data["molecule"]
        pc_dg_solv = data["PC+dG*(solv)(kcal/mol)"]
        expt_solv = y[file_id]

        # Append the data to the result DataFrame
        result_df = pd.concat([result_df, pd.DataFrame({"molecule": [molecule], "PC+dG*(solv)(kcal/mol)": [pc_dg_solv], "expt (solv)(kcal/mol)": [expt_solv]})], ignore_index=True)

    # Save the result DataFrame to a CSV file
    result_df.set_index("molecule", inplace=True)
    result_df.to_csv("results.csv")

    # Calculate the Pearson correlation coefficient
    correlation_coefficient, _ = pearsonr(result_df["PC+dG*(solv)(kcal/mol)"], result_df["expt (solv)(kcal/mol)"])
    print("R=",correlation_coefficient)

    # Create a scatter plot
    plt.figure(figsize=(10, 6))
    plt.scatter(result_df["expt (solv)(kcal/mol)"], result_df["PC+dG*(solv)(kcal/mol)"], label="Data points")

    # Calculate and plot a linear trendline
    trendline = np.polyfit(result_df["expt (solv)(kcal/mol)"], result_df["PC+dG*(solv)(kcal/mol)"], 1)
    plt.plot(result_df["expt (solv)(kcal/mol)"], np.polyval(trendline, result_df["expt (solv)(kcal/mol)"]), color='red', label='Trendline')

    plt.title("Scatter Plot of PC+dG*(solv) vs. expt (R = {pearson})".format(pearson=correlation_coefficient))
    plt.xlabel("expt (ΔGsolv)(kcal/mol)", fontsize=20)
    plt.ylabel("PC+dG*(ΔGsolv)(kcal/mol)", fontsize=20)
    plt.legend()
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