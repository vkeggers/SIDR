"""
Created on Mon 20 May, 2024 at 17:04:09
Updated on Mon 20 May, 2024 at 18:44:17
@author: vkeggers and virallyDanny
@description:

now = datetime.now() # current date and time
date_time = now.strftime("%a %d %B, %Y at %H:%M:%S")
print("date and time:",date_time)

@Usage: `python xgboost.py --Genus <Genus>`
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def plot_origin_counts(input_csv, output_svg):
    # Usage: plot_origin_counts('path/to/inputDataTable', 'path/to/outputSVG')
    # Import the dataframe for the unprocessed data
    preprocessed_data = pd.read_csv(input_csv, sep='\t', header=0)

    # Count the number of unique values in preprocessed_data under the Origin column as a dictionary
    origin_counts = preprocessed_data["Origin"].value_counts().reset_index()
    origin_counts.columns = ["Origin", "Count"]

    # Export the results to a new dataframe
    BLASTresultsDF = pd.DataFrame(origin_counts)

    # Save the unique values of the Origin column to a new list titled 'BLASTresults'
    BLASTresults = preprocessed_data['Origin'].unique()

    # Plotting the count of each Origin using seaborn
    plt.figure(figsize=(20, 6))
    barplot = sns.barplot(x='Origin', y='Count', data=BLASTresultsDF)
    barplot.set_title('Count of each Origin')
    barplot.set_xlabel('Origin')
    barplot.set_ylabel('Count')
    plt.xticks(rotation=45)
    plt.tight_layout()

    # Save the plot
    plt.savefig(output_svg, format='svg')
    plt.show()