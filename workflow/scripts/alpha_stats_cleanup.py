import csv
import sys

import pandas as pd
#import plotly.express as px
import plotly.graph_objects as go
#import plotly.offline as offline

def edit_alpha_stats(input_csv, output_csv, output_html):
    with open(input_csv, 'r') as infile:
        # Open the output CSV file
        with open(output_csv, 'w') as outfile:
            # Create a CSV reader and writer objects
            reader = csv.reader(infile)
            writer = csv.writer(outfile)

            # Iterate through the rows of the input CSV file
            for row in reader:
                # Check if the row contains the line to be removed
                if "Fisher's alpha...loading" in row:
                    # Skip the row if it contains the line to be removed
                    continue

                # Otherwise, write the row to the output CSV file
                writer.writerow(row)
        df = pd.read_csv(output_csv, sep=":")
        # Set the number of significant figures
        significant_figures = 3
        # Round the DataFrame values to the desired significant figures
        df = df.round(significant_figures)
        #df = df.reset_index()
        df = df.reset_index(drop=True)
        # Create the Plotly table
        table = go.Table(
            header=dict(values=list(df.columns)),
            cells=dict(values=[df[col] for col in df.columns]))
        # Set the table layout
        layout = go.Layout(width=1000, 
                           height=1000, 
                           title='Alpha Statistics',
                            font=dict(size=14)  # Increase the font size (change 14 to your desired size)
)
        figure = go.Figure(data=[table], 
                           layout=layout)
        # Save the figure as a self-contained HTML file
        figure.write_html(output_html)
    return 
    
# run the script from service-script/app_taxonomic_classification 
def main(argv):
    inputfile = argv
    input_csv = inputfile[0]
    output_csv = inputfile[1]
    output_html = inputfile[2]
    edit_alpha_stats(input_csv, output_csv, output_html)

if __name__ == "__main__":
    main(sys.argv[1:])
