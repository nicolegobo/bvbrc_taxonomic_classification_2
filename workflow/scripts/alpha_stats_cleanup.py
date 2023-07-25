import csv
import sys
#from turtle import width

import pandas as pd
import plotly.graph_objects as go

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
        df = df.reset_index(drop=True)
        # Create the Plotly table
        table = go.Table(
            header=dict(values=list(df.columns)),
            cells=dict(values=[df[col] for col in df.columns], height = 30)
            )
        # Set the table layout
        layout = go.Layout(autosize=True,
            width=1000, 
            height=1000, 
            title='Alpha Statistics<br><sup> All analyses are run at the species level </sup>',
            font=dict(size=14), # Increase the font size (change 14 to your desired size)
            )

        figure = go.Figure(data=[table], 
                           layout=layout)
        figure.add_annotation(text= 
            "<br> \
            <b>Alpha diversity</b> reflects diversity (how similar or how different) the microbes are within a single sample. <br> \
            <br> \
            <b>Shannon diversity</b> index tells you how diverse the species in a given community are. A <b>higher value</b>  <br> \
            indicates a greater number of species and the evenness of their abundance. If only one species was assigned in the <br> \
            sample, the index would be 0. <br> \
            <br> \
            <b>Simpsons diversity index</b> is a measure of diversity which takes into account the number of species present <br> \
            and their relative abundance. If a sample has a lot of species but only a few prominent taxa, the diversity is <br> \
            still less. Here a <b>higher value</b> indicates <b>lower diversity</b>. It will always be from zero to one. <br> \
            <br> \
            <b>Simpsons reciprocal index </b>  Using the reciprocal of the Simpsons diversity index is more intuitive, <br> \
            a <b>higher value</b> indicates <b>higher diversity</b>. <br> \
            <br> \
            <b>Berger-Parker index </b> shows the proportional importance of the most abundant species. A <b>higher value</b> <br> \
            indicates a <b>larger portion</b> of the sample is assigned to the dominant species. This metric assumes a linear <br> \
            distribution. <br> \
            <br> \
            <b> Fisher's index </b> An index of diversity as a logseries distribution. It corrects for the upward bias of the <br> \
            Laspeyres Price Index and the downward bias of the Paasche Price Index by taking the geometric average <br> \
            of the two weighted indices.<br>",
            align='left',
            valign='middle',
            showarrow=False,
            xref='paper',
            yref='paper',
            xanchor='left',
            x=0
            )
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
