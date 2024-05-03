import re
import sys

import pandas as pd
import plotly.express as px
import plotly.offline as offline

def edit_beta_stats(input_csv, output_csv, output_html):
    table_start = int
    with open(input_csv, 'r') as infile:
        my_lines = [line.strip() for line in infile]
        sample_dict = {}
        for line in my_lines:
            sample_id = ""
            # set up dictionary
            if line.startswith("#"):    
                parts = line.split()
                sample_id = parts[1].rstrip('_')
                sample_id = sample_id.replace("_bracken_output.txt", "")
                parts=sample_id.split("bracken_output/")
                sample_id = parts[1]
                num_pattern = r"#(\d+)"
                num_match = re.search(num_pattern, line)
                numeric_value = int(num_match.group(1))
                sample_dict[numeric_value] = sample_id
            if line.startswith("x"):
                table_start = my_lines.index(line)
    df = pd.read_csv(input_csv, skiprows=table_start, sep='\t', index_col= 0)

    df.columns = df.columns.astype(int)
    #rename columns in DataFrame sample_dict
    df = df.rename(columns = sample_dict)

    # Rename index using the sample_dict
    df = df.rename(index = sample_dict)

    # add description
    description = 'Bray-Curtis dissimilarity [1] examines the abundances of microbes that are shared between two samples, and the number of microbes found in each. Bray-Curtis dissimilarity ranges from 0 to 1. If both samples share the same number of microbes, at the same abundance, their "dissimilarity" will equal zero. If they have absolutely no shared microbes, they will have maximum dissimilarity, i.e. 1. \n'
    with open(output_csv, 'w') as out_csv:
        out_csv.write(description)
        df.to_csv(out_csv, index=True)
    
    # make the heatmap
    fig = px.imshow(df, text_auto=True,
                    labels=dict(x="Sample IDs", y="Sample IDs", color="Bray-Curtis Index"),
                    title ='Beta Diversity (Bray-Curtis Index of Dissimilarity) <br><sup> If two samples have the same microbes at the same abundance the dissimilarity is 0. Conversely, if there are no shared microbes between the dissimilarity is 1. </sup>'
                    )

    fig.update_yaxes(showgrid=True, tickangle=45)
    fig.update_xaxes(showgrid=True, tickangle=45)
    offline.plot(fig, filename=output_html, auto_open=False)

# run the script from service-script/app_taxonomic_classification 
def main(argv):
    inputfile = argv
    input_csv = inputfile[0]
    output_csv = inputfile[1]
    output_html = inputfile[2]
    edit_beta_stats(input_csv, output_csv, output_html)

if __name__ == "__main__":
    main(sys.argv[1:])
