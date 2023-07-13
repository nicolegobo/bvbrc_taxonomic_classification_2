import pandas as pd
import plotly.graph_objects as go
import sys


def alpha_stats_viz(input_csvs, output_csv, output_html):
    main_df = pd.DataFrame()
    dfs = []
    for input_csv in input_csvs:
        df = pd.read_csv(input_csv, sep=':', index_col="Alpha Diversity Tests")
        # Set the number of significant figures
        significant_figures = 3
        # Round the DataFrame values to the desired significant figures
        df = df.round(significant_figures)
        # Append the dataframe to the list
        main_df = main_df.merge(df, how = 'outer', right_index = True, left_index = True)

    # export the dataframe to CSV with the Alpha Statistics column as the index
    main_df.to_csv(output_csv, index_label= "Alpha Diversity Tests")
    main_df = main_df.reset_index()

    # Create the Plotly table
    table = go.Table(
        header=dict(values=list(main_df.columns)),
        cells=dict(values=[main_df[col] for col in main_df.columns])
    )
    # Set the table layout
    layout = go.Layout( autosize=True,
                        margin=dict(l=40, r=0, t=50, b=0),
                        width=1000, 
                        height=1000, 
                        title='Alpha Statistics',
                        font=dict(size=14))  # Increase the font size
    figure = go.Figure(data=[table], 
                        layout=layout)
    figure.add_annotation(text= 
            "<br> \
            <b>Alpha diversity</b> reflects diversity (how similar or how different) the microbes are within a single sample. <br> \
            All analyses are run at the species level. <br> \
            <br> \
            <b>Shannon diversity</b> index tells you how diverse the species in a given community are. A<b>higher value</b> \
            indicates a greater number of species and the evenness of their abundance. If only one species was assigned in the  <br> \
            sample, the index would be 0. <br> \
            <br> \
            <b>Simpsons diversity index</b> is a  measure of diversity which takes into account the number of species present <br> \
            and their relative abundance. If a sample has a lot of species but only a few prominent taxa, the diversity is <br> \
            still less. Here a <b>higher value</b> indicates <b>lower diversity</b>. It will always be from zero to one. <br> \
            <br> \
            <b>Simpsons reciprocal index</b>  Using the reciprocal of the Simpsons diversity index is more intuitive, <br> \
            a <b>higher value</b> indicates <b>higher diversity</b>. <br> \
            <br> \
            <b>Berger-Parker index</b> shows the proportional importance of the most abundant species. A <b>higher value</b> <br> \
            indicates a <b>larger portion</b> of the sample is assigned to the dominant species. This metric assumes a linear <br> \
            distribution. <br> \
            <br> \
            <b>Fisher's index</b> An index of diversity as a logseries distribution. It corrects for the upward bias of the <br> \
            Laspeyres Price Index and the downward bias of the Paasche Price Index by taking the geometric average <br> \
            of the two weighted indices.<br>",
            align='left',
            valign='middle',
            showarrow=False,
            xref='paper',
            yref='paper',
            xanchor='left',
            x=0,
            y=0.35
        )
    # Save the figure as a self-contained HTML file
    figure.write_html(output_html)
    return

if __name__ == "__main__":
    # The first argument is the script name, so we ignore it
    arguments = sys.argv[1:]
    # The last argument is the HTML output path
    output_html = arguments[-1]
    # the second to last argument is the CSV path
    output_csv = arguments[-2]
    # The remaining arguments are input files
    input_csvs = arguments[:-2]
    alpha_stats_viz(input_csvs, output_csv, output_html)