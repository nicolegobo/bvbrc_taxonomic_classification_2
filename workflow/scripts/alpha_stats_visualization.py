import pandas as pd
import plotly.graph_objects as go
import sys


def alpha_stats_viz(input_csvs, output_csv, output_html):
    main_df = pd.DataFrame()
    dfs = []
    for input_csv in input_csvs:
        # df = pd.read_csv(input_csv, sep=':', index_col= 0)
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
    layout = go.Layout(width=600, 
                        height=600, 
                        title='Alpha Statistics',
                        font=dict(size=14))  # Increase the font size
    figure = go.Figure(data=[table], 
                        layout=layout)
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