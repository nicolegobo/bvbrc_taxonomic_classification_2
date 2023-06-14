import csv
import sys

def edit_alpha_stats(input_csv, output_csv):
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
    return 
    
# run the script from service-script/app_taxonomic_classification 
def main(argv):
    inputfile = argv
    input_csv = inputfile[0]
    output_csv = inputfile[1]
    edit_alpha_stats(input_csv, output_csv)

if __name__ == "__main__":
    main(sys.argv[1:])
