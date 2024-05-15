import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

def plot_last_row_histogram(csv_file, output_image):
    # Read the last row from a CSV file, assuming tab-separated values
    try:
        data = pd.read_csv(csv_file, delimiter='\t')
        last_row = data.iloc[-1]
    except Exception as e:
        print(f"Error reading data: {e}")
        return

    # Ensure all data is numeric
    last_row = pd.to_numeric(last_row, errors='coerce')
    last_row = last_row.dropna()
    last_row_values = last_row.tolist()
    if not last_row_values:
        print("No valid numeric data found in the last row.")
        return

    # Calculate percentiles
    percentiles = np.percentile(last_row_values, [10, 50, 90, 99])

    # Plotting the data
    plt.figure(figsize=(10, 6))
    count, bins, ignored = plt.hist(last_row_values, bins=15, log=True, alpha=0.75, color='blue', edgecolor='black')

    # Marking percentiles
    for percentile, description in zip(percentiles, ["10th", "50th", "90th", "99th"]):
        plt.axvline(x=percentile, color='red', linestyle='dashed', linewidth=1)
        plt.text(percentile, plt.ylim()[1]*0.9, f'{description}: {percentile:.2f}', rotation=45)

    plt.xscale('log')
    plt.xlabel('Data Values')
    plt.ylabel('Frequency')
    plt.title('Histogram of the Last Row with Log-Scale X-axis')

    # Save the figure
    plt.savefig(output_image)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Plot the last row of a TSV file as a histogram.")
    parser.add_argument("csv_file", help="The path to the TSV file")
    parser.add_argument("--output", help="The output image file name", default="last_row_histogram.png")

    args = parser.parse_args()

    plot_last_row_histogram(args.csv_file, args.output)

if __name__ == "__main__":
    main()
