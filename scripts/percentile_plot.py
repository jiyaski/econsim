import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import sys


def read_data(file_path, line_number):
    # Read the specified line or last line of data from the tab-delimited CSV file
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        if line_number == -1:
            data = None
            for data in reader:
                pass  # This will end with 'data' being the last line
            if data is None:
                raise ValueError("The file is empty.")
        else:
            for current_line, data in enumerate(reader):
                if current_line == line_number:
                    break
            else:
                raise ValueError(f"Line number {line_number} is out of range.")
        
        data = [float(i) for i in data]  # Convert to float
    return data


def calculate_percentiles(data, num_percentiles=100):
    data.sort()
    percentiles = []
    for i in range(num_percentiles + 1):
        index = int(np.ceil((i / num_percentiles) * (len(data) - 1)))
        percentiles.append(data[index])
    
    percentiles[0] = 0    # bc percentiles are only meaningful from 1-99
    percentiles[100] = 0 
    return percentiles


def plot_histogram(percentiles, plot_filepath):
    plt.figure(figsize=(10, 6))
    plt.bar(range(101), percentiles, width=1.0)
    plt.xlabel('Percentiles') 
    plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(10))
    plt.ylabel('Net Worth')
    plt.yscale('symlog', linthresh=100, linscale=0.5) 
    plt.title('Wealth Distribution Percentiles') 
    plt.grid(alpha=0.7) 
    plt.savefig(plot_filepath)
    plt.close()


def main():
    default_data_filepath = 'data/data.csv'
    default_plot_filepath = 'results/plot.png'
    default_line_number = -1

    # Check command line arguments
    data_filepath = sys.argv[1] if len(sys.argv) > 1 else default_data_filepath
    plot_filepath = sys.argv[2] if len(sys.argv) > 2 else default_plot_filepath
    line_number = int(sys.argv[3]) if len(sys.argv) > 3 else default_line_number

    data = read_data(data_filepath, line_number)
    percentiles = calculate_percentiles(data)
    plot_histogram(percentiles, plot_filepath)


if __name__ == "__main__":
    main()

