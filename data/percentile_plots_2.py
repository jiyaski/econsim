import csv
import numpy as np
import matplotlib.pyplot as plt
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


def calculate_percentiles(data, num_bins=100):
    data.sort()
    bin_size = len(data) // num_bins
    bins = [data[i*bin_size:(i+1)*bin_size] for i in range(num_bins)]
    bin_averages = [np.mean(bin) for bin in bins]
    return bin_averages


def plot_histogram(bin_averages):
    plt.figure(figsize=(10, 6))
    plt.bar(range(1, 101), bin_averages, width=1.0)
    plt.xlabel('Percentiles')
    plt.ylabel('Average Value')
    plt.yscale('log') 
    plt.title('Histogram of Percentile Averages')
    plt.savefig('temp_plot.png')
    plt.close()


def main():
    default_file_path = 'data.csv'
    default_line_number = -1

    # Check command line arguments
    file_path = sys.argv[1] if len(sys.argv) > 1 else default_file_path
    line_number = int(sys.argv[2]) if len(sys.argv) > 2 else default_line_number

    data = read_data(file_path, line_number)
    bin_averages = calculate_percentiles(data)
    plot_histogram(bin_averages)



if __name__ == "__main__":
    main()
