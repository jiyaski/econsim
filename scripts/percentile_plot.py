
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import argparse



def read_csv_line(file_path, line_number): 
    df = pd.read_csv(file_path, sep='\t') 
    if line_number < -len(df) or line_number > len(df): 
        raise ValueError("Line number {line_number} is out of range.") 
    data = df.iloc[line_number].tolist() 
    data = [float(i) for i in data] 
    return data


def read_csv_vertical(file_path): 
    df = pd.read_csv(file_path, sep='\t', engine='python') 
    values = df.iloc[:,1].tolist() 
    print(df)
    return [0] + values + [0]  # must add padding for 0th & 100th percentiles 



def calculate_percentiles(data, num_percentiles=100):
    data.sort()
    percentiles = []
    for i in range(num_percentiles + 1):
        index = int(np.ceil((i / num_percentiles) * (len(data) - 1)))
        percentiles.append(data[index])
    
    percentiles[0] = 0    # bc percentiles are only meaningful from 1-99
    percentiles[100] = 0 
    return percentiles


def plot_histogram(result_percentiles, plot_filepath, empirical_percentiles=None):
    plt.figure(figsize=(12, 8))
    bar_width = 0.4 if (empirical_percentiles != None) else 0.8
    x = range(101)
    
    plt.bar(x, result_percentiles, width=bar_width, label='Simulated Data', color='#4444ffff')
    
    if empirical_percentiles is not None:
        x_compare = [i + bar_width for i in x]
        plt.bar(x_compare, empirical_percentiles, width=bar_width, label='Empirical Data', color='#ff5555bb')
    
    plt.xlabel('Percentiles')
    plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(10))
    plt.ylabel('Net Worth')
    plt.yscale('symlog', linthresh=100, linscale=0.5)
    plt.title('Wealth Distribution Percentiles')
    plt.grid(alpha=0.7)
    plt.legend()
    plt.savefig(plot_filepath)
    plt.close()


def main():
    default_results_data_filepath = 'data/results.csv'
    default_plot_filepath = 'results/plot.png'
    default_line_number = -1


    # Check command line arguments
    parser = argparse.ArgumentParser(description="Process some CSV data and plot histograms.")
    parser.add_argument('-d', '--data', type=str, default=default_results_data_filepath, help='Path to the results data CSV file')
    parser.add_argument('-p', '--plot', type=str, default=default_plot_filepath, help='Path to save the plot image')
    parser.add_argument('-l', '--line', type=int, default=default_line_number, help='Line number to read from the CSV file')
    parser.add_argument('-c', '--compare', type=str, help='Path to CSV containing empirical data to plot alongside results')
    args = parser.parse_args()
    
    results_data = read_csv_line(args.data, args.line)
    result_percentiles = calculate_percentiles(results_data)

    if args.compare: 
        empirical_percentiles = read_csv_vertical(args.compare) 
        plot_histogram(result_percentiles, args.plot, empirical_percentiles)
    else: 
        plot_histogram(result_percentiles, args.plot)



if __name__ == "__main__":
    main()

