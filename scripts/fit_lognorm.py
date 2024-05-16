import sys
import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import os

# input CSV is given as CLI arg 
if len(sys.argv) < 2:
    print("Usage: python3 fit_lognorm.py /path/to/your_data.csv")
    sys.exit(1) 

file_path = sys.argv[1]
base_name = os.path.splitext(os.path.basename(file_path))[0]

data = pd.read_csv(file_path, sep='\t')
percentile_column = data.columns[0]
data_column = data.columns[1]

# Fit the log-normal distribution to the salary data
shape, loc, scale = stats.lognorm.fit(data[data_column], floc=0)

# The scale parameter is e^mu in the log-normal distribution
mu = np.log(scale)
# The shape parameter is the sigma in the log-normal distribution
sigma = shape

print(f"Estimated mu: {mu}")
print(f"Estimated sigma: {sigma}")

# Plot the CDF with the estimated parameters
percentile_ranks = data[percentile_column] / 100
x = np.linspace(min(data[data_column]), max(data[data_column]), 1000)
cdf = stats.lognorm.cdf(x, sigma, scale=scale)

# Plot the actual salaries against the estimated CDF
plt.scatter(data[data_column], percentile_ranks, label="Data", alpha=0.5)
plt.plot(x, cdf, label="Log-normal fit")
plt.xlabel(data_column)
plt.ylabel('Cumulative Probability')
plt.title('Log-normal Fit to ' + data_column + ' Data')
plt.legend()

# Save the plot with a filename based on the input CSV filename
plot_file_name = base_name + '_plot.png' 
plt.savefig(plot_file_name) 
print(f"Plot saved to {plot_file_name}") 

