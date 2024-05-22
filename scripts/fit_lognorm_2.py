
import math
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import mean_absolute_error


# Data for wealth percentiles for ages 18-25 in the U.S., in 2022
# Source: https://personalfinancedata.com/networth-percentile-calculator/?min_age=18&max_age=25&networth=#results
ranks = [10, 20, 30, 40, 50, 60, 70, 80, 90]
z_scores = [-1.282, -0.842, -0.524, -0.253, 0, 0.253, 0.524, 0.842, 1.282]  # From ranks, based on normal distribution
percentiles = [-11270, -270, 720, 6420, 10700, 18501, 26600, 68800, 206400]

# rightshift all datapoints to make them positive, then take natural log 
percentiles = [p - 2 * percentiles[0] for p in percentiles]  # *1 - see bottom of file 
percentiles = [math.log(p) for p in percentiles]

# Fit a linear regression 
z_scores_np = np.array(z_scores)
percentiles_np = np.array(percentiles)
coefficients = np.polyfit(z_scores_np, percentiles_np, 1)
linear_fit = np.poly1d(coefficients)
fitted_values = linear_fit(z_scores_np)
mae = mean_absolute_error(percentiles_np, fitted_values)

# Plot z-scores (x-axis) vs. log-transformed percentiles (y-axis)
plt.figure(figsize=(10, 6))
plt.plot(z_scores, percentiles, marker='o', linestyle='-', color='b', label='Data')
plt.plot(z_scores, fitted_values, color='r', linestyle='--', label=f'Fit: y = {coefficients[0]:.2f}x + {coefficients[1]:.2f}')

# Add title and labels
plt.title('Z-Scores vs. Log-Transformed Percentiles with Linear Regression')
plt.xlabel('Z-Scores')
plt.ylabel('Log-Transformed Percentiles')
plt.grid(True)
plt.legend()
plt.savefig("results/fit_lognorm_2_plot.png")


# print output 
print("Fit for the underlying normal distribution: ")
print(f"Mean: {percentiles_np[z_scores_np == 0][0]}")
print(f"Standard deviation: {coefficients[0]}")
print(f"Mean absolute error: {mae}")



# *1 - we shift because you can't take log of negative values. However, if we only shift such that the min element 
#      becomes 1, its log will be 0 - much less than the next smaller value. This is a distortion to the otherwise 
#      largely linear relationship. To minimize this distortion, we shift by double the amount, so that the 
#      magnitude of the min value is the same as in the original data 

