import numpy as np
import matplotlib.pyplot as plt
import math as m
from astropy.visualization import hist
import pandas as pd
from astropy.stats import bayesian_blocks

np.random.seed(42)

# Define the function with fluctuations and plateaus
# def f(x):
#     y = -x + 1
#     y += np.random.normal(0.5, 100, len(x)) # Add Gaussian noise with mean 0 and standard deviation 30
#     plateau_x_values = [0.2, 0.5, 0.9]
#     plateau_heights = [0.1, 1.5, 0.2]
#     for x_val, height in zip(plateau_x_values, plateau_heights):
#         y[np.abs(x - x_val) < 0.05] += height
#     return y

# Generate random numbers in the interval [0, 1]
# num_samples = 10000
# samples = np.random.uniform(0, 1, num_samples)
#
# # Apply the modified function to the samples
# sampled_values = f(samples)
#
# bins = np.linspace(0, 1, int(m.sqrt(num_samples)/2))

sampled_values = 100 * np.random.random(100)
x = np.exp(-0.5 * (sampled_values - 50) ** 2)
sigma = 0.1
x_obs = np.random.normal(x, sigma)
edges = bayesian_blocks(sampled_values, x_obs, sigma, fitness='measures')

val, hbins, _ = hist(sampled_values, bins=25, color='r', density=False, label='Cumulative', histtype='stepfilled', alpha=0.1, cumulative=False)
# plt.close()
bin_centers = (hbins[:-1] + hbins[1:]) / 2

# plt.scatter(bin_centers, val, color='b', marker='.')

data = np.array([])

for i in range(len(bin_centers)-1):
    tp = np.ones(int(val[i]))*bin_centers[i]
    data = np.append(data, tp)

plt.hist(data, bins=bin_centers, color='tab:blue', density=False, alpha=0.3, label='data')

valb, binb, _ = hist(sampled_values, bins='blocks', color='black', histtype='step', density=True, label='bb', alpha=0.8)

plt.xlabel('f(x)')
plt.ylabel('Cumulative Count')
plt.title('Cumulative Distribution with Bayesian Blocks')
plt.legend()

# Invert the x-axis
plt.gca().invert_xaxis()
# plt.ylim(0, 2.5)

# print(valb)
# print(binb)

# Assuming val, hbins, and data are arrays of different lengths
# Fill val and hbins with zeros to match the length of data
# max_length = max(len(val), len(hbins), len(data))
# val_n = np.pad(val, (0, max_length - len(val)), mode='constant', constant_values=np.nan)
# hbins_n = np.pad(hbins, (0, max_length - len(hbins)), mode='constant', constant_values=np.nan)
#
# # Create DataFrame
# dataframe = pd.DataFrame({'val': val_n, 'hbins': hbins_n, 'data': data})
#
# # Save DataFrame to CSV
# dataframe.to_csv('./data_bbh.csv', sep=',', index=False)

plt.show()