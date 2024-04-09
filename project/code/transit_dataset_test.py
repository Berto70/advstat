import numpy as np
import matplotlib.pyplot as plt
import astropy.visualization as vis

# s1 = np.random.normal(loc=4, scale=0.8, size=1000)
# s2 = np.random.normal(loc=7, scale=0.6, size=1000)
# s3 = np.random.normal(loc=2, scale=0.01, size=1000)
# noise = np.random.uniform(0, 8, 10_000)
#
# data = np.concatenate((s1, s2, s3, noise), axis=0)
#
# # Plot histogram
# plt.hist(data, bins=50, alpha=0.7, label='Data', density=True)
# vis.hist(data, bins='blocks', alpha=0.3, label='bb', color='red', density=True, histtype='step')
#
# plt.legend(loc='upper right')
# plt.ylim(0, 0.8)
# plt.show()

def generate_fake_transit_data_counts(duration=20, depth=0.03, period=10.0, num_points=1000, noise_level=0.8):
    time = np.linspace(0, period, num_points)
    counts = np.ones_like(time) * 1000  # Starting counts

    # Simulate transit
    ingress_duration = 0.1 * period  # Duration of ingress/egress
    ingress_mask = np.logical_and(time >= period / 2 - ingress_duration / 2, time <= period / 2 + ingress_duration / 2)
    egress_mask = np.logical_and(time >= period / 2 + duration - ingress_duration / 2,
                                 time <= period / 2 + duration + ingress_duration / 2)
    transit_mask = np.logical_or(ingress_mask, egress_mask)
    counts[transit_mask] -= depth * 1000  # Reduce counts during transit

    # Add noise
    noise = np.random.poisson(lam=noise_level*1000, size=num_points)  # Poisson noise for counts
    counts += noise

    return time, counts


# Generate fake transit data as counts
time_counts, counts = generate_fake_transit_data_counts()

# Reshape time_counts and counts to have shape (n, 1)
# time_counts, counts = time_counts.reshape(len(time_counts), 1), counts.reshape(len(counts), 1)

new_data = np.zeros(shape=(int(np.sum(counts)), 1))

k = 0
for i in range(counts.shape[0]):
    l = counts[i]  # Extract scalar value from counts[i]
    if l != 0:
        repeated_values = np.repeat(time_counts[i], l)[:, np.newaxis]  # Add new axis
        new_data[k:k+int(l)] = repeated_values
        k += int(l)


# Plot the histogram
plt.figure(figsize=(10, 5))
plt.hist(new_data, bins=50, color='skyblue', edgecolor='black', alpha=0.7, density=True)
vis.hist(new_data, bins='blocks', histtype='step', color='red', density=True)
plt.xlabel('Time')
plt.ylabel('Counts')
plt.title('Histogram of Simulated Exoplanet Transit Light Curve (Counts)')
plt.ylim(0.090, 0.105)
plt.show()

# Plot the light curve
plt.figure(figsize=(10, 5))
plt.plot(time_counts, counts, color='black', lw=0.5)
plt.xlabel('Time')
plt.ylabel('Counts')
plt.title('Simulated Exoplanet Transit Light Curve (Counts)')
plt.grid(True)
plt.show()

