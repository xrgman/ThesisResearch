import numpy as np
import matplotlib.pyplot as plt

# Sample data
data = [(1, 'A', 10), (2, 'B', 15), (1, 'C', 12), (2, 'D', 8), (3, 'E', 20)]

x_values = [x for x in range(1, 4)]

bar_width = 0.1
total_group_width = bar_width * len(x_values)

# Set the width of the plot for nicer results:
plt.figure(figsize=(13, 6), dpi=300)

# Calculate the x-coordinates for each group of bars
x_positions = np.arange(len(x_values)) - (total_group_width / 2) + (bar_width / 2)

# Plot the bar charts for each group
for i in range(1, 4):
    values = [x[2] for x in data if x[0] == i+1]


    plt.bar(x_positions + i * bar_width, values, width=bar_width, label=f'{symbol_bits_values[i]} samples')


plt.xticks(np.arange(len(x_values)), x_values)

plt.rcParams['text.antialiased'] = True
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.tight_layout = True


plt.xlabel("SNR (dB)")
plt.ylabel("Correctly decoded messages (%)")
plt.title("Correct message decoding noise vs number of samples")
plt.legend()