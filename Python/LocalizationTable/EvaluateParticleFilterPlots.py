import matplotlib.pyplot as plt
from collections import Counter
import numpy as np
import scipy

filename = "errors.txt"

PLOT_PROBABILITIES = True

# Plot giving error probability:
if PLOT_PROBABILITIES:
    # Open file and read error values:
    errors = []
    with open(filename, 'r') as file:
        for line in file:
            error = float(line.strip())
            errors.append(error)

    # Calculating probabilities:
    errors_sorted = np.sort(errors)
    p = 1. * np.arange(len(errors)) / (len(errors) - 1)

    plt.plot(errors_sorted, p)
    plt.show()



