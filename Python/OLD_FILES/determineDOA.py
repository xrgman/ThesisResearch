import numpy as np
from scipy.signal import oaconvolve

num_channels = 6
fs = 22050

# arrival_times = [9198, 9202, 9197, 9201, 9190, 9193]
#arrival_times = [4422, 4415, 4417, 4420, 4425, 4428]  # 80deg
# arrival_times = [4470, 4463, 4472, 4465, 4470, 4473]  # 90deg
#arrival_times = [9055, 9049, 9054, 9091, 9054, 9057]
#arrival_times = [5466, 5466, 5463, 5457, 5454, 5463]  # 180deg
#arrival_times = [9426, 9432, 9437, 9434, 9428, 9432]  # 270deg
# arrival_times = [5469, 5754, 5480, 5478, 5471, 5478]  # 270deg - 45cm -> FAIL!
#arrival_times = [9443, 9445, 9432, 9422, 9431, 9441]  # 270deg - 30cm

#arrival_times = [8089, 8089, 8085, 8084, 8084, 8083] # 180deg
#arrival_times = [7036, 7029, 7029, 7029, 7038, 7040]
arrival_times = [7344, 7342, 7341, 7337, 7339, 7341] # 180deg
arrival_times = [57284, 57284, 57287, 57289, 57291, 57287]
arrival_times = [61412, 61417, 61418, 61418, 61414, 61413]
arrival_times = [132801, 132798, 132799, 132798, 132802, 132915]

#arrival_times = [x / fs for x in arrival_times]

alpha = 60
c = 343  # This is a gues, since its based on the temperature
length = 0.0465#0.0476

# Dimensions of array:
#distances = [47, 47, 47, 47, 47, 47]

def angle_mean(angles):
    """Computes the average angle given an array of angles in degrees."""
    complex_sum = 0
    for angle in angles:
        angleRadians = np.radians(angle)
        comp = np.exp(1j*angleRadians)

        complex_sum += comp
    return np.angle(complex_sum, deg=True) % 360


# Step 1 calculate actual time differences
tdoa = np.zeros(num_channels)

for m in range(0, num_channels):
    tdoa[m] = (arrival_times[m] - arrival_times[(m - 1) % num_channels]) / fs

# for m in range(0, num_channels):
#     tdoa[num_channels + m] = (arrival_times[m] - arrival_times[(m - 3) % num_channels]) / fs

theta = np.zeros(num_channels)

for i in range(0, num_channels):
    angle_pair = np.degrees(np.arcsin(np.clip(tdoa[i] * c / length, -1, 1)))  # compute the DOA based on one microphone pair
    index = i - 1
    if index == -1:
        index = num_channels - 1
    if tdoa[index] >= 0:
        angle_pair = 180 - angle_pair
    theta[i] = (angle_pair - (i - 1) * alpha) % 360

bo = angle_mean(theta)

doa = round(360 - angle_mean(theta), 3)

print("DOA: " + str(doa))







# def centered(array, new_length):
#     current_length = len(array)
#
#     start = (current_length - new_length) // 2
#     end = start + new_length
#
#     return array[start:end]
#
#
# def new_approach(in1, in2):
#     N = 11 # len(in1) + len(in2) - 1
#
#     # 1. Transform both to the frequency domain
#     fft_in1 = np.fft.fft(in1, n=N)
#     fft_in2 = np.fft.fft(in2, n=N)
#
#     # 2. Perform element wise multiplication:
#     fft_result = fft_in1 * fft_in2
#
#     # 3. Perform inverse FFT:
#     res = np.fft.ifft(fft_result, n=N)
#
#     # 4. Take centered output:
#     res_centered = centered(res, len(in1))
#
#     return res_centered
#
#
#
# # Test oaconcolbe shizzel:
# in1 = [0.0, 0.06256294442579424, 0.08877834406567583, 0.06320383312479018, 0.0005798516800439467, -0.06259346293527024]
# in2 = [0.019076552693702675, -0.08039789221575043, -0.05389123679777112, 0.05739130386915245, 0.078903003952558, -0.02346880621817295]
#
# res_original = oaconvolve(in1, in2, mode="same")
# res_own = new_approach(in1, in2)
#
# t = 10











