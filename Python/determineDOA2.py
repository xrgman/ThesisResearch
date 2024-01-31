import numpy as np

fs = 22050
n = 6
alpha = 60
c = 343  # This is a gues, since its based on the temperature
length = 0.0475#0.0465#0.0476


#arrival_times = [7344, 7342, 7338, 7337, 7339, 7341]  # 180deg
#arrival_times = [114604, 114602, 114598, 114597, 114599, 114601]  # 180deg
#arrival_times = [57284, 57284, 57287, 57289, 57291, 57287]  # 0 deg
#arrival_times = [61412, 61417, 61418, 61418, 61414, 61413]
#arrival_times = [132801, 132798, 132799, 132798, 132802, 132915]

# 90 deg:
#arrival_times = [89792, 89788, 89789, 89789, 89793, 89909]
#arrival_times = [132801, 132798, 132799, 132798, 132802, 132915]
#arrival_times = [173761, 173758, 173759, 173758, 173761, 173762]
arrival_times = [57982, 57979, 57978, 57979, 57983, 57982]
#arrival_times = [115326, 115323, 115322, 115323, 115326, 115326]

# 0. Determine actual arrival times:
#arrival_times = [x / fs for x in arrival_times]

def doa_estimator(tdoa, l, c, M):
    """Estimates the DOA given a set of TDOA values, length between the
    microphones, speed of sound, and number of microphones."""
    alpha = 360 / M      # The angle between microphones
    theta = np.zeros(M)
    for i in range(M):
        test = np.arcsin(np.clip(tdoa[i] * c / l, -1, 1)) # Radians result

        angle_pair = np.degrees(np.arcsin(np.clip(tdoa[i] * c / l, -1, 1)))  # compute the DOA based on one microphone pair

        index = i - 1
        if index == -1:
            index = M-1
        if tdoa[M + index] >= 0:
            angle_pair = 180 - angle_pair
        theta[i] = (angle_pair - (i-1) * alpha) % 360

    return round(360 - angle_mean(theta), 3)


def angle_mean(angles):
    """Computes the average angle given an array of angles in degrees."""
    complex_sum = 0
    for angle in angles:
        complex_sum += np.exp(1j*np.radians(angle))
    return np.angle(complex_sum, deg=True) % 360


def determine_doa(arr_times):
    tdoa = np.zeros(n * 2)

    for m in range(n):
        tdoa[m] = (arr_times[m] - arr_times[(m - 1) % n]) / fs

    for m in range(n):
        tdoa[n + m] = (arr_times[m] - arr_times[(m - 3) % n]) / fs

    return doa_estimator(tdoa, length, c, n)

# doa = determine_doa(arrival_times)
#
# print("Doa: " + str(doa))