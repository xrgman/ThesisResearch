% Read data from WAV file:
[x, fs] = audioread('encoding.wav');

% Set the desired SNR (in dB)
desired_snr = -6

% Add Gaussian noise to the audio signal
noisy_signal = awgn(x, desired_snr, 'measured');

% Plot the original and noisy signals
t = (0:length(x)-1) / fs; % Time vector
figure;

% Plot the original signal
subplot(2, 1, 1);
plot(t, x);
title('Original Audio Signal');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;

% Plot the noisy signal
subplot(2, 1, 2);
plot(t, noisy_signal);
title(['Noisy Audio Signal (SNR = ' num2str(desired_snr) ' dB)']);
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;


% % Write noisy audio data to new file:
audiowrite('encoding_snr-6.wav', noisy_signal, fs);