%% Load ECG Signal
addpath('C:\Users\PolaNasser\Downloads\wfdb-app-toolbox-0-10-0\mcode')
samplingFrequency = 360;
totalSamples = 3600;
[rawSignal, samplingFrequency, timeVector] = rdsamp('100', [], totalSamples);
cleanEcgSignal = rawSignal(:, 1);
timeAxis = timeVector(:);
cleanEcgSignal = cleanEcgSignal(:);

%% Plot Raw ECG
figure;
plot(timeAxis, cleanEcgSignal);
xlabel('Time (s)');
ylabel('mV');
title('Raw ECG Signal');
grid on;

%% Add Synthetic Noise
baselineWanderNoise = 0.3 * sin(2 * pi * 0.3 * timeAxis);
powerlineInterferenceNoise = 0.2 * sin(2 * pi * 50 * timeAxis);
emgRaw = 0.05 * randn(size(timeAxis));
[bEmg, aEmg] = butter(4, [20 150]/(samplingFrequency/2), 'bandpass');

emgMusclNoise = filtfilt(bEmg, aEmg, emgRaw);noisyEcgSignal = cleanEcgSignal + baselineWanderNoise + powerlineInterferenceNoise + emgMusclNoise;

%% Plot Clean vs Noisy
figure;
subplot(2, 1, 1);
plot(timeAxis, cleanEcgSignal);
title('Clean ECG');
xlabel('Time (s)');
ylabel('mV');
grid on;
subplot(2, 1, 2);
plot(timeAxis, noisyEcgSignal);
title('Noisy ECG');
xlabel('Time (s)');
ylabel('mV');
grid on;

%% PSD of Noisy Signal
[powerSpectralDensity, frequencyAxis] = pwelch(noisyEcgSignal, [], [], [], samplingFrequency);
figure;
semilogy(frequencyAxis, powerSpectralDensity);
xlabel('Frequency (Hz)');
ylabel('PSD');
title('PSD of Noisy ECG');
grid on;
xline(0.5, 'r--', 'Baseline cutoff');
xline(50, 'g--', 'Powerline 50 Hz');
xline(100, 'b--', 'EMG cutoff');

%% Save workspace so other files can use it
save('ecg_data.mat', 'cleanEcgSignal', 'noisyEcgSignal', 'timeAxis', 'samplingFrequency');
