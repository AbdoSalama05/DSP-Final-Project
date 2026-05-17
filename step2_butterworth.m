%% Load ECG data
load('ecg_data.mat');

%% Design Butterworth High-Pass Filter
[butterworthHpNumerator, butterworthHpDenominator] = butter(4, 0.5 / (samplingFrequency / 2), 'high');

%% Design Notch Filter at 50 Hz
normalizedNotchFreq = 50 / (samplingFrequency / 2);
notchBandwidth = normalizedNotchFreq / 35;
[notchNumerator, notchDenominator] = iirnotch(normalizedNotchFreq, notchBandwidth);

%% Design Butterworth Low-Pass Filter
[butterworthLpNumerator, butterworthLpDenominator] = butter(5, 100 / (samplingFrequency / 2), 'low');

%% Plot Frequency Responses
figure;
freqz(butterworthHpNumerator, butterworthHpDenominator, 1024, samplingFrequency);
title('Butterworth High-Pass (Frequency Response)');
figure;
freqz(notchNumerator, notchDenominator, 1024, samplingFrequency);
title('Notch Filter 50 Hz (Frequency Response)');
figure;
freqz(butterworthLpNumerator, butterworthLpDenominator, 1024, samplingFrequency);
title('Butterworth Low-Pass (Frequency Response)');

%% Plot Impulse Responses
impulseResponseLength = 500;
impulseTimeAxis = (0:impulseResponseLength - 1) / samplingFrequency;
highPassImpulseResponse = impz(butterworthHpNumerator, butterworthHpDenominator, impulseResponseLength);
notchImpulseResponse = impz(notchNumerator, notchDenominator, impulseResponseLength);
lowPassImpulseResponse = impz(butterworthLpNumerator, butterworthLpDenominator, impulseResponseLength);
figure;
subplot(3, 1, 1);
plot(impulseTimeAxis, highPassImpulseResponse);
title('Impulse Response (Butterworth High-Pass)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
subplot(3, 1, 2);
plot(impulseTimeAxis, notchImpulseResponse);
title('Impulse Response (Notch 50 Hz)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
subplot(3, 1, 3);
plot(impulseTimeAxis, lowPassImpulseResponse);
title('Impulse Response (Butterworth Low-Pass)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

%% Plot Step Responses
highPassStepResponse = stepz(butterworthHpNumerator, butterworthHpDenominator, impulseResponseLength);
notchStepResponse = stepz(notchNumerator, notchDenominator, impulseResponseLength);
lowPassStepResponse = stepz(butterworthLpNumerator, butterworthLpDenominator, impulseResponseLength);
figure;
subplot(3, 1, 1);
plot(impulseTimeAxis, highPassStepResponse);
title('Step Response (Butterworth High-Pass)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
subplot(3, 1, 2);
plot(impulseTimeAxis, notchStepResponse);
title('Step Response (Notch 50 Hz)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
subplot(3, 1, 3);
plot(impulseTimeAxis, lowPassStepResponse);
title('Step Response (Butterworth Low-Pass)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

%% Pole-Zero Plots
figure;
subplot(1, 3, 1);
zplane(butterworthHpNumerator, butterworthHpDenominator);
title('Pole-Zero (Butterworth HP)');
subplot(1, 3, 2);
zplane(notchNumerator, notchDenominator);
title('Pole-Zero (Notch 50 Hz)');
subplot(1, 3, 3);
zplane(butterworthLpNumerator, butterworthLpDenominator);
title('Pole-Zero (Butterworth LP)');

%% Apply Filters
ecgAfterHighPass = filtfilt(butterworthHpNumerator, butterworthHpDenominator, noisyEcgSignal);
ecgAfterNotch = filtfilt(notchNumerator, notchDenominator, ecgAfterHighPass);
ecgCleanedButterworth = filtfilt(butterworthLpNumerator, butterworthLpDenominator, ecgAfterNotch);

%% Plot Filtering Steps
figure;
subplot(4, 1, 1);
plot(timeAxis, noisyEcgSignal);
title('Noisy ECG');
ylabel('mV');
grid on;
subplot(4, 1, 2);
plot(timeAxis, ecgAfterHighPass);
title('After Butterworth High Pass');
ylabel('mV');
grid on;
subplot(4, 1, 3);
plot(timeAxis, ecgAfterNotch);
title('After Notch Filter');
ylabel('mV');
grid on;
subplot(4, 1, 4);
plot(timeAxis, ecgCleanedButterworth);
title('After Butterworth Low Pass');
ylabel('mV');
xlabel('Time (s)');
grid on;

%% Print Coefficients
disp('Butterworth HP numerator b:');
disp(butterworthHpNumerator);
disp('Butterworth HP denominator a:');
disp(butterworthHpDenominator);
disp('Notch numerator b:');
disp(notchNumerator);
disp('Notch denominator a:');
disp(notchDenominator);
disp('Butterworth LP numerator b:');
disp(butterworthLpNumerator);
disp('Butterworth LP denominator a:');
disp(butterworthLpDenominator);

%% Verify Specs Compliance
[highPassFreqResponse, highPassFreqAxis] = freqz(butterworthHpNumerator, butterworthHpDenominator, 8192, samplingFrequency);
[lowPassFreqResponse, lowPassFreqAxis] = freqz(butterworthLpNumerator, butterworthLpDenominator, 8192, samplingFrequency);

hpPassbandGain = 20 * log10(abs(highPassFreqResponse(find(highPassFreqAxis >= 10, 1))));
hpStopbandAttenuation = 20 * log10(abs(highPassFreqResponse(find(highPassFreqAxis >= 0.1, 1))));
lpPassbandGain = 20 * log10(abs(lowPassFreqResponse(find(lowPassFreqAxis >= 90, 1))));
lpStopbandAttenuation = 20 * log10(abs(lowPassFreqResponse(find(lowPassFreqAxis >= 150, 1))));

% Exact notch attenuation at 50.000 Hz using direct frequency evaluation
notchExactFreqVector = exp(-1j * 2 * pi * 50 / samplingFrequency * (0:2));
notchExactAttenuation = 20 * log10(abs(notchNumerator * notchExactFreqVector') / abs(notchDenominator * notchExactFreqVector'));

fprintf('\nButterworth Specs Verification: \n');
fprintf('HP passband gain at 10 Hz  : %.3f dB\n', hpPassbandGain);
fprintf('HP stopband atten at 0.1 Hz: %.3f dB\n', hpStopbandAttenuation);
fprintf('LP passband gain at 90 Hz  : %.3f dB\n', lpPassbandGain);
fprintf('LP stopband atten at 150 Hz: %.3f dB\n', lpStopbandAttenuation);
fprintf('Notch exact atten at 50 Hz : %.3f dB\n', notchExactAttenuation);

%% Save
save('butter_result.mat', 'ecgCleanedButterworth', 'butterworthHpNumerator', 'butterworthHpDenominator', ...
    'notchNumerator', 'notchDenominator', 'butterworthLpNumerator', 'butterworthLpDenominator');
