%% Load ECG data
load('ecg_data.mat');

%% Design FIR Filters
firOrder = 500;
firHighPassOrder = firOrder;
firStandardOrder = firOrder;

% High-pass FIR filter
firHighPassNumerator = fir1(firOrder, 0.5 / (samplingFrequency / 2), 'high');
firHighPassDenominator = 1;

% Notch FIR band-stop filter at 50 Hz
firNotchNumerator = fir1(firOrder, [49 51] / (samplingFrequency / 2), 'stop');
firNotchDenominator = 1;

% Low-pass FIR filter
firLowPassNumerator = fir1(firOrder, 100 / (samplingFrequency / 2), 'low');
firLowPassDenominator = 1;

%% Plot Frequency Responses
figure;
freqz(firHighPassNumerator, firHighPassDenominator, 1024, samplingFrequency);
title('FIR High-Pass - Frequency Response');
figure;
freqz(firNotchNumerator, firNotchDenominator, 1024, samplingFrequency);
title('IIR Notch 50 Hz - Frequency Response (used in FIR pipeline)');
figure;
freqz(firLowPassNumerator, firLowPassDenominator, 1024, samplingFrequency);
title('FIR Low-Pass - Frequency Response');

%% Plot Impulse Responses
impulseResponseLength = 500;
impulseTimeAxis = (0:impulseResponseLength - 1) / samplingFrequency;
firHighPassImpulseResponse = impz(firHighPassNumerator, firHighPassDenominator, impulseResponseLength);
notchImpulseResponse = impz(firNotchNumerator, firNotchDenominator, impulseResponseLength);
firLowPassImpulseResponse = impz(firLowPassNumerator, firLowPassDenominator, impulseResponseLength);
figure;
subplot(3, 1, 1);
plot(impulseTimeAxis, firHighPassImpulseResponse);
title('Impulse Response - FIR High-Pass');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
subplot(3, 1, 2);
plot(impulseTimeAxis, notchImpulseResponse);
title('Impulse Response - Notch 50 Hz');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
subplot(3, 1, 3);
plot(impulseTimeAxis, firLowPassImpulseResponse);
title('Impulse Response - FIR Low-Pass');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

%% Plot Step Responses
firHighPassStepResponse = stepz(firHighPassNumerator, firHighPassDenominator, impulseResponseLength);
notchStepResponse = stepz(firNotchNumerator, firNotchDenominator, impulseResponseLength);
firLowPassStepResponse = stepz(firLowPassNumerator, firLowPassDenominator, impulseResponseLength);
figure;
subplot(3, 1, 1);
plot(impulseTimeAxis, firHighPassStepResponse);
title('Step Response - FIR High-Pass');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
subplot(3, 1, 2);
plot(impulseTimeAxis, notchStepResponse);
title('Step Response - Notch 50 Hz');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
subplot(3, 1, 3);
plot(impulseTimeAxis, firLowPassStepResponse);
title('Step Response - FIR Low-Pass');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

%% Pole-Zero Plots
figure;
subplot(1, 3, 1);
zplane(firHighPassNumerator, firHighPassDenominator);
title('Pole-Zero - FIR HP');
subplot(1, 3, 2);
zplane(firNotchNumerator, firNotchDenominator);
title('Pole-Zero - Notch 50 Hz');
subplot(1, 3, 3);
zplane(firLowPassNumerator, firLowPassDenominator);
title('Pole-Zero - FIR LP');

%% Apply Filters
ecgAfterHighPass = filtfilt(firHighPassNumerator, firHighPassDenominator, noisyEcgSignal);
ecgAfterNotch = filtfilt(firNotchNumerator, firNotchDenominator, ecgAfterHighPass);
ecgCleanedFir = filtfilt(firLowPassNumerator, firLowPassDenominator, ecgAfterNotch);

%% Plot Filtering Steps
figure;
subplot(4, 1, 1);
plot(timeAxis, noisyEcgSignal);
title('Noisy ECG');
ylabel('mV');
grid on;
subplot(4, 1, 2);
plot(timeAxis, ecgAfterHighPass);
title('After FIR High-Pass');
ylabel('mV');
grid on;
subplot(4, 1, 3);
plot(timeAxis, ecgAfterNotch);
title('After FIR Notch');
ylabel('mV');
grid on;
subplot(4, 1, 4);
plot(timeAxis, ecgCleanedFir);
title('After FIR Low-Pass');
ylabel('mV');
xlabel('Time (s)');
grid on;

%% Print Coefficients
disp('FIR HP first 10 b coefficients:');
disp(firHighPassNumerator(1:10));
disp('Notch numerator b:');
disp(firNotchNumerator);
disp('Notch denominator a:');
disp(firNotchDenominator);
disp('FIR LP first 10 b coefficients:');
disp(firLowPassNumerator(1:10));

%% Verify Specs Compliance
firHpSampleIndices = (0:firHighPassOrder);
firStdSampleIndices = (0:firStandardOrder);
notchSampleIndices = (0:firOrder);

hpPassbandGain = 20 * log10(abs(firHighPassNumerator * exp(-1j * 2 * pi * 10 / samplingFrequency * firHpSampleIndices)'));
hpStopbandAttenuation_01 = 20 * log10(abs(firHighPassNumerator * exp(-1j * 2 * pi * 0.1 / samplingFrequency * firHpSampleIndices)'));
lpPassbandGain = 20 * log10(abs(firLowPassNumerator * exp(-1j * 2 * pi * 90 / samplingFrequency * firStdSampleIndices)'));
lpStopbandAttenuation = 20 * log10(abs(firLowPassNumerator * exp(-1j * 2 * pi * 150 / samplingFrequency * firStdSampleIndices)'));
notchExactFreqVector = exp(-1j * 2 * pi * 50 / samplingFrequency * notchSampleIndices);
notchExactAttenuation = 20 * log10(abs(firNotchNumerator * notchExactFreqVector'));

fprintf('\nFIR Specs Verification: \n');
fprintf('HP passband gain at 10 Hz  : %.3f dB\n', hpPassbandGain);
fprintf('HP stopband atten at 0.1 Hz: %.3f dB\n', hpStopbandAttenuation_01);
fprintf('LP passband gain at 90 Hz  : %.3f dB\n', lpPassbandGain);
fprintf('LP stopband atten at 150 Hz: %.3f dB\n', lpStopbandAttenuation);
fprintf('Notch exact atten at 50 Hz : %.3f dB\n', notchExactAttenuation);

%% Save
save('fir_result.mat', 'ecgCleanedFir', 'firHighPassNumerator', 'firHighPassDenominator', ...
    'firNotchNumerator', 'firNotchDenominator', 'firLowPassNumerator', 'firLowPassDenominator');
