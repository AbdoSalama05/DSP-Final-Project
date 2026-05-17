%% Load ECG data
load('ecg_data.mat');

%% Design Chebyshev Type I High-Pass Filter
passbandRippleDb = 0.25;
[chebyshevHpNumerator, chebyshevHpDenominator] = cheby1(4, passbandRippleDb, 0.5 / (samplingFrequency / 2), 'high');

%% Design Notch Filter at 50 Hz
normalizedNotchFreq = 50 / (samplingFrequency / 2);
notchBandwidth = normalizedNotchFreq / 50;
[chebyshevNotchNumerator, chebyshevNotchDenominator] = iirnotch(normalizedNotchFreq, notchBandwidth);

%% Design Chebyshev Type I Low-Pass Filter
[chebyshevLpNumerator, chebyshevLpDenominator] = cheby1(4, passbandRippleDb, 100 / (samplingFrequency / 2), 'low');

%% Plot Frequency Responses
figure;
freqz(chebyshevHpNumerator, chebyshevHpDenominator, 1024, samplingFrequency);
title('Chebyshev High-Pass (Frequency Response)');
figure;
freqz(chebyshevNotchNumerator, chebyshevNotchDenominator, 1024, samplingFrequency);
title('Notch Filter 50 Hz (Frequency Response)');
figure;
freqz(chebyshevLpNumerator, chebyshevLpDenominator, 1024, samplingFrequency);
title('Chebyshev Low-Pass (Frequency Response)');

%% Plot Impulse Responses
impulseResponseLength = 500;
impulseTimeAxis = (0:impulseResponseLength - 1) / samplingFrequency;
highPassImpulseResponse = impz(chebyshevHpNumerator, chebyshevHpDenominator, impulseResponseLength);
notchImpulseResponse = impz(chebyshevNotchNumerator, chebyshevNotchDenominator, impulseResponseLength);
lowPassImpulseResponse = impz(chebyshevLpNumerator, chebyshevLpDenominator, impulseResponseLength);
figure;
subplot(3, 1, 1);
plot(impulseTimeAxis, highPassImpulseResponse);
title('Impulse Response (Chebyshev High-Pass)');
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
title('Impulse Response (Chebyshev Low-Pass)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

%% Plot Step Responses
highPassStepResponse = stepz(chebyshevHpNumerator, chebyshevHpDenominator, impulseResponseLength);
notchStepResponse = stepz(chebyshevNotchNumerator, chebyshevNotchDenominator, impulseResponseLength);
lowPassStepResponse = stepz(chebyshevLpNumerator, chebyshevLpDenominator, impulseResponseLength);
figure;
subplot(3, 1, 1);
plot(impulseTimeAxis, highPassStepResponse);
title('Step Response (Chebyshev High-Pass)');
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
title('Step Response (Chebyshev Low-Pass)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

%% Pole-Zero Plots
figure;
subplot(1, 3, 1);
zplane(chebyshevHpNumerator, chebyshevHpDenominator);
title('Pole-Zero (Chebyshev HP)');
subplot(1, 3, 2);
zplane(chebyshevNotchNumerator, chebyshevNotchDenominator);
title('Pole-Zero (Notch 50 Hz)');
subplot(1, 3, 3);
zplane(chebyshevLpNumerator, chebyshevLpDenominator);
title('Pole-Zero (Chebyshev LP)');

%% Apply Filters
ecgAfterHighPass = filtfilt(chebyshevHpNumerator, chebyshevHpDenominator, noisyEcgSignal);
ecgAfterNotch = filtfilt(chebyshevNotchNumerator, chebyshevNotchDenominator, ecgAfterHighPass);
ecgCleanedChebyshev = filtfilt(chebyshevLpNumerator, chebyshevLpDenominator, ecgAfterNotch);

%% Plot Filtering Steps
figure;
subplot(4, 1, 1);
plot(timeAxis, noisyEcgSignal);
title('Noisy ECG');
ylabel('mV');
grid on;
subplot(4, 1, 2);
plot(timeAxis, ecgAfterHighPass);
title('After Chebyshev High Pass');
ylabel('mV');
grid on;
subplot(4, 1, 3);
plot(timeAxis, ecgAfterNotch);
title('After Notch Filter');
ylabel('mV');
grid on;
subplot(4, 1, 4);
plot(timeAxis, ecgCleanedChebyshev);
title('After Chebyshev Low-Pass');
ylabel('mV');
xlabel('Time (s)');
grid on;

%% Print Coefficients
disp('Chebyshev HP numerator b:');
disp(chebyshevHpNumerator);
disp('Chebyshev HP denominator a:');
disp(chebyshevHpDenominator);
disp('Notch numerator b:');
disp(chebyshevNotchNumerator);
disp('Notch denominator a:');
disp(chebyshevNotchDenominator);
disp('Chebyshev LP numerator b:');
disp(chebyshevLpNumerator);
disp('Chebyshev LP denominator a:');
disp(chebyshevLpDenominator);

%% Verify Specs Compliance
[highPassFreqResponse, highPassFreqAxis] = freqz(chebyshevHpNumerator, chebyshevHpDenominator, 8192, samplingFrequency);
[lowPassFreqResponse, lowPassFreqAxis] = freqz(chebyshevLpNumerator, chebyshevLpDenominator, 8192, samplingFrequency);

hpPassbandGain = 20 * log10(abs(highPassFreqResponse(find(highPassFreqAxis >= 10, 1))));
hpStopbandAttenuation = 20 * log10(abs(highPassFreqResponse(find(highPassFreqAxis >= 0.1, 1))));
lpPassbandGain = 20 * log10(abs(lowPassFreqResponse(find(lowPassFreqAxis >= 90, 1))));
lpStopbandAttenuation = 20 * log10(abs(lowPassFreqResponse(find(lowPassFreqAxis >= 150, 1))));

notchExactFreqVector = exp(-1j * 2 * pi * 50 / samplingFrequency * (0:2));
notchExactAttenuation = 20 * log10(abs(chebyshevNotchNumerator * notchExactFreqVector') / abs(chebyshevNotchDenominator * notchExactFreqVector'));

fprintf('\nChebyshev Specs Verification: \n');
fprintf('HP passband gain at 10 Hz  : %.3f dB\n', hpPassbandGain);
fprintf('HP stopband atten at 0.1 Hz: %.3f dB\n', hpStopbandAttenuation);
fprintf('LP passband gain at 90 Hz  : %.3f dB\n', lpPassbandGain);
fprintf('LP stopband atten at 150 Hz: %.3f dB\n', lpStopbandAttenuation);
fprintf('Notch exact atten at 50 Hz : %.3f dB\n', notchExactAttenuation);

%% Save
save('cheby_result.mat', 'ecgCleanedChebyshev', 'chebyshevHpNumerator', 'chebyshevHpDenominator', ...
    'chebyshevNotchNumerator', 'chebyshevNotchDenominator', 'chebyshevLpNumerator', 'chebyshevLpDenominator');
