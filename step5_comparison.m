%% Load all results
load('ecg_data.mat');
load('butter_result.mat');
load('cheby_result.mat');
load('fir_result.mat');

%% Compare Clean vs All Filtered (Time Domain)
figure;
subplot(4, 1, 1);
plot(timeAxis, noisyEcgSignal);
title('Noisy ECG');
ylabel('mV');
grid on;
subplot(4, 1, 2);
plot(timeAxis, ecgCleanedButterworth);
title('Butterworth Filtered ECG');
ylabel('mV');
grid on;
subplot(4, 1, 3);
plot(timeAxis, ecgCleanedChebyshev);
title('Chebyshev Filtered ECG');
ylabel('mV');
grid on;
subplot(4, 1, 4);
plot(timeAxis, ecgCleanedFir);
title('FIR Filtered ECG');
ylabel('mV');
xlabel('Time (s)');
grid on;

%% PSD Comparison
[psdNoisySignal, psdFreqAxis] = pwelch(noisyEcgSignal, hamming(512), 256, 1024, samplingFrequency);
[psdButterworthSignal, ~] = pwelch(ecgCleanedButterworth, hamming(512), 256, 1024, samplingFrequency);
[psdChebyshevSignal, ~] = pwelch(ecgCleanedChebyshev, hamming(512), 256, 1024, samplingFrequency);
[psdFirSignal, ~] = pwelch(ecgCleanedFir, hamming(512), 256, 1024, samplingFrequency);
figure;
plot(psdFreqAxis, 10 * log10(psdNoisySignal), 'k');
hold on;
plot(psdFreqAxis, 10 * log10(psdButterworthSignal), 'b');
plot(psdFreqAxis, 10 * log10(psdChebyshevSignal), 'r');
plot(psdFreqAxis, 10 * log10(psdFirSignal), 'g');
legend('Noisy', 'Butterworth', 'Chebyshev', 'FIR');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('PSD Comparison: All Filters');
xlim([0 180]);
xline(0.5, 'k--', '0.5 Hz');
xline(50, 'k--', '50 Hz');
xline(100, 'k--', '100 Hz');
grid on;

%% Spectrograms
figure;
subplot(2, 2, 1);
spectrogram(noisyEcgSignal, hamming(256), 128, 512, samplingFrequency, 'yaxis');
title('Spectrogram (Noisy ECG)');
subplot(2, 2, 2);
spectrogram(ecgCleanedButterworth, hamming(256), 128, 512, samplingFrequency, 'yaxis');
title('Spectrogram (Butterworth)');
subplot(2, 2, 3);
spectrogram(ecgCleanedChebyshev, hamming(256), 128, 512, samplingFrequency, 'yaxis');
title('Spectrogram (Chebyshev)');
subplot(2, 2, 4);
spectrogram(ecgCleanedFir, hamming(256), 128, 512, samplingFrequency, 'yaxis');
title('Spectrogram (FIR)');

%% SNR Calculation (Correct Method)
% Input SNR: ratio of clean signal power to added noise power
addedNoise = noisyEcgSignal - cleanEcgSignal;
inputSnrDb = 10 * log10(sum(cleanEcgSignal.^2) / sum(addedNoise.^2));

% Apply the same filter pipeline to the clean signal as a reference

% Butterworth reference
cleanAfterButterworthHp = filtfilt(butterworthHpNumerator, butterworthHpDenominator, cleanEcgSignal);
cleanAfterButterworthNotch = filtfilt(notchNumerator, notchDenominator, cleanAfterButterworthHp);
cleanReferenceButter = filtfilt(butterworthLpNumerator, butterworthLpDenominator, cleanAfterButterworthNotch);

% Chebyshev reference
cleanAfterChebyshevHp = filtfilt(chebyshevHpNumerator, chebyshevHpDenominator, cleanEcgSignal);
cleanAfterChebyshevNotch = filtfilt(chebyshevNotchNumerator, chebyshevNotchDenominator, cleanAfterChebyshevHp);
cleanReferenceCheby = filtfilt(chebyshevLpNumerator, chebyshevLpDenominator, cleanAfterChebyshevNotch);

% FIR reference
cleanAfterFirHp = filtfilt(firHighPassNumerator, firHighPassDenominator, cleanEcgSignal);
cleanAfterFirNotch = filtfilt(firNotchNumerator, firNotchDenominator, cleanAfterFirHp);
cleanReferenceFir = filtfilt(firLowPassNumerator, firLowPassDenominator, cleanAfterFirNotch);

% Output SNR: ratio of filtered clean signal power to residual noise power
residualNoiseButter = cleanReferenceButter - ecgCleanedButterworth;
residualNoiseCheby = cleanReferenceCheby - ecgCleanedChebyshev;
residualNoiseFir = cleanReferenceFir - ecgCleanedFir;

outputSnrButterworth = 10 * log10(sum(cleanReferenceButter.^2) / sum(residualNoiseButter.^2));
outputSnrChebyshev = 10 * log10(sum(cleanReferenceCheby.^2) / sum(residualNoiseCheby.^2));
outputSnrFir = 10 * log10(sum(cleanReferenceFir.^2) / sum(residualNoiseFir.^2));

snrImprovementButterworth = outputSnrButterworth - inputSnrDb;
snrImprovementChebyshev = outputSnrChebyshev - inputSnrDb;
snrImprovementFir = outputSnrFir - inputSnrDb;

fprintf('\nSNR Before Filtering (Input) : %.3f dB\n', inputSnrDb);
fprintf('SNR After Butterworth        : %.3f dB\n', outputSnrButterworth);
fprintf('SNR After Chebyshev Type I   : %.3f dB\n', outputSnrChebyshev);
fprintf('SNR After FIR                : %.3f dB\n', outputSnrFir);

fprintf('\nSNR Improvement Butterworth  : %+.3f dB\n', snrImprovementButterworth);
fprintf('SNR Improvement Chebyshev    : %+.3f dB\n', snrImprovementChebyshev);
fprintf('SNR Improvement FIR          : %+.3f dB\n', snrImprovementFir);

figure;
snrBarData = [inputSnrDb, outputSnrButterworth, outputSnrChebyshev, outputSnrFir];
bar(snrBarData);
set(gca, 'XTickLabel', {'Input (Noisy)', 'Butterworth', 'Chebyshev', 'FIR'});
ylabel('SNR (dB)');
title('SNR Comparison: Before and After Filtering');
grid on;
for barIndex = 1:4
    text(barIndex, snrBarData(barIndex) + 0.3, sprintf('%.2f dB', snrBarData(barIndex)), ...
        'HorizontalAlignment', 'center', 'FontSize', 9);
end

figure;
snrImprovementData = [snrImprovementButterworth, snrImprovementChebyshev, snrImprovementFir];
bar(snrImprovementData);
set(gca, 'XTickLabel', {'Butterworth', 'Chebyshev', 'FIR'});
ylabel('\DeltaSNR (dB)');
title('SNR Improvement Over Noisy Input');
grid on;
for barIndex = 1:3
    text(barIndex, snrImprovementData(barIndex) + 0.1, sprintf('%+.2f dB', snrImprovementData(barIndex)), ...
        'HorizontalAlignment', 'center', 'FontSize', 9);
end

%% Phase Response Comparison
freqResponsePoints = 8192;

[butterworthHpFreqResponse, highPassFreqAxis] = freqz(butterworthHpNumerator, butterworthHpDenominator, freqResponsePoints, samplingFrequency);
[chebyshevHpFreqResponse, ~] = freqz(chebyshevHpNumerator, chebyshevHpDenominator, freqResponsePoints, samplingFrequency);
[firHpFreqResponse, ~] = freqz(firHighPassNumerator, 1, freqResponsePoints, samplingFrequency);

[butterworthLpFreqResponse, lowPassFreqAxis] = freqz(butterworthLpNumerator, butterworthLpDenominator, freqResponsePoints, samplingFrequency);
[chebyshevLpFreqResponse, ~] = freqz(chebyshevLpNumerator, chebyshevLpDenominator, freqResponsePoints, samplingFrequency);
[firLpFreqResponse, ~] = freqz(firLowPassNumerator, 1, freqResponsePoints, samplingFrequency);

[notchFreqResponse, notchFreqAxis] = freqz(notchNumerator, notchDenominator, freqResponsePoints, samplingFrequency);

figure;
subplot(2, 1, 1);
plot(highPassFreqAxis, 20 * log10(abs(butterworthHpFreqResponse)), 'b', 'LineWidth', 1.4); hold on;
plot(highPassFreqAxis, 20 * log10(abs(chebyshevHpFreqResponse)), 'r', 'LineWidth', 1.4);
plot(highPassFreqAxis, 20 * log10(abs(firHpFreqResponse)), 'g', 'LineWidth', 1.4);
xline(0.5, 'k--', '0.5 Hz'); xline(0.1, 'm--', '0.1 Hz');
ylim([-80 5]); xlim([0 5]);
legend('Butterworth HP', 'Chebyshev HP', 'FIR HP', 'Location', 'southeast');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('High-Pass: Magnitude Response Comparison');
grid on;
subplot(2, 1, 2);
plot(highPassFreqAxis, unwrap(angle(butterworthHpFreqResponse)) * 180 / pi, 'b', 'LineWidth', 1.4); hold on;
plot(highPassFreqAxis, unwrap(angle(chebyshevHpFreqResponse)) * 180 / pi, 'r', 'LineWidth', 1.4);
plot(highPassFreqAxis, unwrap(angle(firHpFreqResponse)) * 180 / pi, 'g', 'LineWidth', 1.4);
xline(0.5, 'k--', '0.5 Hz'); xlim([0 10]);
legend('Butterworth HP', 'Chebyshev HP', 'FIR HP', 'Location', 'northeast');
xlabel('Frequency (Hz)'); ylabel('Phase (degrees)');
title('High-Pass: Phase Response Comparison');
grid on;

figure;
subplot(2, 1, 1);
plot(lowPassFreqAxis, 20 * log10(abs(butterworthLpFreqResponse)), 'b', 'LineWidth', 1.4); hold on;
plot(lowPassFreqAxis, 20 * log10(abs(chebyshevLpFreqResponse)), 'r', 'LineWidth', 1.4);
plot(lowPassFreqAxis, 20 * log10(abs(firLpFreqResponse)), 'g', 'LineWidth', 1.4);
xline(100, 'k--', '100 Hz'); xline(150, 'm--', '150 Hz');
ylim([-100 5]); xlim([0 180]);
legend('Butterworth LP', 'Chebyshev LP', 'FIR LP', 'Location', 'southwest');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Low-Pass: Magnitude Response Comparison');
grid on;
subplot(2, 1, 2);
plot(lowPassFreqAxis, unwrap(angle(butterworthLpFreqResponse)) * 180 / pi, 'b', 'LineWidth', 1.4); hold on;
plot(lowPassFreqAxis, unwrap(angle(chebyshevLpFreqResponse)) * 180 / pi, 'r', 'LineWidth', 1.4);
plot(lowPassFreqAxis, unwrap(angle(firLpFreqResponse)) * 180 / pi, 'g', 'LineWidth', 1.4);
xline(100, 'k--', '100 Hz'); xlim([0 180]);
legend('Butterworth LP', 'Chebyshev LP', 'FIR LP', 'Location', 'southwest');
xlabel('Frequency (Hz)'); ylabel('Phase (degrees)');
title('Low-Pass: Phase Response Comparison');
grid on;

figure;
subplot(2, 1, 1);
plot(notchFreqAxis, 20 * log10(abs(notchFreqResponse)), 'k', 'LineWidth', 1.4);
xline(50, 'r--', '50 Hz'); xlim([40 60]); ylim([-80 5]);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Notch Filter: Magnitude Response - Zoomed');
grid on;
subplot(2, 1, 2);
plot(notchFreqAxis, unwrap(angle(notchFreqResponse)) * 180 / pi, 'k', 'LineWidth', 1.4);
xline(50, 'r--', '50 Hz'); xlim([40 60]);
xlabel('Frequency (Hz)'); ylabel('Phase (degrees)');
title('Notch Filter: Phase Response - Zoomed');
grid on;

figure;
subplot(2, 1, 1);
[groupDelayButterworthHp, groupDelayFreqAxis] = grpdelay(butterworthHpNumerator, butterworthHpDenominator, freqResponsePoints, samplingFrequency);
[groupDelayChebyshevHp, ~] = grpdelay(chebyshevHpNumerator, chebyshevHpDenominator, freqResponsePoints, samplingFrequency);
[groupDelayFirHp, ~] = grpdelay(firHighPassNumerator, 1, freqResponsePoints, samplingFrequency);
plot(groupDelayFreqAxis, groupDelayButterworthHp, 'b', 'LineWidth', 1.4); hold on;
plot(groupDelayFreqAxis, groupDelayChebyshevHp, 'r', 'LineWidth', 1.4);
plot(groupDelayFreqAxis, groupDelayFirHp, 'g', 'LineWidth', 1.4);
xline(0.5, 'k--', '0.5 Hz'); xlim([0 10]);
legend('Butterworth HP', 'Chebyshev HP', 'FIR HP');
xlabel('Frequency (Hz)'); ylabel('Group Delay (samples)');
title('High-Pass: Group Delay Comparison');
grid on;
subplot(2, 1, 2);
[groupDelayButterworthLp, groupDelayLpFreqAxis] = grpdelay(butterworthLpNumerator, butterworthLpDenominator, freqResponsePoints, samplingFrequency);
[groupDelayChebyshevLp, ~] = grpdelay(chebyshevLpNumerator, chebyshevLpDenominator, freqResponsePoints, samplingFrequency);
[groupDelayFirLp, ~] = grpdelay(firLowPassNumerator, 1, freqResponsePoints, samplingFrequency);
plot(groupDelayLpFreqAxis, groupDelayButterworthLp, 'b', 'LineWidth', 1.4); hold on;
plot(groupDelayLpFreqAxis, groupDelayChebyshevLp, 'r', 'LineWidth', 1.4);
plot(groupDelayLpFreqAxis, groupDelayFirLp, 'g', 'LineWidth', 1.4);
xline(100, 'k--', '100 Hz'); xlim([0 180]);
legend('Butterworth LP', 'Chebyshev LP', 'FIR LP');
xlabel('Frequency (Hz)'); ylabel('Group Delay (samples)');
title('Low-Pass: Group Delay Comparison');
grid on;

%% Removed Noise Plots
removedNoiseButter = noisyEcgSignal - ecgCleanedButterworth;
removedNoiseCheby = noisyEcgSignal - ecgCleanedChebyshev;
removedNoiseFir = noisyEcgSignal - ecgCleanedFir;
figure;
subplot(3, 1, 1);
plot(timeAxis, removedNoiseButter);
title('Removed Noise - Butterworth');
xlabel('Time (s)'); ylabel('mV');
grid on;
subplot(3, 1, 2);
plot(timeAxis, removedNoiseCheby);
title('Removed Noise - Chebyshev');
xlabel('Time (s)'); ylabel('mV');
grid on;
subplot(3, 1, 3);
plot(timeAxis, removedNoiseFir);
title('Removed Noise - FIR');
xlabel('Time (s)'); ylabel('mV');
grid on;
