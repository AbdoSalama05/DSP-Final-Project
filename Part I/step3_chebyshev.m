%% Step 3: Chebyshev Type I Filter Design and Application
load('ecg_data.mat');
%% Design Filters
passbandRippleDb=1;
[chebyshevHpNumerator,chebyshevHpDenominator]=cheby1(4,passbandRippleDb,0.5/(samplingFrequency/2),'high');

normalizedNotchFreq=50/(samplingFrequency/2);
notchBandwidth=normalizedNotchFreq/50;

[chebyshevNotchNumerator,chebyshevNotchDenominator]=iirnotch(normalizedNotchFreq,notchBandwidth);
[chebyshevLpNumerator,chebyshevLpDenominator]=cheby1(4,passbandRippleDb,100/(samplingFrequency/2),'low');
%% Frequency Responses
figure;freqz(chebyshevHpNumerator,chebyshevHpDenominator,1024,samplingFrequency);
title('Chebyshev High-Pass (Frequency Response)');
figure;freqz(chebyshevNotchNumerator,chebyshevNotchDenominator,1024,samplingFrequency);
title('Notch Filter 50 Hz (Frequency Response)');
figure;freqz(chebyshevLpNumerator,chebyshevLpDenominator,1024,samplingFrequency);
title('Chebyshev Low-Pass (Frequency Response)');
%% Impulse Responses
impulseResponseLength=500;
impulseTimeAxis=(0:impulseResponseLength-1)/samplingFrequency;
figure;
subplot(3,1,1);plot(impulseTimeAxis,impz(chebyshevHpNumerator,chebyshevHpDenominator,impulseResponseLength));
title('Impulse Response (Chebyshev High-Pass)');xlabel('Time (s)');ylabel('Amplitude');grid on;
subplot(3,1,2);plot(impulseTimeAxis,impz(chebyshevNotchNumerator,chebyshevNotchDenominator,impulseResponseLength));
title('Impulse Response (Notch 50 Hz)');xlabel('Time (s)');ylabel('Amplitude');grid on;
subplot(3,1,3);plot(impulseTimeAxis,impz(chebyshevLpNumerator,chebyshevLpDenominator,impulseResponseLength));
title('Impulse Response (Chebyshev Low-Pass)');xlabel('Time (s)');ylabel('Amplitude');grid on;
%% Step Responses
figure;
subplot(3,1,1);plot(impulseTimeAxis,stepz(chebyshevHpNumerator,chebyshevHpDenominator,impulseResponseLength));
title('Step Response (Chebyshev High-Pass)');xlabel('Time (s)');ylabel('Amplitude');grid on;
subplot(3,1,2);plot(impulseTimeAxis,stepz(chebyshevNotchNumerator,chebyshevNotchDenominator,impulseResponseLength));
title('Step Response (Notch 50 Hz)');xlabel('Time (s)');ylabel('Amplitude');grid on;
subplot(3,1,3);plot(impulseTimeAxis,stepz(chebyshevLpNumerator,chebyshevLpDenominator,impulseResponseLength));
title('Step Response (Chebyshev Low-Pass)');xlabel('Time (s)');ylabel('Amplitude');grid on;
%% Pole-Zero Plots
figure;
subplot(1,3,1);zplane(chebyshevHpNumerator,chebyshevHpDenominator);title('Pole-Zero (Chebyshev HP)');
subplot(1,3,2);zplane(chebyshevNotchNumerator,chebyshevNotchDenominator);title('Pole-Zero (Notch 50 Hz)');
subplot(1,3,3);zplane(chebyshevLpNumerator,chebyshevLpDenominator);title('Pole-Zero (Chebyshev LP)');
%% Apply Filters
ecgAfterHighPass=filtfilt(chebyshevHpNumerator,chebyshevHpDenominator,noisyEcgSignal);
ecgAfterNotch=filtfilt(chebyshevNotchNumerator,chebyshevNotchDenominator,ecgAfterHighPass);
ecgCleanedChebyshev=filtfilt(chebyshevLpNumerator,chebyshevLpDenominator,ecgAfterNotch);
%% Time-Domain Filtering Steps
figure;
subplot(4,1,1);plot(timeAxis,noisyEcgSignal);title('Noisy ECG');ylabel('mV');grid on;
subplot(4,1,2);plot(timeAxis,ecgAfterHighPass);title('After Chebyshev High-Pass');ylabel('mV');grid on;
subplot(4,1,3);plot(timeAxis,ecgAfterNotch);title('After Notch Filter');ylabel('mV');grid on;
subplot(4,1,4);plot(timeAxis,ecgCleanedChebyshev);title('After Chebyshev Low-Pass');ylabel('mV');xlabel('Time (s)');grid on;
%% PSD: Noisy vs Filtered
[psdNoisySignal,psdFreqAxis]=pwelch(noisyEcgSignal,hamming(512),256,1024,samplingFrequency);
[psdChebyshevSignal,~]=pwelch(ecgCleanedChebyshev,hamming(512),256,1024,samplingFrequency);
figure;
plot(psdFreqAxis,10*log10(psdNoisySignal),'r');hold on;
plot(psdFreqAxis,10*log10(psdChebyshevSignal),'b');
legend('Noisy ECG','Chebyshev Filtered');
xlabel('Frequency (Hz)');ylabel('Power/Frequency (dB/Hz)');
title('PSD: Noisy vs Chebyshev Filtered');xlim([0 180]);grid on;
xline(0.5,'--m','0.5 Hz');xline(50,'--g','50 Hz');xline(100,'--c','100 Hz');
%% Spectrogram
figure;
subplot(1,2,1);spectrogram(noisyEcgSignal,hamming(256),128,512,samplingFrequency,'yaxis');title('Spectrogram (Noisy ECG)');
subplot(1,2,2);spectrogram(ecgCleanedChebyshev,hamming(256),128,512,samplingFrequency,'yaxis');title('Spectrogram (Chebyshev Filtered)');
%% SNR
removedNoiseChebyshev=noisyEcgSignal-ecgCleanedChebyshev;
outputSnrChebyshev=10*log10(sum(ecgCleanedChebyshev.^2)/sum(removedNoiseChebyshev.^2));
fprintf('SNR indicator : %.2f dB\n',outputSnrChebyshev);
%% Removed Noise Plot
figure;
plot(timeAxis,removedNoiseChebyshev);
title('Removed Noise - Chebyshev');xlabel('Time (s)');ylabel('mV');grid on;
%% Print Coefficients
disp('Chebyshev HP numerator b:');disp(chebyshevHpNumerator);
disp('Chebyshev HP denominator a:');disp(chebyshevHpDenominator);
disp('Notch numerator b:');disp(chebyshevNotchNumerator);
disp('Notch denominator a:');disp(chebyshevNotchDenominator);
disp('Chebyshev LP numerator b:');disp(chebyshevLpNumerator);
disp('Chebyshev LP denominator a:');disp(chebyshevLpDenominator);
%% Verify Specs
[highPassFreqResponse,highPassFreqAxis]=freqz(chebyshevHpNumerator,chebyshevHpDenominator,8192,samplingFrequency);
[lowPassFreqResponse,lowPassFreqAxis]=freqz(chebyshevLpNumerator,chebyshevLpDenominator,8192,samplingFrequency);
notchExactFreqVector=exp(-1j*2*pi*50/samplingFrequency*(0:2));
fprintf('\nChebyshev Specs Verification:\n');
fprintf('HP passband gain at 10 Hz   : %.3f dB\n',20*log10(abs(highPassFreqResponse(find(highPassFreqAxis>=10,1)))));
fprintf('HP stopband atten at 0.1 Hz : %.3f dB\n',20*log10(abs(highPassFreqResponse(find(highPassFreqAxis>=0.1,1)))));
fprintf('LP passband gain at 90 Hz   : %.3f dB\n',20*log10(abs(lowPassFreqResponse(find(lowPassFreqAxis>=90,1)))));
fprintf('LP stopband atten at 150 Hz : %.3f dB\n',20*log10(abs(lowPassFreqResponse(find(lowPassFreqAxis>=150,1)))));
fprintf('Notch atten at 50 Hz        : %.3f dB\n',20*log10(abs(chebyshevNotchNumerator*notchExactFreqVector')/abs(chebyshevNotchDenominator*notchExactFreqVector')));
%% Save
save('cheby_result.mat','ecgCleanedChebyshev','chebyshevHpNumerator','chebyshevHpDenominator','chebyshevNotchNumerator','chebyshevNotchDenominator','chebyshevLpNumerator','chebyshevLpDenominator');
