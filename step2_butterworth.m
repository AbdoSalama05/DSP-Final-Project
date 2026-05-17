%% Step 2: Butterworth Filter Design and Application
load('ecg_data.mat');
%% Design Filters
[butterworthHpNumerator,butterworthHpDenominator]=butter(4,0.5/(samplingFrequency/2),'high');
normalizedNotchFreq=50/(samplingFrequency/2);

notchBandwidth=normalizedNotchFreq/50;
[notchNumerator,notchDenominator]=iirnotch(normalizedNotchFreq,notchBandwidth);

[butterworthLpNumerator,butterworthLpDenominator]=butter(5,100/(samplingFrequency/2),'low');
%% Frequency Responses
figure;freqz(butterworthHpNumerator,butterworthHpDenominator,1024,samplingFrequency);
title('Butterworth High-Pass (Frequency Response)');
figure;freqz(notchNumerator,notchDenominator,1024,samplingFrequency);
title('Notch Filter 50 Hz (Frequency Response)');
figure;freqz(butterworthLpNumerator,butterworthLpDenominator,1024,samplingFrequency);
title('Butterworth Low-Pass (Frequency Response)');
%% Impulse Responses
impulseResponseLength=500;
impulseTimeAxis=(0:impulseResponseLength-1)/samplingFrequency;
figure;
subplot(3,1,1);plot(impulseTimeAxis,impz(butterworthHpNumerator,butterworthHpDenominator,impulseResponseLength));
title('Impulse Response (Butterworth High-Pass)');xlabel('Time (s)');ylabel('Amplitude');grid on;
subplot(3,1,2);plot(impulseTimeAxis,impz(notchNumerator,notchDenominator,impulseResponseLength));
title('Impulse Response (Notch 50 Hz)');xlabel('Time (s)');ylabel('Amplitude');grid on;
subplot(3,1,3);plot(impulseTimeAxis,impz(butterworthLpNumerator,butterworthLpDenominator,impulseResponseLength));
title('Impulse Response (Butterworth Low-Pass)');xlabel('Time (s)');ylabel('Amplitude');grid on;
%% Step Responses
figure;
subplot(3,1,1);plot(impulseTimeAxis,stepz(butterworthHpNumerator,butterworthHpDenominator,impulseResponseLength));
title('Step Response (Butterworth High-Pass)');xlabel('Time (s)');ylabel('Amplitude');grid on;
subplot(3,1,2);plot(impulseTimeAxis,stepz(notchNumerator,notchDenominator,impulseResponseLength));
title('Step Response (Notch 50 Hz)');xlabel('Time (s)');ylabel('Amplitude');grid on;
subplot(3,1,3);plot(impulseTimeAxis,stepz(butterworthLpNumerator,butterworthLpDenominator,impulseResponseLength));
title('Step Response (Butterworth Low-Pass)');xlabel('Time (s)');ylabel('Amplitude');grid on;
%% Pole-Zero Plots
figure;
subplot(1,3,1);zplane(butterworthHpNumerator,butterworthHpDenominator);title('Pole-Zero (Butterworth HP)');
subplot(1,3,2);zplane(notchNumerator,notchDenominator);title('Pole-Zero (Notch 50 Hz)');
subplot(1,3,3);zplane(butterworthLpNumerator,butterworthLpDenominator);title('Pole-Zero (Butterworth LP)');
%% Apply Filters
ecgAfterHighPass=filtfilt(butterworthHpNumerator,butterworthHpDenominator,noisyEcgSignal);
ecgAfterNotch=filtfilt(notchNumerator,notchDenominator,ecgAfterHighPass);
ecgCleanedButterworth=filtfilt(butterworthLpNumerator,butterworthLpDenominator,ecgAfterNotch);
%% Time-Domain Filtering Steps
figure;
subplot(4,1,1);plot(timeAxis,noisyEcgSignal);title('Noisy ECG');ylabel('mV');grid on;
subplot(4,1,2);plot(timeAxis,ecgAfterHighPass);title('After Butterworth High-Pass');ylabel('mV');grid on;
subplot(4,1,3);plot(timeAxis,ecgAfterNotch);title('After Notch Filter');ylabel('mV');grid on;
subplot(4,1,4);plot(timeAxis,ecgCleanedButterworth);title('After Butterworth Low-Pass');ylabel('mV');xlabel('Time (s)');grid on;
%% PSD: Noisy vs Filtered
[psdNoisySignal,psdFreqAxis]=pwelch(noisyEcgSignal,hamming(512),256,1024,samplingFrequency);
[psdButterworthSignal,~]=pwelch(ecgCleanedButterworth,hamming(512),256,1024,samplingFrequency);
figure;
plot(psdFreqAxis,10*log10(psdNoisySignal),'r');hold on;
plot(psdFreqAxis,10*log10(psdButterworthSignal),'b');
legend('Noisy ECG','Butterworth Filtered');
xlabel('Frequency (Hz)');ylabel('Power/Frequency (dB/Hz)');
title('PSD: Noisy vs Butterworth Filtered');xlim([0 180]);grid on;
xline(0.5,'--m','0.5 Hz');xline(50,'--g','50 Hz');xline(100,'--c','100 Hz');
%% Spectrogram
figure;
subplot(1,2,1);spectrogram(noisyEcgSignal,hamming(256),128,512,samplingFrequency,'yaxis');title('Spectrogram (Noisy ECG)');
subplot(1,2,2);spectrogram(ecgCleanedButterworth,hamming(256),128,512,samplingFrequency,'yaxis');title('Spectrogram (Butterworth Filtered)');
%% SNR
removedNoiseButterworth=noisyEcgSignal-ecgCleanedButterworth;
outputSnrButterworth=10*log10(sum(ecgCleanedButterworth.^2)/sum(removedNoiseButterworth.^2));
fprintf('SNR indicator : %.2f dB\n',outputSnrButterworth);
%% Removed Noise Plot
figure;
plot(timeAxis,removedNoiseButterworth);
title('Removed Noise - Butterworth');xlabel('Time (s)');ylabel('mV');grid on;
%% Print Coefficients
disp('Butterworth HP numerator b:');disp(butterworthHpNumerator);
disp('Butterworth HP denominator a:');disp(butterworthHpDenominator);
disp('Notch numerator b:');disp(notchNumerator);
disp('Notch denominator a:');disp(notchDenominator);
disp('Butterworth LP numerator b:');disp(butterworthLpNumerator);
disp('Butterworth LP denominator a:');disp(butterworthLpDenominator);
%% Verify Specs
[highPassFreqResponse,highPassFreqAxis]=freqz(butterworthHpNumerator,butterworthHpDenominator,8192,samplingFrequency);
[lowPassFreqResponse,lowPassFreqAxis]=freqz(butterworthLpNumerator,butterworthLpDenominator,8192,samplingFrequency);
notchExactFreqVector=exp(-1j*2*pi*50/samplingFrequency*(0:2));
fprintf('\nButterworth Specs Verification:\n');
fprintf('HP passband gain at 10 Hz   : %.3f dB\n',20*log10(abs(highPassFreqResponse(find(highPassFreqAxis>=10,1)))));
fprintf('HP stopband atten at 0.1 Hz : %.3f dB\n',20*log10(abs(highPassFreqResponse(find(highPassFreqAxis>=0.1,1)))));
fprintf('LP passband gain at 90 Hz   : %.3f dB\n',20*log10(abs(lowPassFreqResponse(find(lowPassFreqAxis>=90,1)))));
fprintf('LP stopband atten at 150 Hz : %.3f dB\n',20*log10(abs(lowPassFreqResponse(find(lowPassFreqAxis>=150,1)))));
fprintf('Notch atten at 50 Hz        : %.3f dB\n',20*log10(abs(notchNumerator*notchExactFreqVector')/abs(notchDenominator*notchExactFreqVector')));
%% Save
save('butter_result.mat','ecgCleanedButterworth','butterworthHpNumerator','butterworthHpDenominator','notchNumerator','notchDenominator','butterworthLpNumerator','butterworthLpDenominator');
