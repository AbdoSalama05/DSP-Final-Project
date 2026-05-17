%% Step 4: FIR Filter Design and Application
load('ecg_data.mat');
%% Design Filters
firOrder=500;
firHighPassNumerator=fir1(firOrder,0.5/(samplingFrequency/2),'high');firHighPassDenominator=1;

firNotchNumerator=fir1(firOrder,[49.5 50.5]/(samplingFrequency/2),'stop');firNotchDenominator=1;

firLowPassNumerator=fir1(firOrder,100/(samplingFrequency/2),'low');firLowPassDenominator=1;
%% Frequency Responses
figure;freqz(firHighPassNumerator,1,1024,samplingFrequency);title('FIR High-Pass (Frequency Response)');
figure;freqz(firNotchNumerator,1,1024,samplingFrequency);title('FIR Notch 50 Hz (Frequency Response)');
figure;freqz(firLowPassNumerator,1,1024,samplingFrequency);title('FIR Low-Pass (Frequency Response)');
%% Impulse Responses
impulseResponseLength=500;
impulseTimeAxis=(0:impulseResponseLength-1)/samplingFrequency;
figure;
subplot(3,1,1);plot(impulseTimeAxis,impz(firHighPassNumerator,1,impulseResponseLength));
title('Impulse Response (FIR High-Pass)');xlabel('Time (s)');ylabel('Amplitude');grid on;
subplot(3,1,2);plot(impulseTimeAxis,impz(firNotchNumerator,1,impulseResponseLength));
title('Impulse Response (FIR Notch 50 Hz)');xlabel('Time (s)');ylabel('Amplitude');grid on;
subplot(3,1,3);plot(impulseTimeAxis,impz(firLowPassNumerator,1,impulseResponseLength));
title('Impulse Response (FIR Low-Pass)');xlabel('Time (s)');ylabel('Amplitude');grid on;
%% Step Responses
figure;
subplot(3,1,1);plot(impulseTimeAxis,stepz(firHighPassNumerator,1,impulseResponseLength));
title('Step Response (FIR High-Pass)');xlabel('Time (s)');ylabel('Amplitude');grid on;
subplot(3,1,2);plot(impulseTimeAxis,stepz(firNotchNumerator,1,impulseResponseLength));
title('Step Response (FIR Notch 50 Hz)');xlabel('Time (s)');ylabel('Amplitude');grid on;
subplot(3,1,3);plot(impulseTimeAxis,stepz(firLowPassNumerator,1,impulseResponseLength));
title('Step Response (FIR Low-Pass)');xlabel('Time (s)');ylabel('Amplitude');grid on;
%% Pole-Zero Plots
figure;
subplot(1,3,1);zplane(firHighPassNumerator,1);title('Pole-Zero (FIR HP)');
subplot(1,3,2);zplane(firNotchNumerator,1);title('Pole-Zero (FIR Notch)');
subplot(1,3,3);zplane(firLowPassNumerator,1);title('Pole-Zero (FIR LP)');
%% Apply Filters
ecgAfterHighPass=filtfilt(firHighPassNumerator,1,noisyEcgSignal);
ecgAfterNotch=filtfilt(firNotchNumerator,1,ecgAfterHighPass);
ecgCleanedFir=filtfilt(firLowPassNumerator,1,ecgAfterNotch);
%% Time-Domain Filtering Steps
figure;
subplot(4,1,1);plot(timeAxis,noisyEcgSignal);title('Noisy ECG');ylabel('mV');grid on;
subplot(4,1,2);plot(timeAxis,ecgAfterHighPass);title('After FIR High-Pass');ylabel('mV');grid on;
subplot(4,1,3);plot(timeAxis,ecgAfterNotch);title('After FIR Notch');ylabel('mV');grid on;
subplot(4,1,4);plot(timeAxis,ecgCleanedFir);title('After FIR Low-Pass');ylabel('mV');xlabel('Time (s)');grid on;
%% PSD: Noisy vs Filtered
[psdNoisySignal,psdFreqAxis]=pwelch(noisyEcgSignal,hamming(512),256,1024,samplingFrequency);
[psdFirSignal,~]=pwelch(ecgCleanedFir,hamming(512),256,1024,samplingFrequency);
figure;
plot(psdFreqAxis,10*log10(psdNoisySignal),'r');hold on;
plot(psdFreqAxis,10*log10(psdFirSignal),'b');
legend('Noisy ECG','FIR Filtered');
xlabel('Frequency (Hz)');ylabel('Power/Frequency (dB/Hz)');
title('PSD: Noisy vs FIR Filtered');xlim([0 180]);grid on;
xline(0.5,'--m','0.5 Hz');xline(50,'--g','50 Hz');xline(100,'--c','100 Hz');
%% Spectrogram
figure;
subplot(1,2,1);spectrogram(noisyEcgSignal,hamming(256),128,512,samplingFrequency,'yaxis');title('Spectrogram (Noisy ECG)');
subplot(1,2,2);spectrogram(ecgCleanedFir,hamming(256),128,512,samplingFrequency,'yaxis');title('Spectrogram (FIR Filtered)');
%% SNR
removedNoiseFir=noisyEcgSignal-ecgCleanedFir;
outputSnrFir=10*log10(sum(ecgCleanedFir.^2)/sum(removedNoiseFir.^2));
fprintf('SNR indicator : %.2f dB\n',outputSnrFir);
%% Removed Noise Plot
figure;
plot(timeAxis,removedNoiseFir);
title('Removed Noise - FIR');xlabel('Time (s)');ylabel('mV');grid on;
%% Print Coefficients
disp('FIR HP first 10 b coefficients:');disp(firHighPassNumerator(1:10));
disp('FIR Notch first 10 b coefficients:');disp(firNotchNumerator(1:10));
disp('FIR LP first 10 b coefficients:');disp(firLowPassNumerator(1:10));
%% Verify Specs
firSampleIndices=0:firOrder;
fprintf('\nFIR Specs Verification:\n');
fprintf('HP passband gain at 10 Hz   : %.3f dB\n',20*log10(abs(firHighPassNumerator*exp(-1j*2*pi*10/samplingFrequency*firSampleIndices)')));
fprintf('HP stopband atten at 0.1 Hz : %.3f dB\n',20*log10(abs(firHighPassNumerator*exp(-1j*2*pi*0.1/samplingFrequency*firSampleIndices)')));
fprintf('LP passband gain at 90 Hz   : %.3f dB\n',20*log10(abs(firLowPassNumerator*exp(-1j*2*pi*90/samplingFrequency*firSampleIndices)')));
fprintf('LP stopband atten at 150 Hz : %.3f dB\n',20*log10(abs(firLowPassNumerator*exp(-1j*2*pi*150/samplingFrequency*firSampleIndices)')));
fprintf('Notch atten at 50 Hz        : %.3f dB\n',20*log10(abs(firNotchNumerator*exp(-1j*2*pi*50/samplingFrequency*firSampleIndices)')));
%% Save
save('fir_result.mat','ecgCleanedFir','firHighPassNumerator','firHighPassDenominator','firNotchNumerator','firNotchDenominator','firLowPassNumerator','firLowPassDenominator');
