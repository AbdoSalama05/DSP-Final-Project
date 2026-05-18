%% Step 1: Load ECG Signal and Add Noise
addpath('C:\Users\PolaNasser\Downloads\wfdb-app-toolbox-0-10-0\mcode')
[rawSignal,samplingFrequency,timeVector]=rdsamp('100',[],3600);
samplingFrequency
rawEcgSignal=rawSignal(:,1);
timeAxis=timeVector(:);
rawEcgSignal=rawEcgSignal(:);
%% Plot Raw ECG
figure;
plot(timeAxis,rawEcgSignal);
xlabel('Time (s)');ylabel('mV');
title('Raw ECG Signal');
grid on;
%% Add Synthetic Noise
baselineWanderNoise=0.3*sin(2*pi*0.3*timeAxis);
powerlineInterferenceNoise=0.2*sin(2*pi*50*timeAxis);
[bEmg,aEmg]=butter(4,[20 150]/(samplingFrequency/2),'bandpass');
emgMusclNoise=filtfilt(bEmg,aEmg,0.05*randn(size(timeAxis)));
noisyEcgSignal=rawEcgSignal+baselineWanderNoise+powerlineInterferenceNoise+emgMusclNoise;
%% Plot Raw vs Noisy
figure;
subplot(2,1,1);plot(timeAxis,rawEcgSignal);
title('Raw ECG');xlabel('Time (s)');ylabel('mV');grid on;
subplot(2,1,2);plot(timeAxis,noisyEcgSignal);
title('Noisy ECG');xlabel('Time (s)');ylabel('mV');grid on;
%% PSD of Noisy Signal
[powerSpectralDensity,frequencyAxis]=pwelch(noisyEcgSignal,[],[],[],samplingFrequency);
figure;
semilogy(frequencyAxis,powerSpectralDensity);
xlabel('Frequency (Hz)');ylabel('PSD');
title('PSD of Noisy ECG');grid on;
xline(0.5,'r--','Baseline 0.5 Hz');
xline(50,'g--','Powerline 50 Hz');
xline(100,'b--','EMG 100 Hz');
%% Save
save('ecg_data.mat','rawEcgSignal','noisyEcgSignal','timeAxis','samplingFrequency');