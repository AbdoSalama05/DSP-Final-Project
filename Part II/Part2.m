% audio input
audioFile = input('Enter audio file name (default: MJ_BJ.wav): ','s');
if isempty(audioFile)
    audioFile = 'MJ_BJ.wav';
end
[x, fs] = audioread(audioFile);
if size(x,2) == 2
    x = mean(x,2);
end
x = x(:);

% output sample rate
fprintf('\nOutput sample rate settings:\n');
fprintf('1) Use original audio sample rate\n');
fprintf('2) Enter custom sample rate\n');

srMode = input('Choose option (1 or 2): ');

if isempty(srMode) || srMode == 1
    outputFs = fs;
else
    outputFs = input('Enter output sample rate: ');
    if isempty(outputFs)
        outputFs = fs;
    end
end

%mode
fprintf('\nBands settings\n');
fprintf('1) Preset speech bands\n');
fprintf('2) Custom frequency bands\n');

mode = input('Choose option (1 or 2): ');

if isempty(mode)
    mode = 1;
end

if mode == 1
    bands = [
        0 100;
        100 300;
        300 800;
        800 2000;
        2000 5000;
        5000 10000;
        10000 20000
    ];
else
    numBands = input('Number of bands (5–10): ');
    if numBands < 5 || numBands > 10
        numBands = 7;
    end

    bands = zeros(numBands,2);

    for i = 1:numBands
        if i == 1
            fprintf('First band forced to start at 0 Hz (overriding input)\n');
            bands(i,1) = 0;
        else
            bands(i,1) = input('Start frequency (Hz): ');
        end

        if i == numBands
            fprintf('Last band forced to end at 20000 Hz (overriding input)\n');
            bands(i,2) = 20000;
        else
            bands(i,2) = input('End frequency (Hz): ');
        end
    end
end

numBands = size(bands,1);

% gains
gains_dB = zeros(numBands,1);
for k = 1:numBands
    fprintf('Band %d (%d-%d Hz)\n',k,bands(k,1),bands(k,2));
    g = input('Gain dB (default 0): ');
    if isempty(g)
        g = 0;
    end
    gains_dB(k) = g;
end

% filter type
fprintf('\nFilter type:\n');
fprintf('1) FIR (finite impulse response)\n');
fprintf('2) IIR (infinite impulse response)\n');

filterChoice = input('Choose filter type (1 or 2): ');

if isempty(filterChoice)
    filterChoice = 1;
end

%FIR window
if filterChoice == 1
    fprintf('\nFIR type:\n');
    fprintf('1) Hamming\n');
    fprintf('2) Hanning\n');
    fprintf('3) Blackman\n');
    winChoice = input('Choose window type: ');
    if isempty(winChoice)
        winChoice = 1;
    end
end

%IIR window
if filterChoice == 2
    fprintf('\nIIR type\n');
    fprintf('1) Butterworth\n');
    fprintf('2) Chebyshev Type I\n');
    fprintf('3) Chebyshev Type II\n');
    iirChoice = input('Choose IIR type: ');
    if isempty(iirChoice)
        iirChoice = 1;
    end
end

% filter order
N = input('Filter order (default 50 FIR / 4 IIR): ');

if isempty(N)
    if filterChoice == 1
        N = 50;
    else
        N = 4;
    end
end

%safety
if filterChoice == 1 && N < 20
    fprintf('FIR order too small.\n');
    fprintf('Using minimum FIR order = 20\n');
    N = 20;
end
if filterChoice == 2 && N < 2
    fprintf('IIR order too small.\n');
    fprintf('Using minimum IIR order = 2\n');
    N = 2;
end

% process
y = zeros(size(x));
b_all = cell(numBands,1);
a_all = cell(numBands,1);

nyq = fs/2;

for k = 1:numBands

    f1 = bands(k,1);
    f2 = min(bands(k,2), nyq-1);

    Wn = [f1 f2]/nyq;

    Wn = max(Wn,1e-4);
    Wn = min(Wn,1-1e-4);

    if Wn(1) >= Wn(2)
        Wn(2) = Wn(1) + 0.01;
    end

    % filter design
    if filterChoice == 1

        if winChoice == 2
            win = hann(N+1);
        elseif winChoice == 3
            win = blackman(N+1);
        else
            win = hamming(N+1);
        end

        if f1 == 0
            b = fir1(N,Wn(2),'low',win);
        else
            b = fir1(N,Wn,'bandpass',win);
        end

        a = 1;

    else

        if iirChoice == 2
            if f1 == 0
                [b,a] = cheby1(N,1,Wn(2),'low');
            else
                [b,a] = cheby1(N,1,Wn,'bandpass');
            end

        elseif iirChoice == 3
            if f1 == 0
                [b,a] = cheby2(N,40,Wn(2),'low');
            else
                [b,a] = cheby2(N,40,Wn,'bandpass');
            end

        else
            if f1 == 0
                [b,a] = butter(N,Wn(2),'low');
            else
                [b,a] = butter(N,Wn,'bandpass');
            end
        end
    end

    b_all{k} = b;
    a_all{k} = a;

    % filter signal
    yk = filter(b,a,x);

    % gain
    G = 10^(gains_dB(k)/20);
    yk = yk * G;

    y = y + yk;

    fprintf('Band %d used filter order: %d\n',k,N);
end

% normalize
y = y / (max(abs(y)) + 1e-6);

% time domain
t = (0:length(x)-1)/fs;

figure;
subplot(2,1,1); plot(t,x);
title('Original Speech Signal');
xlabel('Time (s)'); ylabel('Amplitude');

subplot(2,1,2); plot(t,y);
title('Equalized Speech Signal (Multi-Band Output)');
xlabel('Time (s)'); ylabel('Amplitude');

% frequency domain
Nfft = length(x);
f = (-Nfft/2:Nfft/2-1)*(fs/Nfft);

X = fftshift(fft(x));
Y = fftshift(fft(y));

figure;
subplot(2,1,1);
plot(f,20*log10(abs(X)+1e-6));
title('Original Spectrum');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');

subplot(2,1,2);
plot(f,20*log10(abs(Y)+1e-6));
title('Equalized Spectrum');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');

%PSD
figure;
[pxx,fpsd] = pwelch(x,[],[],[],fs);
[pyy,~]    = pwelch(y,[],[],[],fs);
plot(fpsd,10*log10(pxx),fpsd,10*log10(pyy));
legend('Original','Equalized');
title('Power Spectral Density Comparison');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

%Spectrogram
figure;
subplot(2,1,1);
spectrogram(x,256,200,256,fs,'yaxis');
title('Original Spectrogram');

subplot(2,1,2);
spectrogram(y,256,200,256,fs,'yaxis');
title('Equalized Spectrogram');

%filter response
for k = 1:numBands
    b = b_all{k};
    a = a_all{k};

    figure;
    freqz(b,a,1024,fs);
    title(sprintf('Band %d Magnitude & Phase Response (%d–%d Hz)',k,bands(k,1),bands(k,2)));

    figure;
    impz(b,a);
    title(sprintf('Band %d Impulse Response (%d–%d Hz)',k,bands(k,1),bands(k,2)));

    % STEP RESPONSE
    figure;
    step = filter(b,a,[1; zeros(100,1)]);
    plot(step);
    title(sprintf('Band %d Step Response (%d–%d Hz)',k,bands(k,1),bands(k,2)));

    figure;
    zplane(b,a);
    title(sprintf('Band %d Pole-Zero Plot (%d–%d Hz)',k,bands(k,1),bands(k,2)));
end

%resampling
[p,q] = rat(outputFs/fs);
y_out = resample(y,p,q);

% 4x
[p,q] = rat(4);
y_4 = resample(y,p,q);

% half
[p,q] = rat(1/2);
y_half = resample(y,p,q);

% normalize outputs
y_out = y_out / (max(abs(y_out)) + 1e-6);
y_4 = y_4 / (max(abs(y_4)) + 1e-6);
y_half = y_half / (max(abs(y_half)) + 1e-6);

%save
audiowrite('out.wav', y_out, outputFs);
audiowrite('out_4fs.wav', y_4, 4*fs);
audiowrite('out_half.wav', y_half, round(fs/2));

sound(y_out, outputFs);
pause(5);

disp('Done.');