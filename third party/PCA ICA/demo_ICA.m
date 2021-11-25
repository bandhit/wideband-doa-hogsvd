%
% Name:         demo_ICA.m
%               
% Description:  Demonstrates the use of ICA for source separation
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         November 12, 2016
%

%% ICA demo (signals)

rng(42);

% Knobs
n   = 1000;             % # samples
T   = [3, 4, 5];        % # periods for each signal
SNR = 50;               % Signal SNR
d   = 3;                % # mixed observations
r   = 3;                % # independent/principal components

% Generate ground truth
t = @(n,T) linspace(0,1,n) * 2 * pi * T;
Ztrue(1,:) = sin(t(n,T(1)));            % Sinusoid
Ztrue(2,:) = sign(sin(t(n,T(2))));      % Square
Ztrue(3,:) = sawtooth(t(n,T(3)));       % Sawtooth

% Add some noise to make the signals "look" interesting
sigma = @(SNR,X) exp(-SNR / 20) * (norm(X(:)) / sqrt(numel(X)));
Ztrue = Ztrue + sigma(SNR,Ztrue) * randn(size(Ztrue));

% Generate mixed signals
normRows = @(X) bsxfun(@rdivide,X,sum(X,2));
A = normRows(rand(d,3));
Zmixed = A * Ztrue;

% Perform Fast ICA
Zfica = fastICA(Zmixed,r);

% Perform max-kurtosis ICA
Zkica = kICA(Zmixed,r);

% Perform PCA
Zpca = PCA(Zmixed,r);

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cm = hsv(max([3, r, d]));
%cm = linspecer(max([3, r, d]));

figure();

% Truth
subplot(5,1,1);
for i = 1:3
    plot(Ztrue(i,:),'-','Color',cm(i,:)); hold on;
end
title('True signals');
axis tight;

% Observations
subplot(5,1,2);
for i = 1:d
    plot(Zmixed(i,:),'-','Color',cm(i,:)); hold on;
end
title('Observed mixed signals');
axis tight;

% Fast ICA
subplot(5,1,3);
for i = 1:r
    plot(Zfica(i,:),'-','Color',cm(i,:)); hold on;
end
title('Independent components [Fast ICA]');
axis tight;

% Max-kurtosis
subplot(5,1,4);
for i = 1:r
    plot(Zkica(i,:),'-','Color',cm(i,:)); hold on;
end
title('Independent components [max-kurtosis]');
axis tight;

% PCA
subplot(5,1,5);
for i = 1:r
    plot(Zpca(i,:),'-','Color',cm(i,:)); hold on;
end
title('Principal components');
axis tight;
%--------------------------------------------------------------------------

%% ICA demo (audio)
%
% Audio descriptions
%
% source1: siren
% source2: news (stocks)
% source3: foreign 1
% source4: news (food)
% source5: music (classical)
% source6: foreign 2
% source7: music (opera)
% source7: foreign 3
% source9: music (pop)
%
% voice1: talking (linear algebra)
% voice2: talking (sports)
%

rng(42);

% Knobs
Fs = 8000; % Sampling rate of input audio
paths = {'audio/source2.wav','audio/source3.wav','audio/source5.wav'};
%paths = {'audio/voice1.mat','audio/voice2.mat'};
[p, d, r] = deal(numel(paths));
A = randomMixingMatrix(d,p);
%A = [0.6, 0.4; 0.4, 0.6];

% Load audio
Ztrue = loadAudio(paths);

% Generate mixed signals
Zmixed = normalizeAudio(A * Ztrue);

% Fast ICA (negentropy)
Zica1 = normalizeAudio(fastICA(Zmixed,r,'negentropy'));

% Fast ICA (kurtosis)
Zica2 = normalizeAudio(fastICA(Zmixed,r,'kurtosis'));

% Max-kurtosis ICA
Zica3 = normalizeAudio(kICA(Zmixed,r));

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
audio = [];

% True waveforms
for i = 1:p
    audio(end + 1).y = Ztrue(i,:); %#ok
    audio(end).Fs    = Fs;
    audio(end).name  = sprintf('true - %d',i);
end

% Mixed waveforms
for i = 1:d
    audio(end + 1).y = Zmixed(i,:); %#ok
    audio(end).Fs    = Fs;
    audio(end).name  = sprintf('mixed - %d',i);
end

% Fast ICA (negentropy) waveforms
for i = 1:r
    audio(end + 1).y = Zica1(i,:); %#ok
    audio(end).Fs    = Fs;
    audio(end).name  = sprintf('fastICA - neg - %d',i);
end

% Fast ICA (kurtosis) waveforms
for i = 1:r
    audio(end + 1).y = Zica2(i,:); %#ok
    audio(end).Fs    = Fs;
    audio(end).name  = sprintf('fastICA - kur - %d',i);
end

% Max-kurtosis ICA waveforms
for i = 1:r
    audio(end + 1).y = Zica3(i,:); %#ok
    audio(end).Fs    = Fs;
    audio(end).name  = sprintf('kICA - %d',i);
end

% Play audio
PlayAudio(audio);
%--------------------------------------------------------------------------

%% ICA demo (pre-mixed audio)
% TODO(brimoor): how to get this to work?
%
% Audio descriptions
%
% rsm2_mA/B: music + counting english
%  rss_mA/B: counting english + counting spanish
%  rssd_A/B: talking english + talking spanish
%

rng(42);

% Knobs
Fs    = 16000; % Sampling rate of input audio
paths = {'audio/rsm2_mA.wav','audio/rsm2_mB.wav'};
%paths = {'audio/rss_mA.wav','audio/rss_mB.wav'};
%paths = {'audio/rssd_A.wav','audio/rssd_B.wav'};
[d, r] = deal(numel(paths));
r = r + 1;

% Load audio
Zmixed = loadAudio(paths);

% Fast ICA (negentropy)
Zica1 = normalizeAudio(fastICA(Zmixed,r,'negentropy'));

% Fast ICA (kurtosis)
Zica2 = normalizeAudio(fastICA(Zmixed,r,'kurtosis'));

% Max-kurtosis ICA
Zica3 = normalizeAudio(kICA(Zmixed,min(r,d)));

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
audio = [];

% Mixed waveforms
for i = 1:d
    audio(end + 1).y = Zmixed(i,:); %#ok
    audio(end).Fs    = Fs;
    audio(end).name  = sprintf('mixed - %d',i);
end

% Fast ICA (negentropy) waveforms
for i = 1:r
    audio(end + 1).y = Zica1(i,:); %#ok
    audio(end).Fs    = Fs;
    audio(end).name  = sprintf('fastICA - neg - %d',i);
end

% Fast ICA (kurtosis) waveforms
for i = 1:r
    audio(end + 1).y = Zica2(i,:); %#ok
    audio(end).Fs    = Fs;
    audio(end).name  = sprintf('fastICA - kur - %d',i);
end

% Max-kurtosis ICA waveforms
for i = 1:min(r,d)
    audio(end + 1).y = Zica3(i,:); %#ok
    audio(end).Fs    = Fs;
    audio(end).name  = sprintf('kICA - %d',i);
end

% Play audio
PlayAudio(audio);
%--------------------------------------------------------------------------
