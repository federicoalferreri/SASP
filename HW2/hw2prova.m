clc
close all
clear

% Load mixture signal y1
[y1, fs] = audioread("y1.wav");

% Load mixture signal y2
[y2, ~]=audioread("y2.wav");


% Define ULA configuration
d = 0.09; % distance between microphones in meters
N = 3; %number of sources;
M = 2;% number of microphones
pos = (0:N-1)*d; % microphone positions in meters
c= 340;
% Define DOAs and distances of sources
theta1 = pi/6; % DOA of source 1 in radians
theta2 = 17*pi/36; % DOA of source 2 in radians 85°
theta3 = -2*pi/9; % DOA of source 3 in radians 40°
d1 = 0.75; % distance of source 1 from reference microphone in meters
d2 = 0.75; % distance of source 2 from reference microphone in meters
d3 = 0.75; % distance of source 3 from reference microphone in meters

% Compute STFT of mixture signals
win = hann(1024); % window
hop = 256; % hop size
nfft=4096;
y1=y1(:);
y2=y2(:);
winlen=length(win);
y1len=length(y1);
y2len=length(y2);

% stft matrix size estimation and preallocation
NUP = ceil((1+nfft)/2);     % number of ffts
L = 1+fix((y1len-winlen)/hop); % number of frames
Y1 = zeros(NUP, L);% preallocate the stft matrix
Y2 = zeros(NUP, L);
% STFT (via time-localized FFT)
for l = 0:L-1
    % windowing
    y1w = y1(1+l*hop : winlen+l*hop).*win;
    y2w = y2(1+l*hop : winlen+l*hop).*win;
   
    % FFT
     Y1fft= fft(y1w, nfft);
     Y2fft= fft(y2w, nfft);
    
    % update of the stft matrix
    Y1(:, 1+l) = Y1fft(1:NUP);
    Y2(:, 1+l) = Y2fft(1:NUP);
end
% calculation of the time and frequency vectors
% t = (wlen/2:hop:wlen/2+(L-1)*hop)/fs;
% f = (0:NUP-1)*fs/nfft;


% Y1 = stft(y1,fs, Window=win, OverlapLength=hop);
% Y2 = stft(y2,fs, Window=win, OverlapLength=hop);

% Compute the amplitude ratio A and time delay T for each time-frequency bin
freq = size(Y1,1);
time = size(Y1,2);
A1=zeros(freq, time);
A2=zeros(freq, time);
P=zeros(freq, time);
for f = 1:freq
    for t = 1:time
        % Compute the amplitude ratio for each feature vector
        A1(f,t)=(abs(Y1(f,t))/sqrt(abs(Y1(f,t))^2+abs(Y2(f,t))^2));
        A2(f,t)=(abs(Y2(f,t))/sqrt(abs(Y1(f,t))^2+abs(Y2(f,t))^2));
        % Compute the time delay for feature vector
        P(f,t)=1/2*pi*(angle(Y2(f,t)/Y1(f,t)));
    end
end
scatter(A1,A2);
colorbar;
% Concatenate amplitude ratio and time delay vectors into a single feature vector
Avec = reshape([A1; A2; P], [3, freq*time]);

% Run k-means clustering with k=3
[idx, C] = kmeans(Avec', 3);

% Reshape the cluster indices into the time-frequency shape
idx = reshape(idx, [freq, time]);

% Create binary masks based on the cluster indices
mask1 = (idx == 1);
mask2 = (idx == 2);
mask3 = (idx == 3);
% imagesc(mask1);
% colormap("gray");
% Apply the binary masks to the STFTs to separate the sources
S1_hat = Y1 .* mask1;
S2_hat = Y1 .* mask2;
S3_hat = Y1 .* mask3;

% Compute the time-domain signals of the estimated sources
% y1_hat = istft(Y1_hat, fs, Window=win, OverlapLength=hop);
% y2_hat = istft(Y2_hat, fs, Window=win, OverlapLength=hop);
% y3_hat = istft(Y3_hat, fs, Window=win, OverlapLength=hop);

s1len = winlen + (L-1)*hop;    % estimate the length of the signal vector
s1 = zeros(1, s1len);         % preallocate the signal vector
% reconstruction of the whole spectrum
if rem(nfft, 2)             
    % odd nfft excludes Nyquist point
    S1 = [S1_hat; conj(flipud(S1_hat(2:end, :)))];
else                        
    % even nfft includes Nyquist point
    S1 = [S1_hat; conj(flipud(S1_hat(2:end-1, :)))];
end
% columnwise IFFT on the STFT-matrix
s1w = real(ifft(S1));
s1w = s1w(1:winlen, :);
% Weighted-OLA
for l = 1:L
    s1(1+(l-1)*hop : winlen+(l-1)*hop) = s1(1+(l-1)*hop : winlen+(l-1)*hop) + ...
                                      (s1w(:, l).*win)';
end
% scaling of the signal
W0 = sum(win.*win);                  
s1 = s1.*hop/W0;   

s1=audioplayer(s1, fs);
play(s1);
