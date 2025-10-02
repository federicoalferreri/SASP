clc
close all
clear

% Load mixture signal y1
[y1, fs] = audioread("y1.wav");

% Load mixture signal y2
[y2, ~]=audioread("y2.wav");

y=[y1 y2];

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
Y1 = stft(y1,fs, Window=win, OverlapLength=hop);
Y2 = stft(y2,fs, Window=win, OverlapLength=hop);
theta=[theta1 theta2 theta3];

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

% Apply the binary masks to the STFTs to separate the sources
Y1_hat = Y1 .* mask1;
Y2_hat = Y1 .* mask2;
Y3_hat = Y1 .* mask3;

% Compute the time-domain signals of the estimated sources
y1_hat = istft(Y1_hat, fs, Window=win, OverlapLength=hop);
y2_hat = istft(Y2_hat, fs, Window=win, OverlapLength=hop);
y3_hat = istft(Y3_hat, fs, Window=win, OverlapLength=hop);

y3_hat=audioplayer(y3_hat, fs);
play(y3_hat);
