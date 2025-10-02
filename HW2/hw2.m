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
Y= cat(2, Y1, Y2);


% Compute steering vectors
theta = linspace(-pi/2, pi/2, 180); % angles in radians
A1 = exp(1j*2*pi*pos(1)*(sin(theta-theta1)/c)*fs); % steering vector for source 1
A2 = exp(1j*2*pi*pos(1)*(sin(theta-theta2)/c)*fs); % steering vector for source 2
A3 = exp(1j*2*pi*pos(1)*(sin(theta-theta3)/c)*fs); % steering vector for source 3

% Compute spatial covariance matrix R
F = size(Y, 1);
M = size(Y, 2);
R = zeros(2, 2, size(Y, 2), size(Y, 3)); % initialize R with the correct dimensions

for m = 1:size(Y, 3)
    for f = 1:size(Y, 2)
        R(:, :, f, m) = Y(:, f, m)*Y(:, f, m)';
    end
end

% Permute the first two dimensions of R
R = permute(R, [1, 2, 4, 3]);

% Compute k-means clustering
K = 3; % number of clusters
[~, idx] = kmeans(R(:,:), K);

% Reshape the cluster indices into a matrix
idx_mat = reshape(idx, M, N);

% Compute binary masks
mask1 = (idx_mat == 1);
mask2 = (idx_mat == 2);
mask3 = (idx_mat == 3);


% Compute separated signals
S1 = Y1.*repmat(mask1, [F, 1, 1]);
S2 = Y1.*repmat(mask2, [F, 1, 1]);
S3 = Y1.*repmat(mask3, [F, 1, 1]);

% Compute inverse STFT to obtain separated signals
s1 = istft(S1, win, hop);
s2 = istft(S2, win, hop);
s3 = istft(S3, win, hop);


s1=audioplayer(s1, fs);
play(s1);




