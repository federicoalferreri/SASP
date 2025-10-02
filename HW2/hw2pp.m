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
zpad=4;
winlen=length(win);
ylen=length(y1);
n_win = floor((ylen-winlen)/hop)+1;
lenframe=winlen*zpad;

% preallocate the stft matrix
Y1 = zeros(lenframe, n_win);
Y2 = zeros(lenframe, n_win);

for i=0:n_win-1
    frame_y1=win.*(y1(i*hop + 1 : i*hop+winlen));
    frame_y2=win.*(y2(i*hop + 1 : i*hop+winlen));

    frame_y1fft=fft(frame_y1, zpad*length(frame_y1));
    frame_y2fft=fft(frame_y2, zpad*length(frame_y2));

    Y1(:,1+i)=frame_y1fft;
    Y2(:,1+i)=frame_y2fft;
end


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
        P(f,t)=1/2*pi*(real(angle(Y2(f,t)/Y1(f,t))));
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
% imagesc([0;1200],[20;120],mask1);
% colormap("gray");
% Apply the binary masks to the STFTs to separate the sources
S1_hat = Y1 .* mask1;
S2_hat = Y1 .* mask2;
S3_hat = Y1 .* mask3;

s1_hat= zeros(ylen*zpad,1);
for i=1:n_win-1
    s1_hat_ifft= real(ifft(S1_hat(:,i)));
    s1_hat(i*hop+1 : i*hop + lenframe) = s1_hat(i*hop+1 : i*hop + lenframe) + s1_hat_ifft;

end

s1_hat=audioplayer(s1_hat, fs);
play(s1_hat); 
