clc
close all
clear
%% load mixture signals

% Load mixture signal y1
[y1, fs] = audioread("y1.wav");

% Load mixture signal y2
[y2, ~]=audioread("y2.wav");

%% load true source signals

% Load true source signal s1
[s1, fs_s1] = audioread("s1.wav");

% Load true source signal s2
[s2, ~]=audioread("s2.wav");

% Load true source signal s3
[s3, ~] = audioread("s3.wav");

%% ULA configuration parameters
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

%% parameters for STFT of mixture signals
win = hann(1024); % window
hop = 256; % hop size

zpad=2; %zero padding factor for stft computation


winlen=length(win);
ylen=length(y1);

%number of frames 
n_win = floor((ylen-winlen)/hop)+1;

%chosen frame length 
lenframe=winlen*zpad;

%output signals [vectors initialization]
s1_hat= zeros(ylen*zpad,1);
s2_hat= zeros(ylen*zpad,1);
s3_hat= zeros(ylen*zpad,1);

%binary masks initialization
mask1=zeros(lenframe,n_win);
mask2=zeros(lenframe,n_win);
mask3=zeros(lenframe,n_win);

%% mixture signals spectrogram plot

noverlap=128; %number of samples of overlap for the spectrogram function
nfft = 512; %number of fft points
[y1_sp, y1_f, y1_t] = spectrogram(y1, win, noverlap, nfft, fs); 

figure(1);
subplot(2,1,1);
imagesc(y1_t, y1_f, pow2db(abs(y1_sp)));
colorbar;
axis xy;
xlabel ('time [s]');
ylabel ('frequency [Hz]');
title ('signal source y1');

[y2_sp, y2_f, y2_t] = spectrogram(y2, win, noverlap, nfft, fs);

subplot(2,1,2);
imagesc(y2_t, y2_f, pow2db(abs(y2_sp)));
colorbar;
axis xy;
xlabel ('time [s]');
ylabel ('frequency [Hz]');
title ('signal source y2');

%% true source signals spectrograms 

figure(2);

[s1_sp, s1_f, s1_t] = spectrogram(s1, win, noverlap, nfft, fs_s1);

subplot(3,1,1);
imagesc(s1_t, s1_f, pow2db(abs(s1_sp))); %log amplitude spectrogram plot
colorbar; 
axis xy;
xlabel ('time [s]');
ylabel ('frequency [Hz]');
title ('true source signal s1');

[s2_sp, s2_f, s2_t] = spectrogram(s2, win, noverlap, nfft, fs_s1);

subplot(3,1,2);
imagesc(s2_t, s2_f, pow2db(abs(s2_sp)));
colorbar;
axis xy;
xlabel ('time [s]');
ylabel ('frequency [Hz]');
title ('true source signal s2');

[s3_sp, s3_f, s3_t] = spectrogram(s3, win, noverlap, nfft, fs_s1);

subplot(3,1,3);
imagesc(s3_t, s3_f, pow2db(abs(s3_sp)));
colorbar;
axis xy;
xlabel ('time [s]');
ylabel ('frequency [Hz]');
title ('true source signal s3');

%% stft computation
 
for i=0:n_win-1

    %signal subdivisions in frames
    frame_y1=win.*(y1(i*hop + 1 : i*hop+winlen));
    frame_y2=win.*(y2(i*hop + 1 : i*hop+winlen));

    %frequency domain
    frame_y1fft=fft(frame_y1, lenframe);
    frame_y2fft=fft(frame_y2, lenframe);
    
   %% Compute the features A1, A2, P for each one of the frames
    
    A1=zeros(length(frame_y1fft));
    A2=zeros(length(frame_y1fft));
    P=zeros(length(frame_y1fft));

    % compute the amplitude ratio for each feature vector
    A1=(abs(frame_y1fft)./sqrt((abs(frame_y1fft)).^2+(abs(frame_y2fft)).^2));
    A2=(abs(frame_y2fft)./sqrt((abs(frame_y1fft)).^2+(abs(frame_y2fft)).^2));
   
    % compute the phase difference
    P=(real(angle(frame_y2fft./frame_y1fft)))./(2*pi);
        
   % concatenate amplitude ratio and phase difference vectors into a single feature vector
    
    feature_vector = [A1 A2 P];

    % k-means clustering with k=3
    [idx, C] = kmeans(feature_vector, 3);

    %% Binary masks 
   
    % Create binary masks based on the cluster indices
    mask1(:, i+1) = (idx == 1);
    mask2(:, i+1) = (idx == 2);
    mask3(:, i+1) = (idx == 3);

    % Apply the binary masks to the STFTs to separate the sources
    S1_hat = frame_y1fft .* mask1(:, i+1);
    S2_hat = frame_y1fft .* mask2(:, i+1);
    S3_hat = frame_y1fft .* mask3(:, i+1);

    %% time domain
    s1_hat_ifft= real(ifft(S1_hat));
    s2_hat_ifft= real(ifft(S2_hat));
    s3_hat_ifft= real(ifft(S3_hat));

    %% OLA
    s1_hat(i*hop+1 : i*hop + lenframe) = s1_hat(i*hop+1 : i*hop + lenframe) + s1_hat_ifft;
    s2_hat(i*hop+1 : i*hop + lenframe) = s2_hat(i*hop+1 : i*hop + lenframe) + s2_hat_ifft;
    s3_hat(i*hop+1 : i*hop + lenframe) = s3_hat(i*hop+1 : i*hop + lenframe) + s3_hat_ifft;

end

%% [A1 A2] correlation       
       
% Compute the histogram
 figure(3);
 nbins = n_win; % Number of bins for the histogram
 subplot (2,2,1);
 hist(gca, A1, nbins);
 ylabel('A1 histogram');

 ynbins = 200; % Number of bins for the histogram
 subplot (2,2,4);
 hist(gca, A2, nbins);  
 set(gca, 'view', [90 -90]);
 xlabel('A2 histogram');
 
 % Estimate a continuous pdf from the discrete data
 [pdfx_1, xi_1]= ksdensity(feature_vector(:,1)); %A1
 [pdfy_1, yi_1]= ksdensity(feature_vector(:,2)); %A2
       
 % Create 2-d grid of coordinates and function values, suitable for 3-d plotting
 [xxi_1,yyi_1]     = meshgrid(xi_1,yi_1); 
 [pdfxx_1,pdfyy_1] = meshgrid(pdfx_1,pdfy_1);
        
 % Calculate combined pdf, under assumption of independence
 pdfxy_1 = pdfxx_1.*pdfyy_1;
       
 % Plot the results
 subplot(2,2,3);
 contourf(xxi_1,yyi_1,pdfxy_1);


 %% [A1 P] correlation

 % Compute the histogram
 figure(4);
 nbins = 200; % Number of bins for the histogram
 subplot (2,2,1);
 hist(gca, A1, nbins);
 ylabel('A1 histogram');

 subplot (2,2,4);
 hist(gca, P, nbins);
 set(gca, 'view', [90 -90]);
 xlabel('P histogram');

 % Estimate a continuous pdf from the discrete data 
 [pdfx_2, xi_2]= ksdensity(feature_vector(:,1)); %A1
 [pdfy_2, yi_2]= ksdensity(feature_vector(:,3)); %P
        
 % Create 2-d grid of coordinates and function values, suitable for 3-d plotting
 [xxi_2,yyi_2]     = meshgrid(xi_2,yi_2);
 [pdfxx_2,pdfyy_2] = meshgrid(pdfx_2,pdfy_2);
        
 % Calculate combined pdf, under assumption of independence
 pdfxy_2 = pdfxx_2.*pdfyy_2;
 % Plot the results
 subplot(2,2,3);
 contourf(xxi_2,yyi_2,pdfxy_2);


 %% [A2 P] correlation

 % Compute the histogram
 figure(5);
 nbins = 200; % Number of bins for the histogram
 subplot (2,2,1);
 hist(gca, A2, nbins);
 ylabel('A2 histogram');
 subplot (2,2,4);
 hist(gca, P, nbins);
 set(gca, 'view', [90 -90]);
 xlabel('P histogram');

 % Estimate a continuous pdf from the discrete data 
 [pdfx_3, xi_3]= ksdensity(feature_vector(:,2)); %A2
 [pdfy_3, yi_3]= ksdensity(feature_vector(:,3)); %P
        
 % Create 2-d grid of coordinates and function values, suitable for 3-d plotting
 [xxi_3,yyi_3]     = meshgrid(xi_3,yi_3);
 [pdfxx_3,pdfyy_3] = meshgrid(pdfx_3,pdfy_3);
        
 % Calculate combined pdf, under assumption of independence
 pdfxy_3 = pdfxx_3.*pdfyy_3;
 
 % Plot the results
 subplot(2,2,3);
 contourf(xxi_3,yyi_3,pdfxy_3);


 %%  estimated source signals spectrograms

figure(6);

[s1_hat_sp, s1_hat_f, s1_hat_t] = spectrogram(s1_hat, win, noverlap, nfft, fs);

subplot(3,1,1);
imagesc(s1_hat_t, s1_hat_f, pow2db(abs(s1_hat_sp)));
colorbar;
axis xy;
xlim([0 10]);
xlabel ('time [s]');
ylabel ('frequency [Hz]');
title ('1st estimated source signal');

[s2_hat_sp, s2_hat_f, s2_hat_t] = spectrogram(s2_hat, win, noverlap, nfft, fs);

subplot(3,1,2);
imagesc(s2_hat_t, s2_hat_f, pow2db(abs(s2_hat_sp)));
colorbar;
axis xy;
xlim([0 10]);
xlabel ('time [s]');
ylabel ('frequency [Hz]');
title ('2nd estimated source signal');

[s3_hat_sp, s3_hat_f, s3_hat_t] = spectrogram(s3_hat, win, noverlap, nfft, fs);

subplot(3,1,3);
imagesc(s3_hat_t, s3_hat_f, pow2db(abs(s3_hat_sp)));
colorbar;
axis xy;
xlim([0 10]);
xlabel ('time [s]');
ylabel ('frequency [Hz]');
title ('3rd estimated source signal');


%% black and white mixing masks plot
figure(7);

subplot(3,1,1);
imagesc(mask1, [0 1]);
colormap(gray);
colorbar;
axis xy;
xlabel ('time [frame]');
ylabel ('frequency [Hz]');
title ('Binary Mask 1');

subplot(3,1,2);
imagesc(mask2, [0 1]);
colormap(gray);
colorbar;
axis xy;
xlabel ('time [frame]');
ylabel ('frequency [Hz]');
title ('Binary Mask 2');

subplot(3,1,3);
imagesc(mask3, [0 1]);
colormap(gray);
colorbar;
axis xy;
xlabel ('time [frame]');
ylabel ('frequency [Hz]');
title ('Binary Mask 3');

%% audio play and file creation

s1_play=audioplayer(s1_hat, fs);
%play(s1_play); 

s2_play=audioplayer(s2_hat, fs);
%play(s2_play); 

s3_play=audioplayer(s3_hat, fs);
%play(s3_play); 

%audiowrite('Coletta_Ferreri_s1_hat.wav', s1_hat, fs);
%audiowrite('Coletta_Ferreri_s2_hat.wav', s2_hat, fs);
%audiowrite('Coletta_Ferreri_s3_hat.wav', s3_hat, fs);

