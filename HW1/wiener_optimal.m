clc
close all
clear

% Load music signal
[piano, Fs] = audioread("piano.wav");

% Load speech signal
[speech, ~]=audioread("speech.wav");

%LPC order of speech signal
lpc_speech = round(Fs/1000); 

%LPC order of music signal
lpc_piano = round(Fs/1000);  

%number for the padding operation
pad_filters=max(lpc_piano, lpc_speech);

%window for the OLA algorithm

win_len = 1024; %arbitrary value

%speech and music have the same duration, output will hypotetically last
%the same amount of time

L = min (length(speech), length(piano)); %output length
output = zeros (L+win_len+pad_filters+1,1);
output_steep = zeros (L,1);

%COLA condition: hop size = window size/M

hop_size = win_len/4;

%Hann windowing 

window = hann(win_len);

%number of windows
n_win = floor((L-win_len)/hop_size)+1;

%For each frame of the signals: 
for i=0:n_win-1
    %fragmentation of both piano and speech signals

    frame_piano=window.*(piano(i*hop_size + 1 : i*hop_size+win_len));
    frame_speech=window.*(speech(i*hop_size + 1 : i*hop_size+win_len));
    
    frame_piano=padarray(frame_piano, [pad_filters+1, 0], 'post');
    frame_speech=padarray(frame_speech, [pad_filters+1, 0], 'post');
    
    %frequency domain signal frames 
    frame_piano_freq=fft(frame_piano, length(frame_piano)); 
    frame_speech_freq=fft(frame_speech, length(frame_piano));

    %% Music Signal 

    %compute an estimate of the autocorrelation of the music signal

    [autocorr_piano, lags_piano] = xcorr(frame_piano);

    % autocorr. vector

    r_p= autocorr_piano(lags_piano>=1); %positive lags only: 1,...,lpc_piano
    r_p=r_p(1:lpc_piano); % autocorrelation values

    % Build autocorrelation matrix
    autocorr_piano = autocorr_piano(lags_piano >= 0); %take lags: 0,...,lpc_piano
    autocorr_piano = autocorr_piano(1:lpc_piano); 
    R_piano = toeplitz(autocorr_piano); %autocorrelation matrix using the toeplitz definition

    %Wiener-Hopf closed form solution

    WH_piano = R_piano\r_p; %or inv(R_piano)*r_p;

    A_piano = [1;-WH_piano]; %inverse filter coeff, used in the whitening filter

    %zero padded fft of inverse filter coefficients
    A_piano=padarray(A_piano, [win_len, 0], 'post');
    A_piano_freq=fft(A_piano, length(frame_piano)); 

    %Whitening filter in the frequency domain of the piano signal 

    e_M_freq = frame_piano_freq.*A_piano_freq; %error prediction

    %% Speech Signal - same as piano

    %compute an estimate of the autocorrelation of the speech signal

    [autocorr_speech, lags_speech] = xcorr(frame_speech);

    %autocorr. vector 

    r_s = autocorr_speech(lags_speech>=1); %positive lags only
    r_s = r_s (1:lpc_speech); %autocorrelation values

    %Build autocorr. matrix

    autocorr_speech = autocorr_speech (lags_speech>=0); %lags >=0 
    autocorr_speech = autocorr_speech (1:lpc_speech);

    R_speech = toeplitz(autocorr_speech); 

    %Wiener-Hopf closed form solution

    WH_speech = R_speech\r_s; %or inv(R_speech)*rs

    H_speech = [1;-WH_speech];

    H_speech=padarray(H_speech, [win_len,0], 'post');

    H_speech_freq=fft(H_speech, length(frame_piano));
    
    %Shaping filter of the error prediction
    x_MS_freq = e_M_freq.*(1./H_speech_freq); 

    %Overlap and Add (frame by frame elaboration): summing the output frame by frame
    
    x_MS=real(ifft(x_MS_freq)); %back to the time domain to perform OLA
    
   output(i*hop_size+1 : i*hop_size + (length(frame_piano))) = output (i*hop_size+1 : i*hop_size + (length(frame_piano))) + x_MS;  

     
end 


output=output(1:length(speech));

%player object
output_play=audioplayer(output, Fs);

%play
play(output_play);

%write the output on a file
%r_p and r_s: p_piano and p_speech in the slides

%audiowrite('output_wiener.wav', output, Fs);
