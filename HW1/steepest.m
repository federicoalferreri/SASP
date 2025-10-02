clc
close all
clear

% Load speech signal
[speech, Fs_speech]=audioread("speech.wav");

% Load music signal
[piano, Fs] = audioread("piano.wav");

piano=piano*0.3;
speech=speech*0.3;

%LPC order of speech signal

lpc_speech=round(Fs/1000);

%LPC order of music signal
lpc_piano = round(Fs/1000); 

%number for the padding operation
pad_filters=max(lpc_piano, lpc_speech);

%window for the OLA algorithm

win_size=1024;

hop_size = win_size/4;

window=hanning(win_size);

L = min (length(speech), length(piano)); %output length

output_steep = zeros (L+win_size+pad_filters+1,1);

%number of used windows
n_win = floor((L-win_size)/hop_size)+1;

for i=0:n_win-1

    frame_piano=window.*(piano(i*hop_size + 1 : i*hop_size+win_size)); 
    frame_speech=window.*(speech(i*hop_size + 1 : i*hop_size+win_size));

    frame_piano=padarray(frame_piano, [pad_filters+1, 0], 'post');
    frame_speech=padarray(frame_speech, [pad_filters+1, 0], 'post');

    frame_piano_freq=fft(frame_piano, length(frame_piano)); 
    frame_speech_freq=fft(frame_speech, length(frame_piano));

    [autocorr_piano, lags_piano] = xcorr(frame_piano);

    % autocorr. vector

    r_p= autocorr_piano(lags_piano>=1); %positive lags only: 1,...,M_piano
    r_p=r_p(1:lpc_piano); %autocorrelation values
    
     % Build autocorrelation matrix
    autocorr_piano = autocorr_piano(lags_piano >= 0); %take lags: 0,...,M_piano
    autocorr_piano = autocorr_piano (1:lpc_piano); 
    R_piano = toeplitz(autocorr_piano); %

    %% piano signal steepest algorithm
    eigvalues_piano=eig(R_piano);
    a_p = 0.5; %up to 0.9

    mu_p=a_p*2/max(real(eigvalues_piano)); %step-size factor
     
    w_steep_piano=zeros(lpc_piano, 1);%initial guess: null vector;
   
    n_iter=300; %arbitrary number
    
    for j = 1:n_iter
        w_steep_piano_new=w_steep_piano + mu_p * (r_p - R_piano * w_steep_piano);
    end

    A_steep_piano=[1; -w_steep_piano_new];

    A_steep_piano_pad=padarray(A_steep_piano, [win_size, 0], 'post');
    
    A_steep_piano_freq=fft(A_steep_piano_pad, length(frame_piano)); 

    %error prediction

    e_M_steep_freq = frame_piano_freq .* (A_steep_piano_freq);

    

    %% steepest coefficients: speech

    %compute an estimate of the autocorrelation of the speech signal

    [autocorr_speech, lags_speech] = xcorr(frame_speech);

    %autocorr. vector 

    r_s = autocorr_speech(lags_speech>=1); %positive lags only
    r_s = r_s (1:lpc_speech); %autocorrelation values

    %Build autocorr. matrix

    autocorr_speech = autocorr_speech (lags_speech>=0); %lags >=0 
    autocorr_speech = autocorr_speech (1:lpc_speech);

    R_speech = toeplitz(autocorr_speech); 

    eigvalues_speech=eig(R_speech);
    
    a_s = 0.7; %up to 0.9

    mu_s=a_s*2/max(eigvalues_speech); %step-size factor
     
    w_steep_speech=zeros(lpc_speech, 1);%initial guess: optimum Wiener coeff.;
     
    for j = 1:n_iter
        w_steep_speech_new=w_steep_speech + mu_s * (r_s - R_speech * w_steep_speech);
    end

    H_steep_speech=[1; -w_steep_speech_new];

    H_steep_speech_pad=padarray(H_steep_speech, [win_size,0], 'post');

    H_steep_speech_freq=fft(H_steep_speech_pad, length(frame_speech));

    x_MS_steep_freq= e_M_steep_freq.*(1./(H_steep_speech_freq));

    %Overlap and Add (frame by frame elaboration): summing the output frame by frame
    
    x_MS_steep=real(ifft(x_MS_steep_freq));
    
    output_steep(i*hop_size+1 : i*hop_size+(length(frame_piano))) = output_steep (i*hop_size+1 : i*hop_size + (length(frame_piano))) + x_MS_steep;

end 
output_steep=output_steep(1:length(speech));

%output_steep_play=audioplayer(output_steep, Fs);
%play(output_steep_play);

%write the output on a file

audiowrite('output_steepest.wav', output_steep, Fs);
