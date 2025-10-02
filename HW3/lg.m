%------------------------------------------%
%        *** SSSP - HOMEWORK #3 ***        %
%------------------------------------------%
%     Emulation of the Leslie Speaker      %
%------------------------------------------%
% Name:                                    %
% Student ID:                              %
%------------------------------------------%

clear; close all; clc;

%% modulation speed
mod_speed = 'chorale';

%% Read the input file
[x, Fs] = audioread('HammondRef.wav');
x=x(:,1);           % take left channel only

%% FX parameters 
switch lower(mod_speed)
    
    case {'chorale'}
        freq=2;
            
    case {'tremolo'}
        freq=6;
        
    otherwise
        error('mod_speed \"%s\" not found.', mod_speed)

end

%% Apply FX
N = length(x);
n = 1:N;

% global modulator parameters
alpha = 0.9;
% tremble spectral delay filter parameter 
Ms_t = 0.2;
Mb_t = -0.75;
N_sdf_t = 4;
% bass spectral delay filter parameter 
Ms_b = 0.04;
Mb_b = -0.92;
N_sdf_b = 3;

% cross-over network design
fc = 800; % cutoff frequency

%filter order
N_f = 4;

%Coefficients for the two 4th order butterworth filters
%with cutoff frequency fc
[b_lp, a_lp] = butter(N_f, fc/(Fs/2), 'low'); %LPF design
[b_hp, a_hp] = butter(N_f, fc/(Fs/2), 'high'); %HPF design

% allocate input and output buffers for IIR filters
% hp filter buffers
hpf.state = zeros(N_f,1);
hpf.in = zeros(N_f,1);
% lp filter buffers
lpf.state = zeros(N_f,1);
lpf.in = zeros(N_f,1);
% % treble sdf filter buffers
sdf_t.state = zeros(N_sdf_t + 1, 1);
sdf_t.in = zeros(N_sdf_t + 1, 1);
% % bass sdf filter buffers
sdf_b.state = zeros(N_sdf_b + 1, 1);
sdf_b.in = zeros(N_sdf_b + 1, 1);

% modulators
w_bass = 2*pi*freq;
w_treble = 2*pi*(freq+0.1);
t = n/Fs;

s_bass = sin(w_bass*t);
s_treble = sin(w_treble*t);

% Weighted modulators
m_bass = Ms_b * s_bass + Mb_b; % bass modulator
m_treble = Ms_t * s_treble + Mb_t; % treble modulator

% Initializing outputs
y_hpf = zeros(N, 1);
y_lpf = zeros(N, 1);

y_lp_sdf = zeros(N, 1);
y_hp_sdf = zeros(N, 1);

y_lp_am = zeros(N, 1);
y_hp_am = zeros(N, 1);


%sample processing
for n = 1 : N
    %% HPF
    %y(n) = b(1)*x(n) + z1(n-1);
    y_hpf(n) = b_hp(1)*x(n) + hpf.state(1);
    
    %z1(n) = b(2)*x(n) + z2(n-1) - a(2)y(n)
    hpf.in(1) = b_hp(2)*x(n) + hpf.state(2) - a_hp(2)*y_hpf(n);
    
    %z2(n) = b(3)*x(n) + z3(n-1) - a(3)y(n)
    hpf.in(2) = b_hp(3)*x(n) + hpf.state(3) - a_hp(3)*y_hpf(n);
    
    %z3(n) = b(4)*x(n) + z4(n-1) - a(4)y(n)
    hpf.in(3) = b_hp(4)*x(n) + hpf.state(4) - a_hp(4)*y_hpf(n);
    
    %z4(n) = b(5)*x(n) - a(5)y(n)
    hpf.in(4) = b_hp(5)*x(n) - a_hp(5)*y_hpf(n);
    
    hpf.state = hpf.in;
    
    %% LPF
    y_lpf(n) = b_lp(1)*x(n) + lpf.state(1);
    lpf.in(1) = b_lp(2)*x(n) + lpf.state(2) - a_lp(2)*y_lpf(n);
    lpf.in(2) = b_lp(3)*x(n) + lpf.state(3) - a_lp(3)*y_lpf(n);
    lpf.in(3) = b_lp(4)*x(n) + lpf.state(4) - a_lp(4)*y_lpf(n);
    lpf.in(4) = b_lp(5)*x(n) - a_lp(5)*y_lpf(n);
    
    lpf.state = lpf.in;
    
    %% Treble SDF
    % set last sample to current input x(n)
    sdf_t.in(N_sdf_t + 1) = y_hpf(n);
    
    % calculation of summation
    for i = 0:N_sdf_t
        add = nchoosek(N_sdf_t , i)* m_treble(n)^i * (sdf_t.in(i+1) - sdf_t.state(N_sdf_t+1-i));
        y_hp_sdf(n) = y_hp_sdf(n) + add;
    end
    
    % set last sample to current output y(n)
    sdf_t.state(end) = y_hp_sdf(n);
    % shift buffers to process next samples
    sdf_t.in = circshift(sdf_t.in, -1);
    sdf_t.state = circshift(sdf_t.state, -1);
    % clean last place to compute next sample
    sdf_t.state(end) = 0;
    
    %% Bass SDF
    sdf_b.in(N_sdf_b + 1) = y_lpf(n);
    
    for i = 0:N_sdf_b
        add = nchoosek(N_sdf_b , i)* m_bass(n)^i * (sdf_b.in(i+1) - sdf_b.state(N_sdf_b+1-i));
        y_lp_sdf(n) = y_lp_sdf(n) + add;
    end
    sdf_b.state(end) = y_lp_sdf(n);
    sdf_b.in = circshift(sdf_b.in, -1);
    sdf_b.state = circshift(sdf_b.state, -1);
    sdf_b.state(end) = 0;
    
    %% Modulation
    y_lp_am(n) = (1+alpha * m_bass(n)) * y_lp_sdf(n);
    y_hp_am(n) = (1+alpha * m_treble(n)) * y_hp_sdf(n);
    
    y(n) = y_lp_am(n) + y_hp_am(n);
end



%% Avoid any (possible) clipping
y = rescale(y,-1.,1.);

%% Playback
%audiowrite([mod_speed,'.wav'], y, Fs);
soundsc(y, Fs)

%% Read the reference audio file
dir_name = 'Leslie_ref';
addpath(dir_name);
[y_ref, ~] = audioread(fullfile(dir_name, strcat(mod_speed,'.wav')));

%% Display the MSE
MSE = mean(abs(y.'-y_ref).^2);
MSE_str = sprintf('MSE: %g', MSE);
disp(MSE_str)
