function [y,y_lpf,y_hpf,y_hp_sdf] = leslie(x, Fs, freq)
%Leslie Speaker Emulation
%
% J. Pekonen et al. Computationally Efficient Hammond Organ Synthesis
% in Proc. of the 14th International Conference on Digital Audio
% Effects(DAFx-11), Paris, France, Sept. 19-23, 2011

% length of the input signal
N = length(x);

% global modulator parameters
alpha=0.9;

% tremble spectral delay filter parameter 
Ms_t=0.2;
Mb_t=-0.75;
N_sdf_t=4;

% bass spectral delay filter parameter 
Ms_b=0.04;
Mb_b=-0.92;
N_sdf_b=3;

% cross-over network design
fc=800;                 % cutoff frequency

%TODO: compute the coefficients for the two 4th order butterworth filters
%with cutoff frequency fc
filter_order = 4;
[b_lp, a_lp]=butter(filter_order,fc/(Fs/2),'low');  %LPF design
[b_hp, a_hp]=butter(filter_order, fc/(Fs/2), 'high');  %HPF design

% allocate input and output buffers for IIR filters

% hp filter buffers
hpf.state=zeros(filter_order,1); %buffer to perform filtering
hpf.in=zeros(filter_order,1); %we limit the # of input samples to the filter order.

% lp filter buffers
lpf.state=zeros(filter_order,1); 
lpf.in=zeros(filter_order,1);

% treble sdf filter buffers
sdf_t.state=zeros(N_sdf_t+1,1);
sdf_t.in=zeros(N_sdf_t+1,1);

% bass sdf filter buffers
sdf_b.state=zeros(N_sdf_b+1, 1);
sdf_b.in=zeros(N_sdf_b+1, 1);

% Initializing outputs
y_lpf = zeros(N, 1);
y_hpf = zeros(N, 1);

y_lp_sdf = zeros(N, 1);
y_hp_sdf = zeros(N, 1);

y_lp_am = zeros(N, 1);
y_hp_am = zeros(N, 1);


%% sinusoidal signals for the modulators

j=1:N;
t=j/Fs;

%difference in speeds of the different motors expressed in terms of
%frequency of the sinusoid describing the motion
freq_bass=freq;
freq_treble=freq+0.1;

motor_bass=sin(2*pi*freq_bass*t);
motor_treble=sin(2*pi*freq_treble*t);

% weighted modulators
m_b=Ms_b*motor_bass+Mb_b; % bass modulator
m_t=Ms_t*motor_treble+Mb_t; % tremble modulator

%sample processing
for n=1:N

    % compute crossover network filters outputs. 
    % time domain explicit filtering computed manually for the 4 coeff. 
    
    %% low pass filtering: direct form II

    %y(n)=b1*x(n) + buffer 
    y_lpf(n)= b_lp(1)*x(n) + lpf.state(1);

    lpf.in(1)= b_lp(2)*x(n) - a_lp(2)*y_lpf(n) + lpf.state(2);

    lpf.in(2)= b_lp(3)*x(n) - a_lp(3)*y_lpf(n) + lpf.state(3);

    lpf.in(3) = b_lp(4)*x(n) - a_lp(4)*y_lpf(n) + lpf.state(4);

    lpf.in(4) = b_lp(5)*x(n) - a_lp(5)*y_lpf(n);

    %update buffer: var lp.state; contains past samples of the input signal
    
    lpf.state=lpf.in;

    %% high pass filtering
    
    %y(n)=b1*x(n) + buffer 
    y_hpf(n) = b_hp(1)*x(n) + hpf.state(1);

    hpf.in(1)=b_hp(2)*x(n) - a_hp(2)*y_hpf(n) + hpf.state(2);

    hpf.in(2)=b_hp(3)*x(n) - a_hp(3)*y_hpf(n) + hpf.state(3);

    hpf.in(3)=b_hp(4)*x(n) - a_hp(4)*y_hpf(n) + hpf.state(4);

    hpf.in(4)=b_hp(5)*x(n) - a_hp(5)*y_hpf(n);

    %update buffer: var hp.state; contains past samples of the input signal
    hpf.state=hpf.in;

    %% SDF 

    % compute bass SDF output
    %y_lp_sdf=...
    sdf_b.in(N_sdf_b + 1) = y_lpf(n);
    
    for i = 0:N_sdf_b
        add = nchoosek(N_sdf_b , i)* m_b(n)^i * (sdf_b.in(i+1) - sdf_b.state(N_sdf_b+1-i));
        y_lp_sdf(n) = y_lp_sdf(n) + add;
    end
    sdf_b.state(end) = y_lp_sdf(n);
    sdf_b.in = circshift(sdf_b.in, -1);
    sdf_b.state = circshift(sdf_b.state, -1);
    sdf_b.state(end) = 0;

    % compute treble SDF output
    %y_hp_sdf=...
    
    % set last sample to current input x(n)
    sdf_t.in(N_sdf_t + 1) = y_hpf(n);
    
    % calculation of summation
    for i = 0:N_sdf_t
        add = nchoosek(N_sdf_t , i)* m_t(n)^i * (sdf_t.in(i+1) - sdf_t.state(N_sdf_t+1-i));
        y_hp_sdf(n) = y_hp_sdf(n) + add;
    end
    
    % set last sample to current output y(n)
    sdf_t.state(end) = y_hp_sdf(n);
    % shift buffers to process next samples
    sdf_t.in = circshift(sdf_t.in, -1);
    sdf_t.state = circshift(sdf_t.state, -1);
    % clean last place to compute next sample
    sdf_t.state(end) = 0;
    
    %% Amplitude Modulation 
    
    % implement AM modulation blocks
    y_lp_am(n)= (1+alpha*m_b(n))*y_lp_sdf(n);
    y_hp_am(n)= (1+alpha*m_t(n))*y_hp_sdf(n);

    y(n)= y_lp_am(n)+y_hp_am(n);

end

end

