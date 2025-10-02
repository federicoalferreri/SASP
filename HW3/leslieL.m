function [y,y_lpf,y_hpf,y_hp_sdf] = leslieL(x, Fs, freq)
%Leslie Speaker Emulation
%
% J. Pekonen et al. Computationally Efficient Hammond Organ Synthesis
% in Proc. of the 14th International Conference on Digital Audio
% Effects(DAFx-11), Paris, France, Sept. 19-23, 2011

N=length(x); % length of the input signal
n=1:N;

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
fc=800;     % cutoff frequency
fc_n=fc/(Fs/2); %normalized cutoff frequency

%TODO: compute the coefficients for the two 4th order butterworth filters
Nf=4;       %filter order
%with cutoff frequency fc
[b_lp, a_lp]=butter(Nf, fc_n, 'low'); %LPF design
[b_hp, a_hp]=butter(Nf, fc_n, 'high');  %HPF design

% allocate input and output buffers for IIR filters

% lp filter buffers
lpf.state=zeros(Nf,1);  %Buffer contains last Fs samples of original signal
lpf.in=zeros(Nf,1);     %Buffer contains last FS samples of filtered signal
% hp filter buffers
hpf.state=zeros(Nf,1);  %Buffer contains last Fs samples of original signal
hpf.in=zeros(Nf,1);     %Buffer contains last FS samples of filtered signal

% treble sdf filter buffers
sdf_t.state=zeros(N_sdf_t+1,1);
sdf_t.in=zeros(N_sdf_t+1,1);
% bass sdf filter buffers
sdf_b.state=zeros(N_sdf_b+1,1);
sdf_b.in=zeros(N_sdf_b+1,1);

% modulators
w_bass=2*pi*freq; % bass modulator
w_treble=2*pi*(freq+0.1); % tremble modulator
t=n/Fs;

m_b=sin(w_bass*t);
m_t=sin(w_treble*t);

m_b_p= Ms_b*m_b + Mb_b;     %m'(n)_bass
m_t_p=Ms_t*m_t + Mb_t;      %m'(n)_treble


%Outputs
    y_lpf=zeros(N,1);
    y_hpf=zeros(N,1);
    y_lp_sdf=zeros(N,1);
    y_hp_sdf=zeros(N,1);
    y_lp_am=zeros(N,1);
    y_hp_am=zeros(N,1);
    y(n)=zeros(N,1);


%sample processing
for n=1:N
    
    %% High Pass Filter
    
    %y(n)= b(1)x(n) + z1(n-1)
    y_hpf(n)= b_hp(1)*x(n) + hpf.state(1);  %Non c'è il terzo componente perchè non ho ancora info sui samples passati del segnale filtrato (questo y(n) sarà il valore 1 di hpf.state)
    
    %y(n)=b(n)*x(n) + hpf.state(n) - a(n)*hpf.in(n)
    hpf_in(1)=b_hp(2)*x(n)+ hpf.state(2)-a_hp(2)*y_hpf(n);
    hpf_in(2)=b_hp(3)*x(n)+ hpf.state(3)-a_hp(3)*y_hpf(n);
    hpf_in(3)=b_hp(4)*x(n)+ hpf.state(4)-a_hp(4)*y_hpf(n);
    hpf_in(4)=b_hp(5)*x(n)-a_hp(5)*y_hpf(n);
    
    hpf.state=hpf.in;
    
    
    %% Low Pass Filter
    
    %y(n)= b(1)*x(n) + z1(n-1)(buffer)
    y_lpf(n)= b_lp(1)*x(n) + lpf.state(1);
    
    hpf.in(1)=b_lp(2)*x(n) + lpf.state(2)-a_lp(2)*y_lpf(n);
    hpf.in(2)=b_lp(3)*x(n) + lpf.state(3)-a_lp(3)*y_lpf(n);
    hpf.in(3)=b_lp(4)*x(n) + lpf.state(4)-a_lp(4)*y_lpf(n);
    hpf.in(4)=b_lp(5)*x(n) -a_lp(5)*y_lpf(n);
    
    lpf.state=lpf.in;
    
    %% SDF
    
    % Bass
    
    sdf_b.in(N_sdf_b+1)=y_lpf(n);   %Set the last sample to current x(n) input
    
    for l= 0:N_sdf_b
        %nchoosek computer all the combination of l elements from a group
        %of N_sdf_b elements
        sum= nchoosek(N_sdf_b, l)*m_b_p(n)^l * (sdf_b.in(l+1) - sdf_b.state(N_sdf_b+1-l));
        y_lp_sdf(n)= y_lp_sdf(n) + sum;
        %y(n)= nchoosek(N_sdf_b,l)*m(n)^l*(x(n-(N_sdf_b-l))-y(n-l)
    end
    
    sdf_b.state(end)=y_lp_sdf(n);
    sdf_b.in=circshift(sdf_b.in, -1);  %sposto il primo valore in ultima posizione
    sdf_b.state=circshift(sdf_b.state, -1);
    sdf_b.state(end)=0;
    
    
    % Treble
    
    sdf_t.in(N_sdf_t+1)=y_hpf(n);   %Set the last sample to current x(n) input
    
    for l= 0:N_sdf_t
        %nchoosek computer all the combination of l elements from a group
        %of N_sdf_b elements
        sum= nchoosek(N_sdf_t, l)*m_t_p(n)^l * (sdf_t.in(l+1) - sdf_t.state(N_sdf_t+1-l));
        y_hp_sdf(n)= y_hp_sdf(n) + sum;
        %y(n)= nchoosek(N_sdf_b,l)*m(n)^l*(x(n-(N_sdf_b-l))-y(n-l)
    end
    
    sdf_t.state(end)=y_hp_sdf(n);
    sdf_t.in=circshift(sdf_t.in, -1);  %sposto il primo valore in ultima posizione
    sdf_t.state=circshift(sdf_t.state, -1);
    sdf_t.state(end)=0;
    
    
    %% Modulation (AM)
    
    y_lp_am(n)= (1+alpha*m_b_p(n)).*y_lp_sdf(n);
    y_hp_am(n)= (1+alpha*m_t_p(n)).*y_hp_sdf(n);
    
    %Sum signals
    y(n)=y_lp_am(n) + y_hp_am(n);
    
    
end
end

