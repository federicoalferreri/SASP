clear all
close all
clc

% Federico Ferreri 10666908
% Emma Coletta 10683341


%% Import Input Audio Signal
[Vin,~] = audioread('ExpSweep.wav');

%% LTSpice Files for Ground-Truth
[OutLowSpice,~]=audioread('outlowsweep.wav');
[OutMidSpice,~]=audioread('outmidsweep.wav');
[OutHighSpice,FsLTSpice]=audioread('outhighsweep.wav');
TsLTSpice=1/FsLTSpice;

%% Sampling frequency (to be varied: FsLTSpice/downSampFact, with downSampFact={4,3,2})
downSampFact=2;
Fs =FsLTSpice/downSampFact; 

%% Downsample Input Signal
Vin=Vin([1:downSampFact:end]);

%% Sampling Period
Ts=1/Fs;
%% Number of Samples
Nsamp=length(Vin);
%% Simulated time
tstop=Nsamp*Ts;
%% Parameters of Dynamic Element
L1=0.35*10^(-3);
L2=0.35*10^(-3);
L3=3.5*10^(-3);
L4=3.5*10^(-3);
C1= 2.8*10^(-6);
C2= 2.8*10^(-6);
C3= 28*10^(-6);
C4= 4.7*10^(-6);
C5=28*10^(-6);
C6=47*10^(-6);
%% Resistive Parameters
R1=10;
RspkLow=8;
R2=10;
RspkMid=8;
RspkHigh=8;
%% WDF setting of free parameters (adaptation conditions)

%Low
Z2_l = (2*L4)/Ts;
Z6_l = Ts/(2*C5);
Z9_l = RspkLow;
Z10_l = R2;
Z11_l = Ts/(2*C6);
Z12_l = Z11_l + Z10_l;
Z8_l = Z12_l;
Z7_l = (Z8_l * Z9_l)/(Z8_l + Z9_l);
Z5_l = Z7_l;
Z4_l = (Z5_l * Z6_l)/(Z5_l + Z6_l);
Z1_l = Z4_l;
Z3_l = Z1_l + Z2_l;

%Mid
Z2_m = (2*L2)/Ts;
Z6_m = Ts/(2*C2);
Z8_m = Ts/(2*C3);
Z12_m = (2*L3)/Ts;
Z15_m = RspkMid;
Z17_m = Ts/(2*C4);
Z16_m = R1;
Z18_m = Z17_m + Z16_m;
Z14_m = Z18_m;
Z13_m = (Z14_m*Z15_m)/(Z14_m+Z15_m);
Z11_m = Z13_m;
Z10_m = (Z11_m*Z12_m)/(Z11_m+Z12_m);
Z7_m = Z10_m;
Z9_m = Z8_m + Z7_m;
Z5_m = Z9_m;
Z4_m = (Z5_m*Z6_m)/(Z5_m+Z6_m);
Z1_m = Z4_m;
Z3_m = Z1_m + Z2_m;

%High
Z2_h = Ts/(2*C1);
Z5_h = RspkHigh;
Z6_h = (2*L1)/Ts;
Z4_h = (Z5_h*Z6_h)/(Z5_h+Z6_h);
Z1_h = Z4_h;
Z3_h = Z1_h + Z2_h;


%% Computation of Scattering Matrices

%Low
%Scattering matrix of the series junction with ports 1,2,3
gamma1_l = Z2_l/(Z2_l+Z1_l);
S1_l = [gamma1_l, (gamma1_l-1), (gamma1_l-1);
       -gamma1_l, (1-gamma1_l), -gamma1_l   ;
       -1  ,   -1     ,      0   ;];
%Scattering matrix of the parallel junction with ports 4,5,6
gamma2_l = Z5_l/(Z5_l+Z6_l);
S2_l = [ 0   , (1-gamma2_l), gamma2_l  ;
         1   ,   -gamma2_l , gamma2_l  ;
         1   , (1-gamma2_l),(gamma2_l-1);];
%Scattering matrix of the parallel junction with ports 7,8,9
gamma3_l = Z8_l/(Z8_l+Z9_l);
S3_l = [ 0   , (1-gamma3_l), gamma3_l  ;
         1   ,   -gamma3_l , gamma3_l  ;
         1   , (1-gamma3_l),(gamma3_l-1);];
%Scattering matrix of the series junction with ports 10,11,12
gamma4_l = Z11_l/(Z11_l+Z10_l);
S4_l = [gamma4_l, (gamma4_l-1), (gamma4_l-1);
       -gamma4_l, (1-gamma4_l), -gamma4_l   ;
       -1  ,   -1     ,      0   ;];

%Mid
%Scattering matrix of the series junction with ports 1,2,3
gamma1_m = Z2_m/(Z2_m+Z1_m);
S1_m = [gamma1_m, (gamma1_m-1), (gamma1_m-1);
       -gamma1_m, (1-gamma1_m), -gamma1_m   ;
       -1  ,   -1     ,      0   ;];
%Scattering matrix of the parallel junction with ports 4,5,6
gamma2_m = Z5_m/(Z5_m+Z6_m);
S2_m = [ 0   , (1-gamma2_m), gamma2_m  ;
         1   ,   -gamma2_m , gamma2_m  ;
         1   , (1-gamma2_m),(gamma2_m-1);];
%Scattering matrix of the series junction with ports 7,8,9
gamma3_m = Z8_m/(Z8_m+Z7_m);
S3_m = [gamma3_m, (gamma3_m-1), (gamma3_m-1);
       -gamma3_m, (1-gamma3_m), -gamma3_m   ;
       -1  ,   -1     ,      0   ;];
%Scattering matrix of the parallel junction with ports 10,11,12
gamma4_m = Z11_m/(Z11_m+Z12_m);
S4_m = [ 0   , (1-gamma4_m), gamma4_m  ;
         1   ,   -gamma4_m , gamma4_m  ;
         1   , (1-gamma4_m),(gamma4_m-1);];
%Scattering matrix of the parallel junction with ports 13,14,15
gamma5_m = Z14_m/(Z14_m+Z15_m);
S5_m = [ 0   , (1-gamma5_m), gamma5_m  ;
         1   ,   -gamma5_m , gamma5_m  ;
         1   , (1-gamma5_m),(gamma5_m-1);];
%Scattering matrix of the series junction with ports 16,17,18
gamma6_m = Z17_m/(Z17_m+Z16_m);
S6_m = [gamma6_m, (gamma6_m-1), (gamma6_m-1);
       -gamma6_m, (1-gamma6_m), -gamma6_m   ;
       -1  ,   -1     ,      0   ;];

%High
%Scattering matrix of the series junction with ports 1,2,3
gamma1_h = Z2_h/(Z2_h+Z1_h);
S1_h = [gamma1_h, (gamma1_h-1), (gamma1_h-1);
       -gamma1_h, (1-gamma1_h), -gamma1_h   ;
       -1  ,   -1     ,      0   ;];
       
%Scattering matrix of the parallel junction with ports 4,5,6
gamma2_h = Z5_h/(Z5_h+Z6_h);
S2_h = [ 0   , (1-gamma2_h), gamma2_h  ;
         1   ,   -gamma2_h , gamma2_h  ;
         1   , (1-gamma2_h),(gamma2_h-1);];

%% Initialization of Waves

%Low
b2_l = 0;
b6_l = 0;
b11_l = 0;
a9_l = 0;
a10_l = 0;

%Mid
b2_m = 0;
b6_m = 0;
b8_m = 0;
b12_m = 0;
b17_m = 0;
a15_m = 0;
a16_m = 0;

%High
b2_h = 0;
b6_h = 0;
a5_h = 0;

%% Initialize Output Signals
% Low
VoutLow=zeros(size(Vin));
% Mid
VoutMid=zeros(size(Vin));
% High
VoutHigh=zeros(size(Vin));

ii=0;
while (ii<Nsamp)
    ii=ii+1;

    %% Manage Dynamic Elements
    
    %Low
    a2_l = -b2_l;
    a6_l = b6_l;
    a11_l = b11_l;
    
    %Mid
    a2_m = -b2_m;
    a6_m = b6_m;
    a8_m = b8_m;
    a12_m = -b12_m;
    a17_m = b17_m;

    %High
    a2_h = b2_h;
    a6_h = -b6_h;
    
    %% Forward Scan
    
     %Low
    b12_l = S4_l(3,:)*[a10_l;a11_l;0];
    a8_l = b12_l;
    b7_l = S3_l(1,:)*[0;a8_l;a9_l];
    a5_l = b7_l;
    b4_l = S2_l(1,:)*[0;a5_l;a6_l];
    a1_l = b4_l;
    b3_l = S1_l(3,:)*[a1_l;a2_l;0];
    
    %Mid
    b18_m = S6_m(3,:)*[a16_m;a17_m;0];
    a14_m = b18_m;
    b13_m = S5_m(1,:)*[0;a14_m;a15_m];
    a11_m = b13_m;
    b10_m = S4_m(1,:)*[0;a11_m;a12_m];
    a7_m = b10_m;
    b9_m = S3_m(3,:)*[a7_m;a8_m;0];
    a5_m = b9_m;
    b4_m = S2_m(1,:)*[0;a5_m;a6_m];
    a1_m = b4_m;
    b3_m = S1_m(3,:)*[a1_m;a2_m;0];
    
    %High
    b4_h = S2_h(1,:)*[0;a5_h;a6_h];
    a1_h = b4_h;
    b3_h = S1_h(3,:)*[a1_h;a2_h;0];

    

    %% Local Root Scattering

    %Low
    a3_l = 2*Vin(ii)- b3_l;
    %Mid
    a3_m = 2*Vin(ii)- b3_m;
    %High
    a3_h = 2*Vin(ii)- b3_h;



    %% Backward Scan

    %Low
    b1_l = S1_l(1,:)*[a1_l;a2_l;a3_l];
    b2_l = S1_l(2,:)*[a1_l;a2_l;a3_l];
    a4_l = b1_l;
    b5_l = S2_l(2,:)*[a4_l;a5_l;a6_l];
    b6_l = S2_l(3,:)*[a4_l;a5_l;a6_l];
    a7_l = b5_l;
    b8_l = S3_l(2,:)*[a7_l;a8_l;a9_l];
    b9_l = S3_l(3,:)*[a7_l;a8_l;a9_l];
    a12_l = b8_l;
    b10_l = S4_l(1,:)*[a10_l;a11_l;a12_l];
    b11_l = S4_l(2,:)*[a10_l;a11_l;a12_l];
    
    %Mid
    b1_m = S1_m(1,:)*[a1_m;a2_m;a3_m];
    b2_m = S1_m(2,:)*[a1_m;a2_m;a3_m];
    a4_m = b1_m;
    b5_m = S2_m(2,:)*[a4_m;a5_m;a6_m];
    b6_m = S2_m(3,:)*[a4_m;a5_m;a6_m];
    a9_m = b5_m;
    b7_m = S3_m(1,:)*[a7_m;a8_m;a9_m];
    b8_m = S3_m(2,:)*[a7_m;a8_m;a9_m];
    a10_m = b7_m;
    b11_m = S4_m(2,:)*[a10_m;a11_m;a12_m];
    b12_m = S4_m(3,:)*[a10_m;a11_m;a12_m];
    a13_m = b11_m;
    b14_m = S5_m(2,:)*[a13_m;a14_m;a15_m];
    b15_m = S5_m(3,:)*[a13_m;a14_m;a15_m];
    a18_m = b14_m;
    b16_m = S6_m(1,:)*[a16_m;a17_m;a18_m];
    b17_m = S6_m(2,:)*[a16_m;a17_m;a18_m];

    %High
    b1_h = S1_h(1,:)*[a1_h;a2_h;a3_h];
    b2_h = S1_h(2,:)*[a1_h;a2_h;a3_h];
    a4_h = b1_h;
    b5_h = S2_h(2,:)*[a4_h;a5_h;a6_h];
    b6_h = S2_h(3,:)*[a4_h;a5_h;a6_h];


    %% Read Output
    
    VoutLow(ii) = -((a9_l + b9_l) / 2);
    VoutMid(ii) = ((a15_m + b15_m) / 2);
    VoutHigh(ii) = -((a5_h + b5_h) / 2);

    
end


%% Output Plots
figure
set(gcf, 'Color', 'w');
subplot(311)
plot(TsLTSpice*[1:length(OutLowSpice)],OutLowSpice,'r','Linewidth',2); hold on;
plot(Ts*[1:Nsamp],VoutLow,'b--','Linewidth',1); grid on; xlim([0,tstop]); 
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{outLow}}$ [V]','Fontsize',16,'interpreter','latex');
legend('LTspice','WDF','Fontsize',16,'interpreter','latex');
title('Output Signals','Fontsize',18,'interpreter','latex');
subplot(312)
plot(TsLTSpice*[1:length(OutMidSpice)],OutMidSpice,'r','Linewidth',2); hold on;
plot(Ts*[1:Nsamp],VoutMid,'b--','Linewidth',1); grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{outMid}}$ [V]','Fontsize',16,'interpreter','latex');
subplot(313)
plot(TsLTSpice*[1:length(OutHighSpice)],OutHighSpice,'r','Linewidth',2); hold on;
plot(Ts*[1:Nsamp],VoutHigh,'b--','Linewidth',1); grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{outHigh}}$ [V]','Fontsize',16,'interpreter','latex');

%% Error Plots
figure
set(gcf, 'Color', 'w');
subplot(311)
plot(Ts*[1:Nsamp],OutLowSpice([1:downSampFact:end])-VoutLow,'k','Linewidth',1);grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$E_{\mathrm{outLow}}$ [V]','Fontsize',16,'interpreter','latex');
title(['Error Signals. $F_{\mathrm{s}}=$ ',num2str(Fs),' Hz'],'Fontsize',18,'interpreter','latex');
subplot(312)
plot(Ts*[1:Nsamp],OutMidSpice([1:downSampFact:end])-VoutMid,'k','Linewidth',1);grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$E_{\mathrm{outMid}}$ [V]','Fontsize',16,'interpreter','latex');
subplot(313)
plot(Ts*[1:Nsamp],OutHighSpice([1:downSampFact:end])-VoutHigh,'k','Linewidth',1);grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$E_{\mathrm{outHigh}}$ [V]','Fontsize',16,'interpreter','latex');

