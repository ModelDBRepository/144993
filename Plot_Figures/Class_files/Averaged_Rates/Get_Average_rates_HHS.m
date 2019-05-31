% Simulation of the HH model finds the average rates data and saves it

clear all;
close all;
clc

global dt A phi_HH

phi_HH=2;
I0=[10]; %[microamper]/cm^2 - most of the script below does not work for arrays larger than 1
L_I=length(I0);
Time=50; %[msec] length of simulation
dt=0.01/phi_HH; %[msec] time step
tspan=0:dt:Time;
t_0=1/phi_HH;  %[ms] pulse width
L_A=300; %length of the A_array (next)
A_array=linspace(0,1,L_A); %the different A we will sample
V_peak=zeros(L_A,L_I); %the voltage peak at each trial
peak_time=zeros(L_A,L_I); %the time of the voltage peak (just to check if this is a valid AP)
mean_v=zeros(L_A,L_I); %mean voltage for T after stimulation
mean_gamma=zeros(L_A,L_I); %mean beta( V, 's') for T after stimulation
mean_delta=zeros(L_A,L_I); %mean alpha( V, 's') for T after stimulation
t_start=20; %start of stimulation pulse index - we wait this time so that all HH variable will relax to steady state
start=t_start/dt; %index of start
T=Time-t_start;  %[ms] Period of stimulation
T_H=20; %maximum timescale of AP depolazrization
AP_indices=start+(1:round(T_H/dt)); %the indices in which an AP can create a large depolarization
mean_gamma_AP= zeros(L_A,L_I); %mean beta( V, 's') for tao_AP after stimulation

%stimulation current
clock=mod(tspan-t_start,T); %stim is given at t=0

% %Initial conditions - Langevin simulation
% y=zeros(length(tspan),9);
% y(1,:)=[-65.0506   0.0526 0.0526 0.0526   0.3172 0.3172 0.3172 0.3172    0.5958]; %initial condition when s is approximatly at steady state

%Initial conditions - Deterministic simulation
y=zeros(length(tspan),4);
y(1,:)=[-64.9   0.05926  0.3193  0.5958]; %initial condition when s is approximatly at steady state


%Main Simulation
tic
for jj=1:L_I
    for kk=1:L_A;
        I=I0.*(clock<t_0);  %Intracellular Pulse
        A=A_array(kk);
        for ii=1:(length(tspan)-1)     
            y(ii+1,:) = y(ii,:) + hhx(y(ii,:),I(ii));
        end
        [V_peak(kk,jj), peak_time(kk,jj)]= max(y(:,1));
        mean_v(kk,jj) = mean(y(start:end,1));
        mean_gamma(kk,jj) = mean(1000*beta(y(start:end,1),'s'));%in [Hz]
        mean_delta(kk,jj) = mean(1000*alpha(y(start:end,1),'s')); %in [Hz]
        mean_gamma_AP(kk,jj) = mean(1000*beta(y(AP_indices,1),'s')); %in [Hz]
    end
end

toc
 %% post processing

threshold_voltage=-10; %if V_peak>threshold_voltage the an AP occured
AP_distribution=V_peak>threshold_voltage;
Threshold_index=[ zeros(1,L_I)>1 ; diff(AP_distribution)>0.5]; %find when v_threshold was crossed upwards - so AP started to occur
Threshold=NaN*AP_distribution;
Threshold(Threshold_index)=0;
shift=40; %calculate rates at this distance next the threshold
gamma_plus=mean_gamma(round((find(Threshold_index)+shift)));
gamma_minus=mean_gamma(round(find(Threshold_index)-shift));
gamma_H=mean_gamma_AP(round((find(Threshold_index)+shift)));
gamma_M=mean_gamma_AP(round(find(Threshold_index)-shift));
gamma_L1=(T*gamma_plus-T_H*gamma_H)/(T-T_H);
gamma_L2=(T*gamma_minus-T_H*gamma_M)/(T-T_H);
approx_gamma=gamma_minus+(gamma_plus-gamma_minus)*AP_distribution;
approx_delta=ones(1,L_A).*mean(mean_delta);

theta=A_array(Threshold_index);
delta=mean_delta(Threshold_index);
% save('parameters for gamma threshold=-17mv',...
%     'T_H','theta','delta','gamma_H','gamma_M','gamma_L1','gamma_L2');



save('average rates.mat');

