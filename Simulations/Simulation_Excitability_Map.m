% A simulation of periodical stimulation according to the
% "excitability map" of the HHS model

close all
clear all
clc


%% Input parameters - N,Time,f_in, M, k, p_s,params,alpha_scaling
N=1e6; %number of channels
Time=1e3*3600*55; %[ms]
f_in=20; %[Hz]
I0=7.9; %[microampere]

load(['SHHS_p(S)_N1e' num2str(log10(N),2) '.mat'],'p_s','I_array');
load('params_7.5_7.7_7.9_8.1_8.3_dt=5e-4.mat','params');

M=1; %number of channel types
k=0; 
alpha_scaling=0;

%% Parameter Initialization
L_AP=round(1e-3*Time*f_in);

k_v= (k.^(0:(M-1)))';
N_v=round(N.*k_v.^alpha_scaling);%number of channel
w_v= ones(M,1)/M; %weight of each channel type

ii=find(I_array==I0);

theta=p_s(ii,1); b=p_s(ii,2);
prob=@(x) 0.5*(1+erf((x-theta)/(sqrt(2)*b))); % p_AP(s)

delta=params(ii,3).*k_v;  % [Hz]
% delta_p=delta;  % [Hz] during AP 
% delta_m=delta;  % [Hz] during stimulation
% delta_0=delta;  % [Hz] during rest
tau_AP=1e-3*params(ii,1);  % [Hz]
gamma_p=params(ii,4).*k_v;  % [Hz] during AP 
gamma_m=params(ii,5).*k_v;  % [Hz] during stimulation
gamma_0=params(ii,6).*k_v;  % [Hz] during rest   

%% Simulation parameters - for uncoupled two-state channels only

T_TIME=1/f_in;  %[Sec] - cycle time of Stimulations
H_TIME=tau_AP; %[Sec] - tau_AP
L_TIME=T_TIME-H_TIME; %[Sec] - time interval between Stimulations spikes
                                   
rate_h=gamma_p+delta;
P_inf_h=delta./(gamma_p+delta);
rate_m=gamma_m+delta;
P_inf_m=delta./(gamma_m+delta);
rate_l=gamma_0+delta;
P_inf_l=delta./(gamma_0+delta);

P_H=P_inf_h+(1-P_inf_h).*exp(-H_TIME*rate_h);          %The probability to remain available after H_TIME, with AP
P_I_AP=1-(P_inf_l+(P_H-P_inf_l).*exp(-L_TIME*rate_l));  %The probability to inactivate after T_TIME, with AP

P_H=P_inf_m+(1-P_inf_m).*exp(-H_TIME*rate_m);   %The probability to remain available after H_TIME, with no AP
P_I_NO_AP=1-(P_inf_l+(P_H-P_inf_l).*exp(-L_TIME*rate_l));  %The probability to inactivate after T_TIME, with no AP

P_H=P_inf_h+(0-P_inf_h).*exp(-H_TIME*rate_h);          %The probability to become available after H_TIME, with AP
P_A_AP=P_inf_l+(P_H-P_inf_l).*exp(-L_TIME*rate_l);  %The probability to become available after T_TIME, with AP

P_H=P_inf_m+(0-P_inf_m).*exp(-H_TIME*rate_m);          %The probability to become available after H_TIME, with no AP
P_A_NO_AP=P_inf_l+(P_H-P_inf_l).*exp(-L_TIME*rate_l);  %The probability to become available after T_TIME, with no AP


%% Array Initialization                                   

AP=zeros(L_AP,1);              % A binary vector marking at which stimulus was an AP genetated
AP(1)=1;         
s_num=zeros(L_AP,M);                 % Number of channels in each state
s_num(1,:)=round(N_v'*theta);
s=s_num;              % Number of channels in each state
s(1,:)=s_num(1,:)./N_v';

%% Main Loop
tic

for tt=2:L_AP
    
if AP(tt-1)  % If there was an Action Potential
    P_I=P_I_AP;
    P_A=P_A_AP;
else % If there was no Action Potential
   	P_I=P_I_NO_AP;
    P_A=P_A_NO_AP;
end

lambda_A=P_A'.*(N_v'-s_num(tt-1,:));
lambda_I=P_I'.*s_num(tt-1,:);

seed_A = random('Poisson',lambda_A,1,M);
seed_I = random('Poisson',lambda_I,1,M);

s_num(tt,:)=s_num(tt-1,:)+seed_A-seed_I;
s(tt,:)=s_num(tt,:)./N_v'; %notice that individual channels are weighted inversely to their channel number type

if prob(s(tt,:)*w_v)>rand
    AP(tt)=1;
else
    AP(tt)=0;
end

end

toc

T=T_TIME*1e3;


save('HHS_reduced_sim_I0_7.9_N_1e6.mat','f_in','Time','I0','alpha_scaling','k','M','N','AP','s_num','s','w_v');