% This file combines results of several "get_Prob_s.m" simulations
% and saves result (p_{AP}(s)) to a file in both analytical and numerical form.
% This file is then used in various calculations

close all;
clear all;
clc

counter=1; %running average counter
L=4; %number of files
N=1e6;

load(['SHHS_Threshold&Average_rates_I_8.3_N_1e' num2str(log10(N)) '_dt_5e-3_#1.mat']);

Latency_dist_ALL=Latency_dist;
AP_dist_ALL=AP_dist;
Mean_gamma1_dist_ALL=Mean_gamma1_dist; 
Mean_gamma2_dist_ALL=Mean_gamma2_dist; 
Mean_delta1_dist_ALL=Mean_delta1_dist; 
Mean_delta2_dist_ALL=Mean_delta2_dist; 

for ii=2:L
load(['SHHS_Threshold&Average_rates_I_8.3_N_1e' num2str(log10(N)) '_dt_5e-3_#' num2str(ii) '.mat']);
counter=ii-1;

Latency_dist_ALL=(1-(1/(1/counter)))*Latency_dist_ALL+(1/counter)*Latency_dist;
AP_dist_ALL=(1-(1/counter))*AP_dist_ALL+(1/counter)*AP_dist;
Mean_gamma1_dist_ALL=(1-(1/counter))*Mean_gamma1_dist_ALL+(1/counter)*Mean_gamma1_dist; 
Mean_gamma2_dist_ALL=(1-(1/counter))*Mean_gamma2_dist_ALL+(1/counter)*Mean_gamma2_dist; 
Mean_delta1_dist_ALL=(1-(1/counter))*Mean_delta1_dist_ALL+(1/counter)*Mean_delta1_dist; 
Mean_delta2_dist_ALL=(1-(1/counter))*Mean_delta2_dist_ALL+(1/counter)*Mean_delta2_dist;  

end



%% Average Together all Simulations
Latency_dist=Latency_dist_ALL;
AP_dist=AP_dist_ALL;
Mean_gamma1_dist=Mean_gamma1_dist_ALL;
Mean_gamma2_dist=Mean_gamma2_dist_ALL;
Mean_delta1_dist=Mean_delta1_dist_ALL;
Mean_delta2_dist=Mean_delta2_dist_ALL;

%% Create fit 
fo_ = fitoptions('method','NonlinearLeastSquares','Lower',[-Inf    0]);
st_ = [0.617524930342367 0.80457760195393657 ];
set(fo_,'Startpoint',st_);
ft_ = fittype('0.5*(1+erf((x-a)/(sqrt(2)*b)))',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a', 'b'});

% Fit this model using new data
p_s=zeros(L_I,2);
for ii=1:L_I
     temp= fit(s1_array',squeeze(AP_dist(:,:,ii)),ft_,fo_);
     p_s(ii,1)=temp.a;
     p_s(ii,2)=temp.b;
end
save(['Coupled_HHS_I_8.3_p(S)_N1e' num2str(log10(N)) '_dt5e-3.mat'],'p_s','AP_dist','s1_array','I_array','N','dt');

