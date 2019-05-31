% This file combines results of several "get_Prob_HHSIP.m" simulations
% and saves result (p_{AP}(s)) to a file in both analytical and numerical form.
% This file is then used in various calculations

close all;
clear all;
clc

counter=1; %running average counter
L=4; %number of files
N=1e6;
fontsize=15;


load(['HHSIP_Threshold&Average_rates_I_8.3_N_1e' num2str(log10(N)) '_dt_5e-4_#1.mat']);

Latency_dist_ALL=Latency_dist;
AP_dist_ALL=AP_dist;


for ii=2:L
load(['HHSIP_Threshold&Average_rates_I_8.3_N_1e' num2str(log10(N)) '_dt_5e-4_#' num2str(ii) '.mat']);
counter=ii-1;

Latency_dist_ALL=(1-(1/(1/counter)))*Latency_dist_ALL+(1/counter)*Latency_dist;
AP_dist_ALL=(1-(1/counter))*AP_dist_ALL+(1/counter)*AP_dist;

end


% Average Together all Simulations
Latency_dist=Latency_dist_ALL;
AP_dist=AP_dist_ALL;

s1_array=s1_array(2:(end-1));
s2_array=s2_array(2:(end-1));

[s1_mat s2_mat]=ndgrid(s1_array(pp,:),s2_array(pp,:));

AP_dist=AP_dist(2:(end-1),2:(end-1));

mesh(s1_mat,s2_mat,AP_dist);
xlabel('$s_1$','interpreter','latex','Fontsize',fontsize);
ylabel('$s_2$','interpreter','latex','Fontsize',fontsize);


%% Create fit 
%Define a function which returns the residual between your matrix and your fitted curve
myfun = @(params) 0.5*(1+erf((params(1)*s1_mat+params(2)*s2_mat-params(3))/sqrt(2))) - AP_dist;
%Define initial guesses for parameters a, b
params0 = [1,-1,1];
%Add lots of debugging info
opts = optimset('Display','Iter','TolFun',1e-20);
%Fit
fitparams = lsqnonlin(myfun,params0,[],[],opts)

% Fit this model using new data
p_s=zeros(1,3);
p_s(1)=fitparams(1);
p_s(2)=fitparams(2);
p_s(3)=fitparams(3);     

%plot
mesh(s1_mat,s2_mat,myfun(p_s)+AP_dist);
xlabel('$s_1$','interpreter','latex','Fontsize',fontsize);
ylabel('$s_2$','interpreter','latex','Fontsize',fontsize);

%plot
figure
mesh(s1_mat,s2_mat,myfun(p_s));
xlabel('$s_1$','interpreter','latex','Fontsize',fontsize);
ylabel('$s_2$','interpreter','latex','Fontsize',fontsize);


% save(['HHSIP_p(S)_N1e' num2str(log10(N)) '.mat'],'p_s','AP_dist','s1_mat','s2_mat','I_array','N','dt');

