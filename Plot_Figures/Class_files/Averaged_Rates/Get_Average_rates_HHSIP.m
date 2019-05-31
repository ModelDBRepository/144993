% Finds parameters of HHSIP model and saves them  to file -for different currents

clear all;
close all;
clc

% gets [T_H,thetaS,delta,gamma_H,gamma_M,gamma_L,gamma_plus,gamma_minus] for all I0

I0_array=7.5:0.2:8.3; %[microamper]/cm^2 %Range(Fleidervish1996):  7.3-9.4
L_I=length(I0_array);

cell_params=cell(L_I,1); 

tic
for jj=1:L_I
    I0=I0_array(jj);    
    cell_params(jj,:)={get_params(I0)};
end
toc

save('HHSIP_params_7.5_7.7_7.9_8.1_8.3.mat')