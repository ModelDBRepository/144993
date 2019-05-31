% This script maps the response of a SHHSIP neuron to a fixed length current pulse
% and saves the results to file.
% We examine whether an AP was created, and at what latency, for different
% combinations of stimulation current amplitude (I0), sodium and pottassium availability (s1, s2) .
%units are msec,mV,mMHO,microFarad,microAmper

% Changes gm=0.09

clear all;
close all;
clc

global I0 t_0 s1 s2 shift

load SHHS_p(S)_N1e6.mat

N=1e6; %noise level
theta=p_s(:,1);
sigma_p_s=100/sqrt(N);

%% Simulation Parameters and arrays
Rep=50; %number of repeations
I_array=[7.5 7.7 7.9 8.1 8.3]; % microAmper
L_I=length(I_array);
v_threshold=-10; %mV - an AP happened if this voltage is crossed
L_s1=50+2; 
L_s2=40+2;
s1_array=zeros(L_I,L_s1); s2_array=zeros(L_I,L_s2); 
AP=zeros(L_s1,L_s2,L_I,Rep);
peaktime=zeros(L_s1,L_s2,L_I,Rep);
K=1000;  %for %[kHz] -> [Hz] coversion, later
Vmax=-inf;

%% HH parameters
phi_HH=2;
VNa=50; VK=-77;VL=-54; %[mV]
gNa=120;gK=36;gL=0.3;gM=0.05*gK; %[mS]   % gm=0.09
Cm=1/phi_HH; %[microFarad] - not 1 Cm as in HH, since spike time is shorter
dy2dt=zeros(1,4);
alpha=zeros(1,3); %recovery rates
beta=alpha; %inactivation rates

%% Stimulation Parameters
dt=0.01/phi_HH;
shift=40; %[ms] when does pulse start
Time=shift+30; tspan=0:dt:Time;
AfterPulse_indices=ceil((shift:dt:Time)/dt);
t_0=1/phi_HH; %[ms] width of stimulation pulse


%% Sodium Slow inactivation kinetics parameters
phi_s1=1/20; %s1 is slower, to get relaxation as in Gal 2010, and 
amp=0.0034*3;sigma=1/0.3;Vhalf=-17; % makes inactivation faster to make 
sigma2=15; %change s2 - makes inactivation faster 
sqrt_dt_N=sqrt(dt/N);

%% Main Simulation
tic

for rr=1:Rep
for pp=1:L_I
    I0=I_array(pp); %strength of stimulation pulse
    I=I0.*(tspan-shift<t_0).*(tspan-shift>0);
    y=zeros(length(tspan),9);
    y(1,:)=[ -66.4379 0.0446 0.0446 0.0446 0.2959 0.2959 0.2959 0.2959 0.6451]'; %initial conditions - neuron start from rest (the same rest?!)
    
    s1_min=theta(pp)-sigma_p_s; s1_max=theta(pp)+sigma_p_s;
    s2_min=0.7; s2_max=1; 
    s1_array(pp,:)=[0 linspace(s1_min,s1_max,L_s1-2) 1]; %always complete s1 to reach 0 and 1
    s2_array(pp,:)=[0 linspace(s2_min,s2_max,L_s2-2) 1] ;
     
    for ii=1:L_s1
        for jj=1:L_s2
            s1=s1_array(pp,ii);
            s2=s2_array(pp,jj);
            for kk=1:(length(tspan)-1)
                
              V=y(kk,1);
              beta(1)= phi_HH*0.125.*exp(-(V+65)./80); %n
              beta(2)= phi_HH*4.*exp(-(V+65)./18); %m
              beta(3)= phi_HH/(exp(-0.1*(V+35))+1); %h

              alpha(1)=phi_HH*0.01*(V+55)./(1-exp(-0.1*(V+55)));%n
              alpha(2)=phi_HH*0.1*(V+40)./(1-exp(-0.1.*(V+40)));%m
              alpha(3)=phi_HH*0.07.*exp(-(V+65)./20); %h
             
              dy2dt(1)=(gNa.*(y(kk,2).*y(kk,3).*y(kk,4)).*y(kk,9).*s1.*(VNa-V)+(y(kk,5).*y(kk,6).*y(kk,7).*y(kk,8)).*(gK+gM.*s2).*(VK-V)+gL.*(VL-V)+I(kk))./Cm;% dVm/dt
              dy2dt(2)=(1-y(kk,2)).*alpha(2)-y(kk,2).*beta(2);  % dm/dt
              dy2dt(3)=(1-y(kk,3)).*alpha(2)-y(kk,3).*beta(2);  % dm/dt  
              dy2dt(4)=(1-y(kk,4)).*alpha(2)-y(kk,4).*beta(2);  % dm/dt  
              dy2dt(5)=(1-y(kk,5)).*alpha(1)-y(kk,5).*beta(1); % dn/dt
              dy2dt(6)=(1-y(kk,6)).*alpha(1)-y(kk,6).*beta(1); % dn/dt
              dy2dt(7)=(1-y(kk,7)).*alpha(1)-y(kk,7).*beta(1); % dn/dt
              dy2dt(8)=(1-y(kk,8)).*alpha(1)-y(kk,8).*beta(1); % dn/dt
              dy2dt(9)=(1-y(kk,9)).*alpha(3)-y(kk,9).*beta(3); % dh/dt
%             dy2dt(10)=(1-y(kk,10)).*alpha(4)-y(kk,10).*beta(4);  % ds1/dt  - notice phi!
%             dy2dt(11)=(1-y(kk,11)).*beta(5)-y(kk,11).*alpha(5);  % ds2/dt  - notice the switch in alpha and beta!

              noise(2)=sqrt(((1-y(kk,2)).*alpha(2)+y(kk,2).*beta(2)));  % dm/dt 
              noise(3)=sqrt(((1-y(kk,3)).*alpha(2)+y(kk,3).*beta(2)));  % dm/dt 
              noise(4)=sqrt(((1-y(kk,4)).*alpha(2)+y(kk,4).*beta(2)));  % dm/dt 
              noise(5)=sqrt(((1-y(kk,5)).*alpha(1)+y(kk,5).*beta(1))); % dn/dt
              noise(6)=sqrt(((1-y(kk,6)).*alpha(1)+y(kk,6).*beta(1))); % dn/dt  
              noise(7)=sqrt(((1-y(kk,7)).*alpha(1)+y(kk,7).*beta(1))); % dn/dt  
              noise(8)=sqrt(((1-y(kk,8)).*alpha(1)+y(kk,8).*beta(1))); % dn/dt  
              noise(9)=sqrt(((1-y(kk,9)).*alpha(3)+y(kk,9).*beta(3))); % dh/dt

              y(kk+1,:) = y(kk,:)+dy2dt*dt+sqrt_dt_N*randn(1,9).*noise;  
            end                   
            
            [Vmax peaktime(ii,jj,pp,rr)]=max(y(:,1));
            if Vmax>v_threshold
               AP(ii,jj,pp,rr)=1;             
            end

        end
    end
end

end
% 
toc

Latency_dist=sum((peaktime*dt-shift).*AP,4)./sum(AP,4);
AP_dist=mean(AP,4);


% save('temp.mat');
save('HHSIP_Threshold&Average_rates_I_7.5-8.3_N_1e6_dt5e-3_#4.mat');
% save('SHHS_p(S).mat','AP_dist','s1_array','I_array');



%% Figures
% close all;
% clc

% load('2D_Threshold&Average_rates.mat');
fontsize=15;
scrsz = get(0,'ScreenSize'); %get screen size for figures  
figure('Position',[scrsz(3)*0 scrsz(4)*0 scrsz(3) scrsz(4)]);

for pp=1:L_I
    
[s1_mat s2_mat]=ndgrid(s1_array(pp,:),s2_array(pp,:));
hold all;

subplot(2,1,1);
contourf(s1_mat,s2_mat,squeeze(AP_dist(:,:,pp)));
xlabel('$s_1$','interpreter','latex','Fontsize',fontsize);
ylabel('$s_2$','interpreter','latex','Fontsize',fontsize);


subplot(2,1,2);
mesh(s1_mat,s2_mat,squeeze(Latency_dist(:,:,pp)));
xlabel('$s_1$','interpreter','latex','Fontsize',fontsize);
ylabel('$s_2$','interpreter','latex','Fontsize',fontsize);
zlabel('Latency [ms]','interpreter','latex','Fontsize',fontsize);
xlim([s1_min s1_max]);
ylim([s2_min s2_max]);
zlim([0 1.2*max(max(squeeze(Latency_dist(:,:,pp))))+eps]);


end
 
% 
% set(gcf, 'Color', 'w');
% export_fig HHSAP_Threshold_Latency_AVG_rates_I0_8.tif -painters % cool function - check internt for support