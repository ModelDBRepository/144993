% This script maps the response probability of a Coupled HHS neuron to
% a fixed length current pulse, as a function of of stimulation current
% amplitude (I0) and sodium availability. The result is saved the results to a file.
% To improve quality of results, run several iterations of this file, and
% then combine results using "Combine_Get_Prob_s_Reps.m" to generate p(s)

% units are msec,mV,mMHO,microFarad,microAmper

clear all;
close all;
clc

global I0 t_0 s1 s2 shift

%% Simulation Parameters and arrays
Rep=50; %number of repeations
I_array=8.3%[7.5 7.7 7.9 8.1 8.3]; % microAmper
L_I=length(I_array);
v_threshold=-10; %mV - an AP happened if this voltage is crossed
s1_min=0.8; s1_max=1;
s2_min=0; s2_max=0; 
L_s1=300; s1_array=[0 linspace(s1_min,s1_max,L_s1-1)]; %always complete s1 to reach 0 and 1
L_s2=1; s2_array=linspace(s2_min,s2_max,L_s2);
[s1_mat s2_mat]=ndgrid(s1_array,s2_array);
Mean_Voltage=zeros(L_s1,L_s2,L_I,Rep); %the average voltage during simulation (AP+rest)
Mean_gamma1=zeros(L_s1,L_s2,L_I,Rep); %the average gamma during simulation (AP+rest)
Mean_gamma2=zeros(L_s1,L_s2,L_I,Rep); %the average gamma during simulation (AP+rest)
Mean_delta1=zeros(L_s1,L_s2,L_I,Rep); %the average gamma during simulation (AP+rest)
Mean_delta2=zeros(L_s1,L_s2,L_I,Rep); %the average gamma during simulation (AP+rest)
AP=zeros(L_s1,L_s2,L_I,Rep);
peaktime=zeros(L_s1,L_s2,L_I,Rep);
Vmax=zeros(L_s1,L_s2,L_I,Rep);
K=1000;  %for %[kHz] -> [Hz] coversion, later


%% HH parametersy
phi_HH=2;
VNa=50; VK=-77;VL=-54; %[mV]
gNa=120;gK=36;gL=0.3;gM=0.01*gK; %[mS]
Cm=1/phi_HH; %[microFarad] - not 1 Cm as in HH, since spike time is shorter
dy2dt=zeros(1,4);
alpha=zeros(1,3); %recovery rates
beta=alpha; %inactivation rates
N=1e6; %noise level

%% Stimulation Parameters
dt=0.01/phi_HH;
shift=50; %[ms] when does pulse start
Time=shift+20; tspan=0:dt:Time;
AfterPulse_indices=ceil((shift:dt:Time)/dt);
t_0=1/phi_HH; %[ms] width of stimulation pulse

%% Sodium Slow inactivation kinetics parameters
phi_s1=1/20; %s1 is slower, to get relaxation as in Gal 2010, and 
amp=0.0034*3;sigma=1/0.3;Vhalf=-17; % makes inactivation faster to make 
sigma2=15; %change s2 - makes inactivation faster 
sqrt_dt_N=sqrt(dt/N);
eye4=eye(4);

%% Main Simulation
tic

for rr=1:Rep
for pp=1:L_I
    I0=I_array(pp); %strength of stimulation pulse
    I=I0.*(tspan-shift<t_0).*(tspan-shift>0);
    y=tspan*0;
    V= -66.4379; %initial conditions - neuron start from rest (the same rest?!)
    bn= phi_HH*0.125.*exp(-(V+65)./80); %n
    bm= phi_HH*4.*exp(-(V+65)./18); %m
    bh= phi_HH*1/(exp(-0.1*(V+35))+1); %h

    am=phi_HH*0.01*(V+55)./(1-exp(-0.1*(V+55)));%n
    an=phi_HH*0.1*(V+40)./(1-exp(-0.1.*(V+40)));%m
    ah=phi_HH*0.07.*exp(-(V+65)./20); %h

    n=[4*bn^3*an 6*an^2*bn^2 4*an^3*bn an^4]'/(an+bn)^4;
    temp=[am^3 3*bm*am^2  3*bm^2*am bm^3]'/(am+bm)^3/(ah+bh); 
    m=[temp*ah ; temp*bh];  
    
    for ii=1:L_s1
        for jj=1:L_s2
            s1=s1_array(ii);
            s2=s2_array(jj);
            
            for kk=1:length(tspan)
                  y(kk)=V;
      % rates
                  bn= phi_HH*0.125.*exp(-(V+65)./80); %n
                  bm= phi_HH*4.*exp(-(V+65)./18); %m
                  bh= phi_HH*1/(exp(-0.1*(V+35))+1); %h
 
                  an=phi_HH*0.01*(V+55)./(1-exp(-0.1*(V+55)));%n
                  am=phi_HH*0.1*(V+40)./(1-exp(-0.1.*(V+40)));%m
                  ah=phi_HH*0.07.*exp(-(V+65)./20); %h
                
                  % Voltage step
                  V=V+dt*(s1*gNa.*m(1).*(VNa-V)+n(4).*gK.*(VK-V)+gL.*(VL-V)+I(kk))./Cm;% dVm/dt
                
                % Pottasium channel step - here we used normalization
                  n0=1-sum(n); % unrepresented state - found by normalization
                  AK= [ -3*an-bn 2*bn 0 0 ;
                        3*an -2*an-2*bn 3*bn 0 ;
                        0 2*an -an-3*bn 4*bn;
                        0 0 an -4*bn];
                  aK=4*an*[ 1; 0; 0; 0];
                  
                  SK= [ -sqrt(4*an*n0+bn*n(1)) sqrt(3*an*n(1)+2*bn*n(2)) 0 0 ;
                        0 -sqrt(3*an*n(1)+2*bn*n(2)) sqrt(2*an*n(2)+3*bn*n(3)) 0;
                        0 0 -sqrt(2*an*n(2)+3*bn*n(3)) sqrt(an*n(3)+4*bn*n(4)) ;
                        0 0 0 -sqrt(an*n(3)+4*bn*n(4))];
                  
                  n=n+dt*(AK*n+aK*n0)+sqrt_dt_N*SK*randn(4,1);
                  
                % Sodium channel step - here we did not use normalization
                
                  ANa_m= [ -3*bm am 0 0;
                            3*bm -(am+2*bm) 2*am 0 ;
                            0 2*bm -(2*am+bm) 3*am ;
                            0 0 bm -3*am];
                  ANa=[ (ANa_m-eye4*bh) ah*eye4 ;
                         bh*eye4 (ANa_m-eye4*ah);];
                     
                                 
                  SNa=zeros(8,10); %(number of states) X (number of transition pairs)
                  index=1;
                  for ll=1:8
                      for uu=ll+1:8
                        if ANa(ll,uu)>0
                            temp=sqrt(abs(ANa(ll,uu)*m(uu)+ANa(uu,ll)*m(ll))); %notice abs - must have this here since probabilities get very small
                            SNa(ll,index)=temp;
                            SNa(uu,index)=-temp;
                            index=index+1;
                        end
                      end
                  end
                 
                  m=m+dt*ANa*m+sqrt_dt_N*SNa*randn(10,1);
                
            end
            
            temp_V=y(AfterPulse_indices);
            gamma1_temp= phi_s1*amp./(exp(-(temp_V+Vhalf)/sigma)+1); % s2 - orignally from Fleidervish1996 
            s2_inf_temp=1./(1+exp(-(temp_V+35)/10));
            tao_s2_temp=1e3./(3.3*exp((temp_V+35)/sigma2)+exp(-(temp_V+35)/20));
            gamma2_temp=(1-s2_inf_temp)./tao_s2_temp; %s2- changed from Yamada1989  !!!!
            delta1_temp=phi_s1*0.001.*exp(-(temp_V+85)./30); %s1
            delta2_temp = s2_inf_temp./tao_s2_temp; %0.0034./(exp(-0.1*(V+17))+1);- changed from Yamada1989    
            
            
            [Vmax(ii,jj,pp,rr) peaktime(ii,jj,pp,rr)]=max(y);
            if Vmax(ii,jj,pp,rr)>v_threshold
               AP(ii,jj,pp,rr)=1;             
            end
            Mean_Voltage(ii,jj,pp,rr)=mean(temp_V);
            Mean_gamma1(ii,jj,pp,rr)=mean(gamma1_temp);
            Mean_gamma2(ii,jj,pp,rr)=mean(gamma2_temp);
            Mean_delta1(ii,jj,pp,rr)=mean(delta1_temp); 
            Mean_delta2(ii,jj,pp,rr)=mean(delta2_temp);
        end
    end
end

end
% 
toc

Latency_dist=sum((peaktime*dt-shift).*AP,4)./sum(AP,4);
AP_dist=mean(AP,4);
Mean_gamma1_dist=K*mean(Mean_gamma1,4); %[kHz] -> [Hz]
Mean_gamma2_dist=K*mean(Mean_gamma2,4); %[kHz] -> [Hz]
Mean_delta1_dist=K*mean(Mean_delta1,4); %[kHz] -> [Hz]
Mean_delta2_dist=K*mean(Mean_delta2,4); %[kHz] -> [Hz]

% save('temp.mat');
save('SHHS_Threshold&Average_rates_I_8.3_N_1e6_dt_5e-3_#1.mat');
% save('SHHS_p(S).mat','AP_dist','s1_array','I_array');



%% Figures
% close all;
% clc

% load('2D_Threshold&Average_rates.mat');
fontsize=15;
scrsz = get(0,'ScreenSize'); %get screen size for figures  
figure('Position',[scrsz(3)*0 scrsz(4)*0 scrsz(3) scrsz(4)]);

for pp=1:L_I
hold all;

subplot(2,2,1);
plot(s1_array,squeeze(AP_dist(:,:,pp)));
xlabel('$s_1$','interpreter','latex','Fontsize',fontsize);
ylabel('$p(s)$','interpreter','latex','Fontsize',fontsize);
xlim([s1_min s1_max]);
ylim([0 1]);
hold all;

subplot(2,2,3);
plot(s1_array,squeeze(Latency_dist(:,:,pp)));
xlabel('$s_1$','interpreter','latex','Fontsize',fontsize);
ylabel('Latency [ms]','interpreter','latex','Fontsize',fontsize);
xlim([s1_min s1_max]);
hold all;


% subplot(2,3,1);
% mesh(s1_mat,s2_mat,squeeze(Mean_Voltage(:,:,pp)));
% xlabel('$s_1$','interpreter','latex','Fontsize',fontsize);
% ylabel('$s_2$','interpreter','latex','Fontsize',fontsize);
% zlabel('$\bar V$ $(s_1,s_2)$ [mV]','interpreter','latex','Fontsize',fontsize);

subplot(2,2,2);
plot(s1_array,squeeze(Mean_gamma1_dist(:,:,pp)));
xlabel('$s_1$','interpreter','latex','Fontsize',fontsize);
ylabel('$\bar \gamma_1$ $(s_1,s_2)$ [Hz]','interpreter','latex','Fontsize',fontsize);
xlim([s1_min s1_max]);
ylim([0 1.2*max(max(squeeze(Mean_gamma1_dist(:,:,pp))))]);
hold all;


subplot(2,2,4);
plot(s1_array,squeeze(Mean_delta1_dist(:,:,pp)));
xlabel('$s_1$','interpreter','latex','Fontsize',fontsize);
ylabel('$\bar \delta_1$ $(s_1,s_2)$ [Hz]','interpreter','latex','Fontsize',fontsize);
xlim([s1_min s1_max]);
ylim([0 1.2*max(max(squeeze(Mean_delta1_dist(:,:,pp))))]);
hold all;


end
 
% 
% set(gcf, 'Color', 'w');
% export_fig HHSAP_Threshold_Latency_AVG_rates_I0_8.tif -painters % cool function - check internt for support