function  res = get_params(I0)
% %this function extracts
% res=[T_H,w1,delta1,gamma1_H,gamma1_M,gamma1_L,gamma1_plus,gamma1_minus ;
%      theta,w2,delta2,gamma2_H,gamma2_M,gamma2_L,gamma2_plus,gamma2_minus];
% notice the range of s1_array,s2_array  and adjust to include threshold

phi_HH=2;
VNa=50; VK=-77;VL=-54; %[mV]
gNa=120;gK=36;gL=0.3;gM=0.01*gK; %[mS]
Cm=1/phi_HH; %[microFarad/cm^2]
phi_s1=1/20; %s1 is slower, to get relaxation as in Gal 2010, and 
amp=0.0034*3;sigma=1/0.3;Vhalf=-17;% makes inactivation faster to make 
sigma2=15;  % increase steepnes of inactivation cureve in s2 

Time=70; %[msec] length of simulation
dt=0.001/phi_HH; %[msec] time step
tspan=0:dt:Time;
t_0=1/phi_HH;  %[ms] pulse width
t_start=20; %start of stimulation pulse index - we wait this time so that all HH variable will relax to steady state
start=t_start/dt; %index of start
T_H=15; %maximum timescale of AP depolazrization
AP_indices=start+(1:round(T_H/dt)); %the indices in which an AP can create a large depolarization
last_AP_indices=max(AP_indices); %the last of the AP indices

L_s1=1e3; %length of the A_array (next)
L_s2=1e2;
s1_array=linspace(0.85,0.95,L_s1); %the different A we will sample
s2_array=linspace(0,1,L_s2); %the different A we will sample

V_peak=zeros(L_s1,L_s2); %the voltage peak at each trial
mean_gamma1=zeros(L_s1,L_s2); %mean beta( V, 's1') for T after stimulation
mean_delta1=zeros(L_s1,L_s2); %mean alpha( V, 's1') for T after stimulation
mean_gamma1_AP= zeros(L_s1,L_s2); %mean beta( V, 's1') for H_H after stimulation
mean_gamma1_NO_AP=zeros(L_s1,L_s2); %mean beta( V, 's1') for T>t>T_H 

mean_delta2=zeros(L_s1,L_s2); %mean beta( V, 's2') for T after stimulation
mean_gamma2=zeros(L_s1,L_s2); %mean alpha( V, 's2') for T after stimulation
mean_gamma2_AP= zeros(L_s1,L_s2); %mean beta( V, 's2') for H_H after stimulation
mean_gamma2_NO_AP=zeros(L_s1,L_s2); %mean beta( V, 's2') for T>t>T_H 
scale=1e3; %change scale from kHZ to Hz

%stimulation current
I=I0.*(tspan<t_start+t_0).*(tspan>t_start);  %Intracellular Pulse

%Initial conditions - Langevin simulation
y=zeros(length(tspan),4);
y(1,:)=[-65.0506 0.0526 0.3172 0.5958]; %initial condition when s is approximatly at steady state

%Main Simulation

for kk=1:L_s1;
    for jj=1:L_s2    
        s1=s1_array(kk);
        s2=s2_array(jj);
        
    for ii=1:(length(tspan)-1)
        
    dy2dt=zeros(1,4);beta=zeros(3,1);alpha=beta;
    V=y(ii,1); m=y(ii,2); n=y(ii,3); h=y(ii,4);

      beta(1)= phi_HH*0.125.*exp(-(V+65)./80); %n
      beta(2)= phi_HH*4.*exp(-(V+65)./18); %m
      beta(3)= phi_HH*1/(exp(-0.1*(V+35))+1); %h

      alpha(1)= phi_HH*0.01*(V+55)./(1-exp(-0.1*(V+55)));%n
      alpha(2)= phi_HH*0.1*(V+40)./(1-exp(-0.1.*(V+40)));%m
      alpha(3)= phi_HH*0.07.*exp(-(V+65)./20); %h

      dy2dt(1)=(gNa.*(m.^3).*h.*s1.*(VNa-V)+(gK+gM.*s2).*(n.^4).*(VK-V)+gL.*(VL-V)+I(ii))./Cm;% dVm/dt
      dy2dt(2)=(1-m).*alpha(2)-m.*beta(2);  % dm/dt
      dy2dt(3)=(1-n).*alpha(1)-n.*beta(1); % dn/dt
      dy2dt(4)=(1-h).*alpha(3)-h.*beta(3); % dh/dt
        
    y(ii+1,:) = y(ii,:) + dt*dy2dt;
    end
    
    temp_V=y(:,1);
    gamma1_temp= phi_s1*amp./(exp(-(temp_V-Vhalf)/sigma)+1); % s1 - orignally from Fleidervish1996 
    s2_inf_temp=1./(1+exp(-(temp_V+35)/10));
    tao_s2_temp=1e3./(3.3*exp((temp_V+35)/sigma2)+exp(-(temp_V+35)/20));
    delta2_temp=(1-s2_inf_temp)./tao_s2_temp; %s2- switched from Yamada1989 !!
    delta1_temp=phi_s1*0.001.*exp(-(temp_V+85)./30); %s1
    gamma2_temp = s2_inf_temp./tao_s2_temp; %s2- switched from Yamada1989 !!

    V_peak(kk,jj) = max(y(:,1));
    mean_gamma1(kk,jj) = scale*mean(gamma1_temp(start:end));%in [Hz]
    mean_delta1(kk,jj)  = scale*mean(delta1_temp(start:end)); %in [Hz]
    mean_gamma1_AP(kk,jj)  = scale*mean(gamma1_temp(AP_indices)); %in [Hz]
    mean_gamma1_NO_AP(kk,jj) =  scale*mean(gamma1_temp(last_AP_indices+1:end)); %in [Hz]
    mean_delta2(kk,jj)  = scale*mean(delta2_temp(start:end));%in [Hz]
    mean_gamma2(kk,jj)  = scale*mean(gamma2_temp(start:end)); %in [Hz]
    mean_gamma2_AP(kk,jj)  = scale*mean(gamma2_temp(AP_indices,1)); %in [Hz]
    mean_gamma2_NO_AP(kk,jj) =  scale*mean(gamma2_temp(last_AP_indices+1:end)); %in [Hz]
    end
end

%% Calculate rates
threshold_voltage=-10; %if V_peak>threshold_voltage the an AP occured
AP_distribution=V_peak>threshold_voltage;
Threshold_index=[ zeros(1,L_s2) ; diff(AP_distribution)] >0.5; %find when v_threshold was crossed upwards - so AP started to occur
[ind1_mat ind2_mat]=ndgrid(1:L_s1,1:L_s2);
mid_threshold_index1=round(mean(ind1_mat(Threshold_index)));
mid_threshold_index2=round(mean(ind2_mat(Threshold_index)));
shift=1; %calculate rates at this distance next the threshold

gamma1_plus=mean_gamma1(mid_threshold_index1+shift,mid_threshold_index2);
gamma1_minus=mean_gamma1(mid_threshold_index1-shift,mid_threshold_index2);
gamma1_H=mean_gamma1_AP(mid_threshold_index1+shift,mid_threshold_index2);
gamma1_M=mean_gamma1_AP(mid_threshold_index1-shift,mid_threshold_index2);
gamma1_L=mean_gamma1_NO_AP(mid_threshold_index1,mid_threshold_index2);
delta1=mean_delta1(mid_threshold_index1,mid_threshold_index2);

gamma2_plus=mean_gamma2(mid_threshold_index1+shift,mid_threshold_index2);
gamma2_minus=mean_gamma2(mid_threshold_index1-shift,mid_threshold_index2);
gamma2_H=mean_gamma2_AP(mid_threshold_index1+shift,mid_threshold_index2);
gamma2_M=mean_gamma2_AP(mid_threshold_index1-shift,mid_threshold_index2);
gamma2_L=mean_gamma2_NO_AP(mid_threshold_index1,mid_threshold_index2);
delta2=mean_delta2(mid_threshold_index1,mid_threshold_index2);

w1=1;
theta=w1*s1_array(Threshold_index(:,1)'); %assumin s2_min=0
w2=theta-w1*s1_array(Threshold_index(:,end)'); %this calculation is only true if threshold  touches s2=1;

res=[T_H,w1,delta1,gamma1_H,gamma1_M,gamma1_L,gamma1_plus,gamma1_minus ;
     theta,w2,delta2,gamma2_H,gamma2_M,gamma2_L,gamma2_plus,gamma2_minus];

end
