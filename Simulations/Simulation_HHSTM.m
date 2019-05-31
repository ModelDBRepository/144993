% Generates data file of intermittent mode in the the HHSTM (uncoupled stochastic HH
% model (N=1e6) with synaptic depression a-la Tsodyks & Markram (1996) ) 

clear all;
close all;
clc

% This simulation  does not require any sub-functions
% 1D slow inactivaion (gM=0)+stochasticitys1

% HH rates and Capacitance adjusted, Also s rates (gamma and delta)

intermittent=1; %if intermittent mode is required, set as 1, if transient mode required, set as 0
I_array=70;% [microamper]  % 23 [microamper] is 
f_array=20;%[1 5 10 12.5 15 17.5 20 22.5 25 27.5 30 35 40 45]; % [Hz]
L_I=length(I_array);
T_array=1e3./f_array;   %[ms]
L_f=length(f_array);
Time=1e3*55*3600; %[msec] length of simulation

cell_AP=cell(L_I,L_f);cell_R=cell_AP;cell_s1=cell_AP;
cell_Latency=cell_AP;cell_Amplitude=cell_AP; %sampling arrays
cell_stim_time=cell_AP;

%% HH parameters
N=1e6; %number of conducting channels - determines noise level
N_s=1e6; %number of (?) - determines noise level in synapse
phi_HH=2;%2; % make model faster
VNa=50; VK=-77;VL=-54; %[mV]
gNa=120;gK=36;gL=0.3;gM=0;%0.01; %[mS]   
Cm=1/phi_HH; %[microFarad] - not 1 Cm as in HH, since spike time is shorter
dt=0.01/phi_HH; %[msec] time step
t_0=1/phi_HH;  %[ms] pulse width

HH_flag=1; % have HH model m,n,h kinetics

%% Slow kinetics parameters
phi_s=1/20; %s is slower, to get relaxation as in Gal 2010, and 
amp=0.0034*3;sigma=1/0.3;Vhalf=-17; % makes inactivation faster to make
sigma2=15;
sqrt_dt_N=sqrt(dt/N);

% synapse model
tau_inact=3; %[msec]
tau_rec=800;%[msec]
U_se=0.67;
sqrt_dt_N_s=sqrt(dt/N_s);

% Rin=100e6; %[ohm]
% A_se=250; %[pA]
% tau_mem-50; %[msec]

%% Main Simulation
tic
for jj=1:L_I
    I0=I_array(jj); % [microamper] 
    for ii=1:L_f

        f_in=f_array(ii);
        T=1e3/f_in; %[msec]
        stim_num=round(Time/T);
        stim_time=T*(0:stim_num-1);
        AP=zeros(stim_num,1);R_array=AP; s1=AP;Latency=NaN*AP;Amplitude=NaN*AP; %sampling arrays

        RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock))); % so rand won't generate the same pattern each time...
        defaultStream = RandStream.getDefaultStream;

        %set initial conditions  and intialize arrays
        if intermittent==1
            s_initial=0.8735; %start just above threshold for intermittent mode
            R0=0.06;
            E0=0;
        else 
            s_initial=1;
            R0=1;
            E0=0;
        end
        I0=I_array(jj);
        y0=[  -64.9603    0.0532    0.0534    0.0530    0.3182    0.3182    0.3180    0.3183    0.5954 s_initial 0.77];%initial condition when s is approximatly at steady state
        y=y0; %intial conditions 
        R=R0;
        E=E0;

        dy2dt=zeros(1,11);
        noise=zeros(1,11); %conductance noise
        noise_s=zeros(1,3); %synaptic noise
        alpha=zeros(1,5); %recovery rates
        beta=alpha; %inactivation rates

        V_threshold=-10; %[mV] AP is defined if V crosses V_threshold upwards
        stim_flag=0;
        delay=0;
        Vmax=-inf;

        cycle=round(T/dt); %note that this is an approximation of the resulted T, not the real one
        
        % continuous sample - remove at long simulations
%         v_array=zeros(stim_num*cycle,1);
%         r_array=zeros(stim_num*cycle,1);
%         E_array=zeros(stim_num*cycle,1);
%         s_array=zeros(stim_num*cycle,1);

        for pp=1:stim_num
            R_array(pp)=R;
            s1(pp)=y(10);
            delay=-dt;
            stim_flag=1;

            for kk=1:cycle
                if stim_flag==1  %sample
                    delay=delay+dt;
                    if (y(1)>V_threshold)
                        AP(pp)=1;
                        Vmax=y(1);
                        stim_flag=0;
                    end
                end

                if Vmax>-inf
                    if (y(1)<Vmax) %sample
                        Latency(pp)=delay;
                        Amplitude(pp)=Vmax;
                        Vmax=-inf;
                    else
                        Vmax=y(1);
                        delay=delay+dt;
                    end                
                end          


                %% SO_hhx_Langevin inserted here     
                %fast kinetics
                  V=y(1);
                  beta(1)= HH_flag*phi_HH*0.125.*exp(-(V+65)./80); %n
                  beta(2)= HH_flag*phi_HH*4.*exp(-(V+65)./18); %m
                  beta(3)= HH_flag*phi_HH*1/(exp(-0.1*(V+35))+1); %h
                  
                  alpha(1)=HH_flag*phi_HH*0.01*(V+55)./(1-exp(-0.1*(V+55)));%n
                  alpha(2)=HH_flag*phi_HH*0.1*(V+40)./(1-exp(-0.1.*(V+40)));%m
                  alpha(3)=HH_flag*phi_HH*0.07.*exp(-(V+65)./20); %h
                  
                  %slow kinetics
                  beta(4)= phi_s*amp./(exp(-(V-Vhalf)/sigma)+1); % s1 - different then Fleidervish1996 !!
                  alpha(4)=phi_s*0.001.*exp(-(V+85)./30); %s1
%                  s2_inf=1./(1+exp(-(V+35)/10));
%                  tao_s2=1e3./(3.3*exp((V+35)/sigma2)+exp(-(V+35)/20));
%                  beta(5)=(1-s2_inf)./tao_s2; %s2- orignally from Yamada1989

%                  alpha(5) = s2_inf./tao_s2; %0.0034./(exp(-0.1*(V+17))+1);- orignally from Yamada1989        

                dy2dt(1)=(gNa.*(y(2).*y(3).*y(4)).*y(9).*y(10).*(VNa-y(1))+(y(5).*y(6).*y(7).*y(8)).*(gK+gM.*y(11)).*(VK-y(1))+gL.*(VL-y(1))+I0*E)./Cm;% dVm/dt
                dy2dt(2)=(1-y(2)).*alpha(2)-y(2).*beta(2);  % dm/dt
                dy2dt(3)=(1-y(3)).*alpha(2)-y(3).*beta(2);  % dm/dt  
                dy2dt(4)=(1-y(4)).*alpha(2)-y(4).*beta(2);  % dm/dt  
                dy2dt(5)=(1-y(5)).*alpha(1)-y(5).*beta(1); % dn/dt
                dy2dt(6)=(1-y(6)).*alpha(1)-y(6).*beta(1); % dn/dt
                dy2dt(7)=(1-y(7)).*alpha(1)-y(7).*beta(1); % dn/dt
                dy2dt(8)=(1-y(8)).*alpha(1)-y(8).*beta(1); % dn/dt
                dy2dt(9)=(1-y(9)).*alpha(3)-y(9).*beta(3); % dh/dt
                dy2dt(10)=(1-y(10)).*alpha(4)-y(10).*beta(4);  % ds1/dt  - notice phi!
%               dy2dt(11)=(1-y(11)).*beta(5)-y(11).*alpha(5);  % ds2/dt  - notice the switch in alpha and beta!

                noise(2)=sqrt(((1-y(2)).*alpha(2)+y(2).*beta(2)));  % dm/dt 
                noise(3)=sqrt(((1-y(3)).*alpha(2)+y(3).*beta(2)));  % dm/dt 
                noise(4)=sqrt(((1-y(4)).*alpha(2)+y(4).*beta(2)));  % dm/dt 
                noise(5)=sqrt(((1-y(5)).*alpha(1)+y(5).*beta(1))); % dn/dt
                noise(6)=sqrt(((1-y(6)).*alpha(1)+y(6).*beta(1))); % dn/dt  
                noise(7)=sqrt(((1-y(7)).*alpha(1)+y(7).*beta(1))); % dn/dt  
                noise(8)=sqrt(((1-y(8)).*alpha(1)+y(8).*beta(1))); % dn/dt  
                noise(9)=sqrt(((1-y(9)).*alpha(3)+y(9).*beta(3))); % dh/dt
                noise(10)=sqrt((((1-y(10)).*alpha(4)+y(10).*beta(4))));  % ds1/dt  - notice phi!
%                 noise(11)=sqrt(((1-y(11)).*beta(5)+y(11).*alpha(5)));  % ds2/dt  - notice the switch in alpha and beta!

                y=y+dy2dt*dt+sqrt_dt_N*randn(1,11).*noise;
                
                pulse=((dt*kk<t_0)/t_0); %delta function approximation
                dE2dt=-E/tau_inact+pulse*U_se*R;
                dR2dt=(1-R-E)/tau_rec-pulse*U_se*R;

                noise_s(1)=sqrt_dt_N_s*sqrt(pulse*U_se*R)*randn;
                noise_s(2)=sqrt_dt_N_s*sqrt(abs((1-R-E)/tau_rec))*randn; %abs() was used to avoid numerical problem of imaginary numbers
                noise_s(3)=sqrt_dt_N_s*sqrt(abs(E/tau_inact))*randn; %abs() was used to avoid numerical problem of imaginary numbers
                
                E=E+dE2dt*dt+noise_s(1)+noise_s(3); %noise?
                R=R+dR2dt*dt-noise_s(1)+noise_s(2);  %noise?
               
                % continuous sample - remove at long simulations               
%                  v_array(kk+(pp-1)*cycle)=V;
%                  s_array(kk+(pp-1)*cycle)=y(10);
%                  E_array(kk+(pp-1)*cycle)=E;
%                  r_array(kk+(pp-1)*cycle)=R;

                %% here SO_hhx_Langevin  ends 
            end

        end     

        cell_AP(jj,ii)={AP};cell_R(jj,ii)={R};
        cell_Latency(jj,ii)={Latency};cell_Amplitude(jj,ii)={Amplitude}; %sampling arrays
        cell_stim_time(jj,ii)={stim_time};

    end
end
save('HHSTM_N_s_1e6_N_1e6_Time_55hrs_dt5e-3.mat');
% save('temp.mat');

toc

% time=(1:cycle*stim_num)*dt; %[ms]
% subplot(2,1,1)
% plotyy(time,v_array,time,s_array)
% subplot(2,1,2)
% plotyy(time,E_array,time,r_array)

time=(1:stim_num)/f_in;
plotyy(time,s1,time,R_array);
