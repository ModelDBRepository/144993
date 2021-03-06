% Generates data file of intermittent mode in a stochastic HHS model with N=1e6

clear all;
close all;
clc

% This simulation  does not require any sub-functions
% 1D slow inactivaion (gM=0)+stochasticity
% HH rates and Capacitance adjusted, Also s rates (gamma and delta)

load('params_7.5_7.7_7.9_8.1_8.3_dt=5e-4.mat');

intermittent=1; %if intermittent mode is required, set as 1, if transient mode required, set as 0
I_array=[7.5 7.7 7.9 8.1 8.3] ;% [microamper]
f_array=[1 5 10 12.5 15 17.5 20 22.5 25 27.5 30 35 40 45]; % [Hz]
L_I=length(I_array);
T_array=1e3./f_array;   %[ms]
L_f=length(f_array);
Time=1e6; %[msec] length of simulation

cell_AP=cell(L_I,L_f);cell_s1=cell_AP;cell_s2=cell_AP;
cell_Latency=cell_AP;cell_Amplitude=cell_AP; %sampling arrays
cell_stim_time=cell_AP;

%% HH parameters
N=1e10; %noise level
phi_HH=2; % make model faster
VNa=50; VK=-77;VL=-54; %[mV]
gNa=120;gK=36;gL=0.3;gM=0;%0.01; %[mS]   
Cm=1/phi_HH; %[microFarad] - not 1 Cm as in HH, since spike time is shorter
dt=0.001/phi_HH; %[msec] time step
t_0=1/phi_HH;  %[ms] pulse width

%% Slow kinetics parameters
phi_s=1/20; %s is slower, to get relaxation as in Gal 2010, and 
amp=0.0034*3;sigma=1/0.3;Vhalf=-17; % makes inactivation faster to make
sigma2=15;
sqrt_dt_N=sqrt(dt/N);
                

%% Main Simulation

name=['SHHS_I=7.9_f=20_N_1e10_Time_1e6_dt5e-4' ]; %name of data output
save_flag=f_array(1)*3600; %save every simulation time hour

tic
for jj=3:3
    I0=I_array(jj); % [microamper] 
    for ii=7:7

        f_in=f_array(ii);
        T=1e3/f_in; %[msec]
        stim_num=round(Time/T);
        stim_time=T*(0:stim_num-1);
        AP=zeros(stim_num,1);s1=AP;s2=AP;Latency=NaN*AP;Amplitude=NaN*AP; %sampling arrays

        RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock))); % so rand won't generate the same pattern each time...
        defaultStream = RandStream.getDefaultStream;

        %set initial conditions  and intialize arrays
        if intermittent==1
        s_initial=params(jj,2)+0.001; %start just above threshold for intermittent mode
        else 
            s_initial=1;
        end
        I0=I_array(jj);
        y0=[ -66.4379 0.0446 0.0446 0.0446  0.2959 0.2959 0.2959 0.2959 0.6451 s_initial 0.77];%initial condition when s is approximatly at steady state
        y=y0; %intial conditions 

        dy2dt=zeros(1,11);
        alpha=zeros(1,5); %recovery rates
        beta=alpha; %inactivation rates

        V_threshold=-10; %[mV] AP is defined if V crosses V_threshold upwards
        stim_flag=0;
        delay=0;
        Vmax=-inf;

        cycle=round(T/dt); %note that this is an approximation of the resulted T, not the real one
        Pulse_width=round(t_0/dt);

        for pp=1:stim_num
             s1(pp)=y(10);
             s2(pp)=y(11);
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

                I=I0*(kk<Pulse_width);

                %% SO_hhx_Langevin inserted here     
                  V=y(1);
                  beta(1)= phi_HH*0.125.*exp(-(V+65)./80); %n
                  beta(2)= phi_HH*4.*exp(-(V+65)./18); %m
                  beta(3)= phi_HH*1/(exp(-0.1*(V+35))+1); %h
                  beta(4)= phi_s*amp./(exp(-(V-Vhalf)/sigma)+1); % s1 - different then Fleidervish1996 !!
                 s2_inf=1./(1+exp(-(V+35)/10));
                 tao_s2=1e3./(3.3*exp((V+35)/sigma2)+exp(-(V+35)/20));
                 beta(5)=(1-s2_inf)./tao_s2; %s2- orignally from Yamada1989

                  alpha(1)=phi_HH*0.01*(V+55)./(1-exp(-0.1*(V+55)));%n
                  alpha(2)=phi_HH*0.1*(V+40)./(1-exp(-0.1.*(V+40)));%m
                  alpha(3)=phi_HH*0.07.*exp(-(V+65)./20); %h
                  alpha(4)=phi_s*0.001.*exp(-(V+85)./30); %s1
                  alpha(5) = s2_inf./tao_s2; %0.0034./(exp(-0.1*(V+17))+1);- orignally from Yamada1989        

                dy2dt(1)=(gNa.*(y(2).*y(3).*y(4)).*y(9).*y(10).*(VNa-y(1))+(y(5).*y(6).*y(7).*y(8)).*(gK+gM.*y(11)).*(VK-y(1))+gL.*(VL-y(1))+I)./Cm;% dVm/dt
                dy2dt(2)=(1-y(2)).*alpha(2)-y(2).*beta(2);  % dm/dt
                dy2dt(3)=(1-y(3)).*alpha(2)-y(3).*beta(2);  % dm/dt  
                dy2dt(4)=(1-y(4)).*alpha(2)-y(4).*beta(2);  % dm/dt  
                dy2dt(5)=(1-y(5)).*alpha(1)-y(5).*beta(1); % dn/dt
                dy2dt(6)=(1-y(6)).*alpha(1)-y(6).*beta(1); % dn/dt
                dy2dt(7)=(1-y(7)).*alpha(1)-y(7).*beta(1); % dn/dt
                dy2dt(8)=(1-y(8)).*alpha(1)-y(8).*beta(1); % dn/dt
                dy2dt(9)=(1-y(9)).*alpha(3)-y(9).*beta(3); % dh/dt
                dy2dt(10)=(1-y(10)).*alpha(4)-y(10).*beta(4);  % ds1/dt  - notice phi!
                dy2dt(11)=(1-y(11)).*beta(5)-y(11).*alpha(5);  % ds2/dt  - notice the switch in alpha and beta!

                noise(2)=sqrt(((1-y(2)).*alpha(2)+y(2).*beta(2)));  % dm/dt 
                noise(3)=sqrt(((1-y(3)).*alpha(2)+y(3).*beta(2)));  % dm/dt 
                noise(4)=sqrt(((1-y(4)).*alpha(2)+y(4).*beta(2)));  % dm/dt 
                noise(5)=sqrt(((1-y(5)).*alpha(1)+y(5).*beta(1))); % dn/dt
                noise(6)=sqrt(((1-y(6)).*alpha(1)+y(6).*beta(1))); % dn/dt  
                noise(7)=sqrt(((1-y(7)).*alpha(1)+y(7).*beta(1))); % dn/dt  
                noise(8)=sqrt(((1-y(8)).*alpha(1)+y(8).*beta(1))); % dn/dt  
                noise(9)=sqrt(((1-y(9)).*alpha(3)+y(9).*beta(3))); % dh/dt
                noise(10)=sqrt((((1-y(10)).*alpha(4)+y(10).*beta(4))));  % ds1/dt  - notice phi!
                noise(11)=sqrt(((1-y(11)).*beta(5)+y(11).*alpha(5)));  % ds2/dt  - notice the switch in alpha and beta!

                y=y+dy2dt*dt+sqrt_dt_N*randn(1,11).*noise;    

                %% here SO_hhx_Langevin  ends 
            end
            
            if mod(pp,save_flag)==0 %save data every simulation hour
                save([name '_' num2str(100*(pp)/stim_num,3) '%Completed.mat']);
                if pp>save_flag
                    delete([name '_' num2str(100*(pp-save_flag)/stim_num,3) '%Completed.mat']);
                end
            end

        end     

        cell_AP(jj,ii)={AP};cell_s1(jj,ii)={s1};cell_s2(jj,ii)={s2};
        cell_Latency(jj,ii)={Latency};cell_Amplitude(jj,ii)={Amplitude}; %sampling arrays
        cell_stim_time(jj,ii)={stim_time};

    end
end

toc

delete([name '*.mat']);
save([name '.mat']);
