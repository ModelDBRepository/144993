% Generates data file of 55 hours of Poisson stimulation of HHMS model
% (N=1e6, many different slow inactivation gating variables, increasing noise) 
% can modify alpha,M,k

clear all;
close all;
clc

% This simulation  does not require any sub-functions
% 1D slow inactivaion (gM=0)+stochasticity
% HH rates and Capacitance adjusted, Also s rates (gamma and delta)

intermittent=1; %if intermittent mode is required, set as 1, if transient mode required, set as 0
I_array=7.4;%[7.5 7.7 7.9 8.1 8.3];% [microamper]
D=0;%4e-8; % current diffusion constant [(microA)^2/milisec]
f_array=20;%[1 5 10 12.5 15 17.5 20 22.5 25 27.5 30 32.5 35 37.5 40 42.5 45]; % [Hz]
L_I=length(I_array);
T_array=1e3./f_array;   %[ms]
L_f=length(f_array);
Time=55*3600*1e3; %[msec] length of simulation


% cell_AP=cell(L_I,L_f);cell_s1=cell_AP;cell_s2=cell_AP;cell_I=cell_AP;
% cell_Latency=cell_AP;cell_Amplitude=cell_AP; %sampling arrays
% cell_stim_time=cell_AP;

%% HH parameters
N=1e6; %noise level
phi_HH=2; % make model faster
VNa=50; VK=-77;VL=-54; %[mV]
gNa=120;gK=36;gL=0.3;gM=0*gK; %[mS]   
Cm=1/phi_HH; %[microFarad] - not 1 Cm as in HH, since spike time is shorter
dt=0.01/phi_HH; %[msec] time step
t_0=1/phi_HH;  %[ms] pulse width
I_step=sqrt(D*dt);


%% Slow kinetics parameters
k=0.2; %slowing down factor of kinetic timescales
M=5;
k_v= k.^[zeros(1,9) (0:M-1)];
alpha_scaling=1.4;
sqrt_dt_N=sqrt(dt./(N.*k_v.^alpha_scaling));

w=ones(M,1)/M; 

phi_s=1/20; %s is slower, to get relaxation as in Gal 2010, and 
amp=0.0034*3;sigma=1/0.3;Vhalf=-17; % makes inactivation faster to make 
%% Main Simulation

name=['Poisson_LHHMS_I0=' num2str(I_array,2) '_alpha=' num2str(alpha_scaling) '_k=' num2str(k) '_N=1e' num2str(log10(N),1) ]; %name of data output
save_flag=f_array(1)*3600; %save every simulation time hour

tic
for jj=1:L_I    
    for ii=1:L_f

        f_in=f_array(ii);
        T=1e3/f_in; %[msec]
        stim_num=round(Time/T);
        stim_time=T*(0:stim_num-1);
        AP=zeros(stim_num,1);s1=AP;Latency=NaN*AP;Amplitude=NaN*AP;I_sample=NaN*AP; %sampling arrays
        Intervals=zeros(stim_num,1);

        RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock))); % so rand won't generate the same pattern each time...
        defaultStream = RandStream.getDefaultStream;

        %set initial conditions  and intialize arrays
        if dt==5e-4;
            load('params_7.5_7.7_7.9_8.1_8.3_dt=5e-4.mat','params');
        elseif dt==5e-3;
            load('params_7.5_7.7_7.9_8.1_8.3.mat','params');
        end
        
        if intermittent==1
            s_initial=params(3,2); %start on threshold for intermittent mode. set as 0.8953 for I0=7.9 if dt5e-3;
        else 
            s_initial=1;
        end
        
       I0=I_array(jj);
       y0=[ -66.4379 0.0446 0.0446 0.0446  0.2959 0.2959 0.2959 0.2959 0.6451 s_initial s_initial s_initial s_initial s_initial]; %initial condition when s is approximatly at steady state
      
        y=y0; %intial conditions 

        dy2dt=zeros(1,14);  dy=dy2dt;
        alpha=zeros(1,5); %recovery rates
        beta=alpha; %inactivation rates

        V_threshold=-10; %[mV] AP is defined if V crosses V_threshold upwards
        stim_flag=0;
        delay=0;
        Vmax=-inf;

        Pulse_width=round(t_0/dt);
        cycle=round(T/dt); %note that this is an approximation of the resulted T, not the real one
       

        for pp=1:stim_num
            I_sample(pp)=I0;
             s1(pp)=y(10);
            delay=-dt;
            stim_flag=1;
            current_cycle=geornd(1/cycle); %can be any distribution
            Intervals(pp)=current_cycle*dt; %for dt->0, this becomes exp distribution
            
            for kk=1:current_cycle
                if stim_flag==1  %sample
                    delay=delay+dt;
                    if (y(1)>V_threshold)
                        AP(pp)=1;
                        Vmax=y(1);
                        stim_flag=0;
                    end
                end

                if Vmax>-200
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
                 
                  alpha(1)=phi_HH*0.01*(V+55)./(1-exp(-0.1*(V+55)));%n
                  alpha(2)=phi_HH*0.1*(V+40)./(1-exp(-0.1.*(V+40)));%m
                  alpha(3)=phi_HH*0.07.*exp(-(V+65)./20); %h
                  alpha(4)=phi_s*0.001.*exp(-(V+85)./30); %s1
                 
                dy2dt(1)=(gNa.*(y(2).*y(3).*y(4)).*y(9)*(w(1)*y(10)+w(2)*y(11)+w(3)*y(12)+w(4)*y(13)+w(5).*y(14)).*(VNa-y(1))+(y(5).*y(6).*y(7).*y(8)).*gK.*(VK-y(1))+gL.*(VL-y(1))+I)./Cm;% dVm/dt
                dy2dt(2)=(1-y(2)).*alpha(2)-y(2).*beta(2);  % dm/dt
                dy2dt(3)=(1-y(3)).*alpha(2)-y(3).*beta(2);  % dm/dt  
                dy2dt(4)=(1-y(4)).*alpha(2)-y(4).*beta(2);  % dm/dt  
                dy2dt(5)=(1-y(5)).*alpha(1)-y(5).*beta(1); % dn/dt
                dy2dt(6)=(1-y(6)).*alpha(1)-y(6).*beta(1); % dn/dt
                dy2dt(7)=(1-y(7)).*alpha(1)-y(7).*beta(1); % dn/dt
                dy2dt(8)=(1-y(8)).*alpha(1)-y(8).*beta(1); % dn/dt
                dy2dt(9)=(1-y(9)).*alpha(3)-y(9).*beta(3); % dh/dt
                dy2dt(10)=(1-y(10)).*alpha(4)-y(10).*beta(4);  % ds1/dt 
                dy2dt(11)=k*((1-y(11)).*alpha(4)-y(11).*beta(4));  % ds2/dt 
                dy2dt(12)=k^2*((1-y(12)).*alpha(4)-y(12).*beta(4));  % ds3/dt  
                dy2dt(13)=k^3*((1-y(13)).*alpha(4)-y(13).*beta(4));  % ds4/dt  
                dy2dt(14)=k^4*((1-y(14)).*alpha(4)-y(14).*beta(4));  % ds5/dt  
 

                noise(2)=sqrt(abs((1-y(2)).*alpha(2)+y(2).*beta(2)));  % dm/dt 
                noise(3)=sqrt(abs((1-y(3)).*alpha(2)+y(3).*beta(2)));  % dm/dt 
                noise(4)=sqrt(abs((1-y(4)).*alpha(2)+y(4).*beta(2)));  % dm/dt 
                noise(5)=sqrt(abs((1-y(5)).*alpha(1)+y(5).*beta(1))); % dn/dt
                noise(6)=sqrt(abs((1-y(6)).*alpha(1)+y(6).*beta(1))); % dn/dt  
                noise(7)=sqrt(abs((1-y(7)).*alpha(1)+y(7).*beta(1))); % dn/dt  
                noise(8)=sqrt(abs((1-y(8)).*alpha(1)+y(8).*beta(1))); % dn/dt  
                noise(9)=sqrt(abs((1-y(9)).*alpha(3)+y(9).*beta(3))); % dh/dt
                noise(10)=sqrt(abs(1-y(10)).*alpha(4)+y(10).*beta(4));  % ds1/dt 
                noise(11)=sqrt(abs((1-y(11)).*alpha(4)+y(11).*beta(4)).*k);  % ds2/dt  %notice change!!!!
                noise(12)=sqrt(abs((1-y(12)).*alpha(4)+y(12).*beta(4)).*k^2);  % ds3/dt    %notice change!!!!
                noise(13)=sqrt(abs((1-y(13)).*alpha(4)+y(13).*beta(4)).*k^3);  % ds4/dt    %notice change!!!!
                noise(14)=sqrt(abs((1-y(14)).*alpha(4)+y(14).*beta(4)).*k^4);  % ds5/dt   %notice change!!!!


                I0=I0+I_step*(rand-0.5);
                if I0>8.3
                    I0=8.3;
                end
                y=y+ dy2dt*dt+randn(1,14).*noise.*sqrt_dt_N;    

                %% here SO_hhx_Langevin  ends 

            end

            if mod(pp,save_flag)==0 %save data every simulation hour
                save([name '_' num2str(100*(pp)/stim_num,3) '%Completed.mat']);
                if pp>save_flag
                    delete([name '_' num2str(100*(pp-save_flag)/stim_num,3) '%Completed.mat']);
                end
            end
            
        end     

    end
end

toc

delete([name '*.mat']);
save([name '.mat']);
