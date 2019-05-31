classdef Data_array  < hgsetget 
    %DATA_ARRAYS Summary of this class goes here
    
    % This object is generated ("seeded") by  a "Data_source" object.
    % The Data_source points the object to the source simulation\experiment data file
    % This data object can generate power spectra and arrays for Get_Plot
    % functions.
    
    % !!!!!!!!!!!! Important!!!!!!!!!!!!
    % according to location of simulation data files
    % Update folder location location in GetLocation function
    % also update relevant data file names in update_arrays function
    % !!!!!!!!!!!! Important!!!!!!!!!!!!
    
    %   Detailed explanation goes here
    
    properties (SetAccess = private, Hidden)        
        window
        noverlap         
    end         
    
    properties (SetAccess = private)
        AP
        s   
        Latency
        Intervals
        data_Parameters
        source
    end
    
    methods %constuctor and Dependent variables
        
        function obj=Data_array(Data_source)
           obj.source=Data_source; 
           obj.data_Parameters=Parameters(20,7.7,1e6,5,0.2,1.4,5e-4,Data_source.model_type) ;%Parameters(f_in,I0,N,M,epsilon,alpha,dt) 
           obj.AP=[];
           obj.s=[];
           obj.Latency=[];
           obj.Intervals=[];
           addlistener(Data_source,'SourceChanged',@(src,evnt)update_arrays(obj,src,evnt) );           
           update_arrays(obj,obj.source,0);
        end
        
        function SetSourceAux(obj,aux)  %choose a specific kind of data
            obj.source.aux=aux;
            update_arrays(obj,obj.source,0);            
        end
        
        function prop_hat=get_hat(obj,prop)
            x=obj.(prop);
            if size(x,2)==1
                prop_hat=x-mean(x);
            else
                L=size(x,1);
                prop_hat=x-repmat(mean(x,1),L,1);
            end
        end          
                
        function [prop_PSD f]=get_PSD(obj,prop)            
            
            Update_PSD_params(obj) %set parameter values
            
            x=get_hat(obj,prop);
            if size(x,2)>1  %for s
             x=x(:,1);
            end
            
            [prop_PSD f]=pwelch(x,obj.window,obj.noverlap,[],obj.data_Parameters.f_in); 
            
            
        end  
        
        function [prop_CPSD f shuffled_CPSD]=get_CPSD(obj,prop1,prop2)
            
            
            
            Update_PSD_params(obj) %set parameter values
            
            x=get_hat(obj,prop1);           
            y=get_hat(obj,prop2);
            
            if size(x,2)>1  %for s
             x=x(:,1);
            end
            if size(y,2)>1  %for s
             y=y(:,1);
            end            
            
            shuffled_x=x(randperm(length(x)));    
            shuffled_CPSD=cpsd(shuffled_x,y,obj.window,obj.noverlap,[],obj.data_Parameters.f_in); 
            
            [prop_CPSD f]=cpsd(x,y,obj.window,obj.noverlap,[],obj.data_Parameters.f_in);
            
        end
        
        function [x_axis_time,y_axis_time,Pulses_M]=PulseMatrix(obj) 

            AP=obj.AP;
            t=cumsum(obj.Intervals);
            
            P=floor(sqrt(length(AP)));
            y_axis_time=t(1:P); %sec
            x_axis_time=t(1:P^2)/3600; %hours
            Pulses_M=reshape(AP(1:P^2),P,P);  %#ok<*PROP> % Pulses Matrix
        
        end
       
     end   
    
    methods (Access=private,Hidden)        
               
        function Truncate_data(obj)
            start_cut=obj.source.start_cut;  % remove part of beginning  of data ?
            end_cut=obj.source.end_cut;    % remove part of end of data  ?
            obj.AP=obj.AP(1+round(end*start_cut):round(end*end_cut));
            obj.Intervals=obj.Intervals(1+round(end*start_cut):round(end*end_cut));  
            obj.Latency=obj.Latency(1+round(end*start_cut):round(end*end_cut))/1e3; %also change unites to sec
            obj.Latency(isnan(obj.Latency))=mean(obj.Latency(~isnan(obj.Latency)));
            
            if size(obj.s,2)==1
                obj.s=obj.s(1+round(end*start_cut):round(end*end_cut)); 
            else
                obj.s=obj.s(1+round(end*start_cut):round(end*end_cut),:);
            end
        end       
        
        function Update_PSD_params(obj)
            obj.noverlap=0; %number samples to overlap from adjecnt intervals. Should be 0 to get good results
            L=length(obj.AP);
            
            alpha=obj.data_Parameters.alpha;            
            M=obj.data_Parameters.M;
            if (M==1 )||(alpha==0)
                nd=8; %number of intervals to average spectrum on
                L_window=floor(L/nd);  
                obj.window=2;
                % hamming - pwelch default
                obj.window=hamming(L_window);
            else %long memory
                nd=1;
                L_window=floor(L/nd);  
                %rectwin - If I want to use Behran statistic
                obj.window=rectwin(L_window);
                
                %bell taper window - recommended for long memory series in
                %robinson2003time, page 258
                obj.window=0.5*(1-cos(2*pi*((1:L_window)'+0.5)/L_window)); 
            end
            
            

        end     
        
        function update_arrays(obj,src,~)
            
            %obtain location of data directory
            location=Data_array.GetLocation('import');
             
            model_type=obj.data_Parameters.model_type; %do not delete this line!
                        
            %load data
                if check_source(src,'reduced','HHMS','periodical','long')
                    addpath([ location '\Matlab\Neuron Article\Stochastic_Paper_Figures\Reduced model\']);    
%                     load('HHMS_reduced_sim_N_1e6_k_0.2_alpha_1.5.mat','AP','s','N','f_in','M','k','alpha_scaling'); 
                    load('HHMS_reduced_sim_I0_7.7_N_1e6_k_0.2_alpha_1.4.mat','I0','AP','s','N','f_in','M','k','alpha_scaling'); 
                                       
                    epsilon=k;alpha=alpha_scaling;dt=5e-4;
                    prop_list=properties(obj.data_Parameters);
                    for ii=1:length(prop_list)
                        prop=cell2mat(prop_list(ii));
                        obj.data_Parameters.(prop)=eval(prop);    
                    end
                    
                    obj.AP=AP;
                    obj.s=s;
                    obj.Latency=s*NaN;
                    obj.Intervals=0*obj.AP+obj.data_Parameters.T;                    
                elseif check_source(src,'reduced','HHS','periodical','long')
                    addpath([ location '\Matlab\Neuron Article\Stochastic_Paper_Figures\Reduced model\']);    
                    load('HHS_reduced_sim_I0_7.9_N_1e6.mat','I0','AP','s','N','f_in','M','k','alpha_scaling'); 
                                       
                    epsilon=k;alpha=alpha_scaling;dt=5e-4;
                    prop_list=properties(obj.data_Parameters);
                    for ii=1:length(prop_list)
                        prop=cell2mat(prop_list(ii));
                        obj.data_Parameters.(prop)=eval(prop);    
                    end
                    
                    obj.AP=AP;
                    obj.s=s;
                    obj.Latency=s*NaN;
                    obj.Intervals=0*obj.AP+obj.data_Parameters.T;  
                    
                elseif check_source(src,'simulation','HHS','periodical','long')
                    % HHS Simulation - Periodical Stimuation   
                    
                    addpath([ location '\Matlab\Neuron Article\Stochastic_Paper_Figures\SHHS_intermittent']);    
                    load('SHHS_stat_O1.mat','I0','AP','s1','Latency','N','f_in','dt'); %long simulation
                    
                    M=1;epsilon=0.2;alpha=0;
                    prop_list=properties(obj.data_Parameters);
                    for ii=1:length(prop_list)
                        prop=cell2mat(prop_list(ii));
                        obj.data_Parameters.(prop)=eval(prop);    
                    end
                    
                    obj.AP=AP;
                    obj.s=s1;
                    obj.Latency=Latency;
                    obj.Intervals=0*obj.AP+obj.data_Parameters.T;
                    
               elseif check_source(src,'simulation','HHS','periodical','short') %short simulation
                    %short simulation
                   
                    addpath([ location '\Matlab\Neuron Article\Stochastic_Paper_Figures\SHHS_intermittent']);    
                    if isnumeric(src.aux) %auxilary data on I0 and f_in
                       load('SHHS_intermittent_fixed_N_1e8_Time_1e6_dt_5e-4.mat','cell_AP','cell_s1','cell_Latency','N','f_array','I_array','dt');
                       pp=src.aux(1);
                       ff=src.aux(2);
                    else
                        switch src.aux
                            case 'N=1e4'
                               load('SHHS_intermittent_fixed_N_1e4_Time_1e6_dt_5e-4.mat','cell_AP','cell_s1','cell_Latency','N','f_array','I_array','dt');
                               pp=3; ff=5;  % choose I0 (pp, from I0_array) and f_in (ff from f_array)
                            case 'N=1e6'
                               load('SHHS_intermittent_fixed_N_1e6_Time_1e6_dt_5e-4.mat','cell_AP','cell_s1','cell_Latency','N','f_array','I_array','dt');
                               pp=3; ff=7;  % choose I0 (pp, from I0_array) and f_in (ff from f_array)
                            case 'N=1e8'  %default case 'N=1e8'
                               load('SHHS_intermittent_fixed_N_1e8_Time_1e6_dt_5e-4.mat','cell_AP','cell_s1','cell_Latency','N','f_array','I_array','dt');
                               pp=3; ff=7;  % choose I0 (pp, from I0_array) and f_in (ff from f_array)
                            case 'N=1e10'  %default case 'N=1e8'
                               load('SHHS_I=7.9_f=20_N_1e10_Time_1e6_dt5e-4.mat','cell_AP','cell_s1','cell_Latency','N','f_array','I_array','dt');
                               pp=3; ff=7;  % choose I0 (pp, from I0_array) and f_in (ff from f_array)
                            case 'N=1e12'  %default case 'N=1e8'
                               load('SHHS_I=7.9_f=20_N_1e12_Time_1e6_dt5e-4.mat','cell_AP','cell_s1','cell_Latency','N','f_array','I_array','dt');
                               pp=3; ff=7;  % choose I0 (pp, from I0_array) and f_in (ff from f_array)
                           otherwise  %default case 'N=1e8'
                               load('SHHS_intermittent_fixed_N_1e6_Time_1e6_dt_5e-4.mat','cell_AP','cell_s1','cell_Latency','N','f_array','I_array','dt');
                               pp=3; ff=7;  % choose I0 (pp, from I0_array) and f_in (ff from f_array)
                        end
                    end

                    s1=cell2mat(cell_s1(pp,ff));
                    AP=cell2mat(cell_AP(pp,ff));
                    Latency=cell2mat(cell_Latency(pp,ff));
             
                    f_in=f_array(ff);
                    I0=I_array(pp);
                    M=1;epsilon=0.2;alpha=0;

                    prop_list=properties(obj.data_Parameters);
                    for ii=1:length(prop_list)
                        prop=cell2mat(prop_list(ii));
                        obj.data_Parameters.(prop)=eval(prop);    
                    end

                    obj.AP=AP;
                    obj.s=s1;
                    obj.Latency=Latency;
                    obj.Intervals=0*obj.AP+obj.data_Parameters.T;
                    

                elseif check_source(src,'simulation','HHMS','periodical','long')
                    % LHHMS Simulation - Periodical Stimuation
                    addpath([ location '\Matlab\Neuron Article\Stochastic_Paper_Figures\SHHMS\Linear HHMS\']);    
%                     load('O2_LHHMS_alpha=1.5_k=0.2_N=1e6.mat','I0','AP','s1','Latency','N','f_in','M','k','alpha_scaling','dt'); 
%                     load('O1_LHHMS_f=20_alpha=1.4_k=0.2_N=1e6.mat','I0','AP','s1','Latency','N','f_in','M','k','alpha_scaling','dt'); 

                    switch src.aux
                        
                    case 'alpha=1'
                        load('O2_LHHMS_alpha=1_k=0.2_N=1e6.mat','I0','AP','s1','Latency','N','f_in','M','k','alpha_scaling','dt'); 
                    case 'alpha=1.5'
                        load('O2_LHHMS_alpha=1.5_k=0.2_N=1e6.mat','I0','AP','s1','Latency','N','f_in','M','k','alpha_scaling','dt'); 
                    case 'alpha=2'
                        load('O2_LHHMS_alpha=2_k=0.2_N=1e6.mat','I0','AP','s1','Latency','N','f_in','M','k','alpha_scaling','dt'); 
                     otherwise
                        load('O1_LHHMS_I0=7.7_alpha=1.4_k=0.2_N=1e6.mat','I0','AP','s1','Latency','N','f_in','M','k','alpha_scaling','dt'); 
                    end
                  
                    epsilon=k; alpha=alpha_scaling;
                    
                    prop_list=properties(obj.data_Parameters);
                    for ii=1:length(prop_list)
                        prop=cell2mat(prop_list(ii));
                        obj.data_Parameters.(prop)=eval(prop);    
                    end
                    
                    obj.AP=AP;
                    obj.s=s1;
                    obj.Latency=Latency;
                    obj.Intervals=0*obj.AP+obj.data_Parameters.T;

                elseif check_source(src,'simulation','HHMS','periodical','short')
                    % LHHMS Simulation for predictors - Periodical Stimuation, small dt (5e-4),short time (6 minutes)
                    addpath([ location '\Matlab\Neuron Article\Stochastic_Paper_Figures\Predictors\']);
                    switch src.aux
                        case 'N=1e6'
                            load('LHHMS_dt=5e-4_f=20_alpha=0_k=0.2_N=1e6_T=0.1hr.mat','I0','AP','s','Latency','N','f_in','M','k','alpha_scaling','dt'); 
                        case 'N=1e8'
                            load('LHHMS_dt=5e-4_f=20_alpha=0_k=0.2_N=1e8_T=0.1hr.mat','I0','AP','s','Latency','N','f_in','M','k','alpha_scaling','dt'); 
                        case 'N=1e10'
                            load('LHHMS_dt=5e-4_f=20_alpha=0_k=0.2_N=1e10_T=0.1hr.mat','I0','AP','s','Latency','N','f_in','M','k','alpha_scaling','dt'); 
                        case 'N=1e12'
                            load('LHHMS_dt=5e-4_f=20_alpha=0_k=0.2_N=1e12_T=0.1hr.mat','I0','AP','s','Latency','N','f_in','M','k','alpha_scaling','dt');                                       
                       otherwise                                        
                    load('LHHMS_dt=5e-4_I0=7.7_alpha=1.4_k=0.2_N=1e6_T=0.1hr.mat','I0','AP','s','Latency','N','f_in','M','k','alpha_scaling','dt'); 
                    end
                    epsilon=k; alpha=alpha_scaling;
                    
                    prop_list=properties(obj.data_Parameters);
                    for ii=1:length(prop_list)
                        prop=cell2mat(prop_list(ii));
                        obj.data_Parameters.(prop)=eval(prop);    
                    end
                                        
                    obj.AP=AP;
                    obj.s=s;
                    obj.Latency=Latency;
                    obj.Intervals=0*obj.AP+obj.data_Parameters.T;
                    
                elseif check_source(src,'simulation','HHMS','Poisson','long')
                    % LHHMS Simulation - Poisson Stimuation
                    addpath([ location '\Matlab\Neuron Article\Stochastic_Paper_Figures\Poisson_SHHMS\']);
                    load('Poisson_LHHMS_f=20_alpha=0_k=0.2_N=1e8.mat','I0','AP','Latency','Intervals','s1','N','f_in','M','k','alpha_scaling','dt'); 
                    epsilon=k; alpha=alpha_scaling; Intervals=Intervals/1e3;
                    
                    prop_list=properties(obj.data_Parameters);
                    for ii=1:length(prop_list)
                        prop=cell2mat(prop_list(ii));
                        obj.data_Parameters.(prop)=eval(prop);    
                    end
                                        
                    obj.AP=AP;
                    obj.s=s1;
                    obj.Latency=Latency;
                    obj.Intervals=Intervals;
                    
                elseif check_source(src,'simulation','HHS','Poisson','long')
                    % LHHMS Simulation - Poisson Stimuation
                    addpath([ location '\Matlab\Neuron Article\Stochastic_Paper_Figures\SHHS_intermittent\']);
                    load('Poisson_HHS_f=20_N=1e6.mat','I0','AP','Latency','Intervals','s1','N','f_in','M','k','alpha_scaling','dt'); 
                    epsilon=k; alpha=alpha_scaling; Intervals=Intervals/1e3;
                    
                    prop_list=properties(obj.data_Parameters);
                    for ii=1:length(prop_list)
                        prop=cell2mat(prop_list(ii));
                        obj.data_Parameters.(prop)=eval(prop);    
                    end
                                        
                    obj.AP=AP;
                    obj.s=s1;
                    obj.Latency=Latency;
                    obj.Intervals=Intervals;
                    
                    
                    
                elseif check_source(src,'experiment','HHMS','periodical','long')
                    % Experimental Data - Periodical Stimuation
                    addpath([ location '\Matlab\Asaf Data\']);
%                     load('20090305_16B_20hz.mat');
                    load('20090305_47B_20hz.mat','X');
                    T=min(diff(X(:,1))); %Sampling Rate[Hz]  - remember to change according to data!!!!
                    Time=X(end,1); %End of Experiment time [Sec]
                    t=T:T:Time;
                                        
                    obj.AP=hist(X(:,1),t)';
                    obj.s=obj.AP*NaN;
                    obj.Latency=nan*obj.AP; %don't care for now. data is in X(:,2)
                    obj.Intervals=0*obj.AP+obj.data_Parameters.T;
                    obj.data_Parameters.f_in=1/mean(obj.Intervals);
                    
                elseif check_source(src,'experiment','HHMSIP','periodical','long')
                    % Experimental Data - Periodical Stimuation
                    addpath([ location '\Matlab\Asaf Data\']);
%                     load('20080408E_e10003_25hz.mat','X');
                    load('20080408E_e10013_25hz.mat','X');
%                     load('20090305_37A_20hz.mat','X');
%                     load('20090305_23A_20hz.mat','X');
%                     load('20090305_45A_20hz.mat','X');
%                     load('20090305_12A_20hz.mat');
                    T=min(diff(X(:,1))); %Sampling Rate[Hz]  - remember to change according to data!!!!
                    Time=X(end,1); %End of Experiment time [Sec]
                    t=T:T:Time;
                                        
                    obj.AP=hist(X(:,1),t)';
                    obj.s=obj.AP*NaN;
                    obj.Latency=nan*obj.AP; %don't care for now. data is in X(:,2)
                    obj.Intervals=0*obj.AP+obj.data_Parameters.T;
                    obj.data_Parameters.f_in=1/mean(obj.Intervals);
                    obj.data_Parameters.dt=5e-3; %so I can use HHMSIP for different currents
                    obj.data_Parameters.alpha=1.8; %so I can use HHMSIP for different currents
                    
                elseif check_source(src,'experiment','HHMS','Poisson','long')
                    % Experimental Data - Poisson Stimuation
                    addpath([ location '\PoissDataForDaniel\']);
                %     load('160408_2_poisson_isi_25_neuron43.mat');
                %     load('110408_1_poisson_isi_25_neuron38.mat');
                %     load('160408_2_poisson_isi_25_neuron23.mat');
                %     load('160408_2_poisson_isi_25_neuron11.mat');
                    load('110408_1_poisson_isi_25_neuron23.mat','ndata');
                %     load('160408_2_poisson_isi_25_neuron58.mat');
                    
                    obj.AP=ndata.response01(1:end-1)';
                    obj.s=obj.AP*NaN;
                    obj.Intervals=diff(ndata.stimtime)';
                    obj.Latency=ndata.spikelat(1:end-1); %[s];
                    obj.data_Parameters.f_in=1/mean(obj.Intervals);
                    
                elseif check_source(src,'experiment','HHMS','1/f','long')          
                    % Experimental Data - 1/f Stimuation
                    addpath([ location '\PoissDataForDaniel\1_f\']);
%                      load('e36_s63_response.mat','r','stimtime','lat');
                    load('e65_s63_response.mat','r','stimtime','lat');
%                     load('e82_s63_response.mat','r','stimtime','lat');
%                     load('e15_s63_response.mat','r','stimtime','lat');
                                        
                    obj.AP=r(1:end-1);
                    obj.s=obj.AP*NaN;
                    obj.Intervals=diff(stimtime);
                    obj.Latency=lat(1:end-1);     
                    obj.data_Parameters.f_in=1/mean(obj.Intervals);

                elseif check_source(src,'simulation','HHMS','1/f','long')  
                    % Simulation Data - 1/f Stimuation
                    addpath([ location '\Matlab\Neuron Article\Stochastic_Paper_Figures\SHHMS_1_f\']);
%                     load('LHHMS_1_f_alpha=1.5_k=0.2_N=1e6.mat','I0','AP','s','Intervals','Latency','N','f_in','M','k','alpha_scaling','dt'); 
                    load('LHHMS_dt5e-4_f_alpha=1.4_k=0.2_N=1e6.mat','I0','AP','s','Intervals','Latency','N','f_in','M','k','alpha_scaling','dt'); 
                  
                    epsilon=k; alpha=alpha_scaling;
                    f_in=f_in*1e3; 
                    
                    prop_list=properties(obj.data_Parameters);
                    for ii=1:length(prop_list)
                        prop=cell2mat(prop_list(ii));
                        obj.data_Parameters.(prop)=eval(prop);    
                    end
                                        
                    obj.AP=AP;
                    obj.s=s;
                    obj.Latency=Latency;
                    obj.Intervals=Intervals*1e-3;
                    
                elseif check_source(src,'simulation','Coupled_HHS','periodical','short')  
                    % Simulation Data - Coupled HHS
                    addpath([ location '\Matlab\Neuron Article\Stochastic_Paper_Figures\Coupled Subunits HHS\']);
%                 load('intermittent.mat','AP','s1','Latency','N','f_in','I0','dt');
                    load('Coupled_intermittent_long.mat','AP','s1','Latency','N','f_in','I0','dt');
                    epsilon=0.2; alpha=0; M=1;
                    
                    prop_list=properties(obj.data_Parameters);
                    for ii=1:length(prop_list)
                        prop=cell2mat(prop_list(ii));
                        obj.data_Parameters.(prop)=eval(prop);    
                    end
                                        
                    obj.AP=AP;
                    obj.s=s1;
                    obj.Latency=Latency;
                    obj.Intervals=0*obj.AP+obj.data_Parameters.T;
                    
                elseif check_source(src,'simulation','HHSTM','periodical','long')  
                    % Simulation Data - HHSTM 
                    addpath([ location '\Matlab\Neuron Article\Stochastic_Paper_Figures\HHSTM\']);
%                     load('HHSTM_N_1e6_Time_1e6_dt5e-4.mat','AP','s1','Latency','N','f_in','I0','dt');
                    load('HHSTM_N_s_1e6_N_1e6_Time_55hrs_dt5e-3.mat','AP','s1','Latency','N','f_in','I0','dt');
                  
                    epsilon=0.2; alpha=0; M=1;                   
                    
                    prop_list=properties(obj.data_Parameters);
                    for ii=1:length(prop_list)
                        prop=cell2mat(prop_list(ii));
                        obj.data_Parameters.(prop)=eval(prop);    
                    end
                                        
                    obj.AP=AP;
                    obj.s=s1;
                    obj.Latency=Latency;
                    obj.Intervals=0*obj.AP+obj.data_Parameters.T;
                    
                elseif check_source(src,'simulation','HHSIP','periodical','short')  
                    % Simulation Data - HHSIP 
                    addpath([ location '\Matlab\Neuron Article\Stochastic_Paper_Figures\HHSIP\']);
                    if isnumeric(src.aux) %auxilary data on I0 and f_in
                       load('HHSIP_N=1e6_f=10-50_I_7.5-8.3_T_1e7.mat','cell_AP','cell_s1','cell_s2','N','f_array','I_array','dt');                  
                       pp=src.aux(1);
                       ff=src.aux(2);
                    else
                        switch src.aux
                            case  'N=1e6'
                                load('HHSIP_N=1e6_f=10-50_I_7.5-8.3_T_1e7.mat','cell_AP','cell_s1','cell_s2','N','f_array','I_array','dt');                  
                            case  'N=1e8'
                                load('HHSIP_N=1e8_f=30_I_8.3_Time_1e6.mat','cell_AP','cell_s1','cell_s2','N','f_array','I_array','dt');                  
                            case  'N=1e10'
                                load('HHSIP_N=1e10_f=30_I_8.3_Time_1e6.mat','cell_AP','cell_s1','cell_s2','N','f_array','I_array','dt');                                
                            case  'N=1e12'
                                load('HHSIP_N=1e12_f=30_I_8.3_Time_1e6.mat','cell_AP','cell_s1','cell_s2','N','f_array','I_array','dt');                  
                            case 'Null'                            
                                load('HHSIP_N=1e6_f=10-50_I_7.5-8.3_T_1e7.mat','cell_AP','cell_s1','cell_s2','N','f_array','I_array','dt');                  
                            otherwise
                                error(' inncorrect source.aux!' );
                        end
                      pp=5; ff=3;  % choose I0 (pp, from I0_array) and f_in (ff from f_array)
                    end
                    s1=cell2mat(cell_s1(pp,ff));
                    s2=cell2mat(cell_s2(pp,ff));
                    AP=cell2mat(cell_AP(pp,ff));

                    f_in=f_array(ff);
                    I0=I_array(pp);
                    M=1;epsilon=0.2;alpha=0;

                    prop_list=properties(obj.data_Parameters);
                    for ii=1:length(prop_list)
                        prop=cell2mat(prop_list(ii));
                        obj.data_Parameters.(prop)=eval(prop);    
                    end

                    obj.AP=AP;
                    obj.s=[ s1 s2];
                    obj.Latency=AP*NaN;
                    obj.Intervals=0*obj.AP+obj.data_Parameters.T;
                        
                 elseif check_source(src,'simulation','HHSIP','periodical','long')  
                    % Simulation Data - HHSIP 
                    addpath([ location '\Matlab\Neuron Article\Stochastic_Paper_Figures\HHSIP\']);
                    load('HHSIP_N=1e6_f=30_I_8.3_Time_55hrs.mat','cell_AP','cell_s1','cell_s2','N','f_array','I_array','dt');
                   
                    if isnumeric(src.aux) %auxilary data on I0 and f_in
                        error('no auxilary parameters available for this source!')
                    end
                    pp=5; ff=3;
                    
                    s1=cell2mat(cell_s1(pp,ff));
                    s2=cell2mat(cell_s2(pp,ff));
                    AP=cell2mat(cell_AP(pp,ff));

                    f_in=f_array(ff);
                    I0=I_array(pp);
                    M=1;epsilon=0.2;alpha=0;

                    prop_list=properties(obj.data_Parameters);
                    for ii=1:length(prop_list)
                        prop=cell2mat(prop_list(ii));
                        obj.data_Parameters.(prop)=eval(prop);    
                    end

                    obj.AP=AP;
                    obj.s=[ s1 s2];
                    obj.Latency=AP*NaN;
                    obj.Intervals=0*obj.AP+obj.data_Parameters.T;

                elseif check_source(src,'simulation','HHMSIP','periodical','long')  
                    % Simulation Data - HHMSIP
                    addpath([ location '\Matlab\Neuron Article\Stochastic_Paper_Figures\HHMSIP']);
                    if ischar(src.aux)
                        if strcmpi(src.aux,'Null');
                            ff=3;
                        else
                            error('source.aux must be "Null" or numeric fo HHMSIP');
                        end
                    elseif isnumeric(src.aux)
                        ff=src.aux(2);
                    else 
                        error('source.aux must be "Null" or numeric fo HHMSIP');
                    end
                        
                    f_array=[10 20 30 40 50];
                    
                    load(['O2_HHMSIP_f=' num2str(f_array(ff)) '_alpha=2_k=0.2_N=1e6.mat'],'AP','s1','s2','Latency','N','f_in','I0','dt','M','k','alpha_scaling');
                    epsilon=k; alpha=alpha_scaling;
                    
                    prop_list=properties(obj.data_Parameters);
                    for ii=1:length(prop_list)
                        prop=cell2mat(prop_list(ii));
                        obj.data_Parameters.(prop)=eval(prop);    
                    end
                                        
                    obj.AP=AP;
                    obj.s=[ s1 s2];
                    obj.Latency=Latency;
                    obj.Intervals=0*obj.AP+obj.data_Parameters.T;
                    
                elseif check_source(src,'simulation','MHHMS','periodical','long')  
                    % Simulation Data - Multiplicative HHMS 
                    addpath([ location '\Matlab\Neuron Article\Stochastic_Paper_Figures\SHHMS\Multiplicative HHMS\']);
                    
                    switch src.aux
                    case 'alpha=2'
                        load('SHHMS_O2_k=0.2_alpha=2_I=9_N=1e6.mat','AP','s1','Latency','N','f_in','I0','dt','k');
                        alpha=2;
                    otherwise 
                        load('SHHMS_O1_N=1e8.mat','AP','s1','Latency','N','f_in','I0','dt','k');
                        alpha=1;
                    end
                    
%                     load('SHHMS_O2_k=0.2_alpha=2_I=9_N=1e6.mat','AP','s1','Latency','N','f_in','I0','dt','k');
                    epsilon=k; M=5;% change parameters with simulation (M and alpha_scaling is not defined in simulation)                  
                    
                    prop_list=properties(obj.data_Parameters);
                    for ii=1:length(prop_list)
                        prop=cell2mat(prop_list(ii));
                        obj.data_Parameters.(prop)=eval(prop);    
                    end
                                        
                    obj.AP=AP;
                    obj.s=s1;
                    obj.Latency=Latency;
                    obj.Intervals=0*obj.AP+obj.data_Parameters.T;
                
                else 
                    
                    error('This type of dataset does not exist yet... make one!');
                end 
                
                Truncate_data(obj); 
            
        end
        
        % Dependent variables - with mean removed
     end   
     
    methods (Static)
         
         function location=GetLocation(str)
            if strcmpi(str,'import')
                switch computer
                    case 'PCWIN' %office 1 computer
                        location='C:\Documents and Settings\danielso\My Documents\My Dropbox';
                    otherwise %Home computer
                        location='E:\Dropbox';
                end
            elseif strcmpi(str,'export')
                switch computer
                    case 'PCWIN' %office 1 computer
                        location='E:\Technion\';
                    otherwise %Home computer
                        location='C:\Technion\';
                end
            end
            
            
         end                   
        
     end
   
    
end

