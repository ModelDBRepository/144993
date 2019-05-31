classdef Model < hgsetget 
    %MODEL Summary of this class goes here
    
    % Model object is "seeded" by a Parameters object
    % Needs p_{AP}(s) functions (save in p_s_functions sub-directory)
    % and Average rates (save in Averaged_Rates sub-directory)
    % also, uses the "Rates" class
    % can generate relevant linear model parameters, as well as Spectra
    
    %   Detailed explanation goes here
    
    properties ( Hidden)        
        params_obj
        
        num=1000;   %default size of f vector on PSDs                
        
        %Internal parameters - can be accessed only through specific
        %function, for Predictor.GreyBox  function
        NoiseVar
        a
        A_star
        w        
    end

    properties  
        F=[]  %note - this is not the same F as in paper. Here F=I+T_{star}A_{star}+aw'   (paper does not have aw term)
        Gamma=[]
        H=[]
        Q=[]
        S=[]
        R=[]
        K=[]
        s0=[]
        M=[]
        p_star=[]
        
 
    end
      
    methods
        
       function obj=Model(params_obj)
            
            obj.params_obj=params_obj;
            UpdateModel(obj,obj.params_obj);
            
            % add listener to all model parameters
            prop_list=properties(params_obj);
            for ii=1:length(prop_list)
                prop=cell2mat(prop_list(ii));
                addlistener(obj.params_obj,prop,'PostSet',@(src,evnt)UpdateModel(obj,src,evnt) );           
            end          
        end 
        
       function [prop_PSD f]=get_PSD(obj,prop,input,f_min)     
            F=obj.F;
            w_v=obj.H;
            Q=obj.Q;
            S=obj.S;
            R=obj.R;
            d=obj.Gamma;
            M=obj.M; 
            f_in=obj.params_obj.f_in;  
            e_v=zeros(M,1); e_v(1)=1;
            
            [PSD_input f]=Input_PSD(obj,input,f_min);
            
            omega=2*pi*f/f_in;
            prop_PSD=0*omega;
        
            if strcmpi(prop,'s')
                for ii=1:length(omega)
                    H_c_inv=exp(1i*omega(ii))*eye(M)-F;                     
                    prop_PSD(ii)=(e_v.'/ H_c_inv)*(Q+d*d.'*PSD_input(ii)*f_in)*(H_c_inv'\e_v); 
                end
            elseif strcmpi(prop,'AP')
                for ii=1:length(omega)
                    H_c_inv=exp(1i*omega(ii))*eye(M)-F; 
                    vec= [ (w_v.'/H_c_inv) , 1].';
                    prop_PSD(ii)=vec'*[ (Q+d*d.'*PSD_input(ii)*f_in) S ; S.' R]*vec;
                end
            elseif strcmpi(prop,'Intervals')
                prop_PSD=PSD_input*f_in;
            else
                error('second argument must be "AP" or "s" or "Intervals"');
            end
            
            prop_PSD=2*real(prop_PSD)/f_in; %convert from discrete two sided PSD to a continuous one-sided PSD, and remove imaginary part (numeric error)
           
        end  
        
       function [prop_CPSD f]=get_CPSD(obj,prop,input,f_min)  %does not calculate S_{Ys1}          
            F=obj.F;
            w_v=obj.H;
            Q=obj.Q;
            S=obj.S;
            R=obj.R;
            d=obj.Gamma;
            M=obj.M; 
            f_in=obj.params_obj.f_in;  
            e_v=zeros(M,1); e_v(1)=1;
            
            [PSD_input f]=Input_PSD(obj,input,f_min);
            
            omega=2*pi*f/f_in;
            prop_CPSD=0*omega;
        
            if strcmpi(prop,'s')
                for ii=1:length(omega)
                    H_c_inv=exp(1i*omega(ii))*eye(M)-F;                     
                    prop_CPSD(ii)=(e_v.'/ H_c_inv)*(d*PSD_input(ii)*f_in); 
                end
            elseif strcmpi(prop,'AP')
                for ii=1:length(omega)
                    H_c_inv=exp(1i*omega(ii))*eye(M)-F;                     
                    prop_CPSD(ii)=(w_v.'/ H_c_inv)*(d*PSD_input(ii)*f_in); 
                end
            else
                error('second argument must be "AP" or "s"');
            end
            
            prop_CPSD=2*prop_CPSD/f_in; %convert from discrete two sided PSD to continuous, ones sided PSD
        end       
        
       function Sigma_s=GetSigma_s(obj,var_T) % Get steady state s noise
            M=obj.M; F=obj.F;d=obj.Gamma; Q=obj.Q;
            D=Q+d*d'*var_T;
            D_vec=reshape(D,M^2,1);
            V_vec=( eye(M^2)-kron(F,F) )\ D_vec;
            Sigma_s=reshape(V_vec,M,M);  
        end
        
       function res=GetInternalParam(obj,param_name) % Get Internal parameter
           switch param_name
               case 'NoiseVar'
                    res=obj.NoiseVar; %innovation noise variance! (?)
               case 'a'
                   res=obj.a;
               case 'A_star'
                   res=obj.A_star;
               case 'w'
                   res=obj.w;
               otherwise
                   error('Parameter must be "NoiseVar", "a", "A_star" or "w"');
           end
        end
     

    end
    
  
    methods (Access=private,Hidden)
        
        function [res f]=Input_PSD(obj,input,f_min)
            f_in=obj.params_obj.f_in;
            
            if strcmpi(input,'periodical')
%               Periodical Stimuation
                f=logspace(log10(f_min),log10(0.5*f_in),obj.num);
                res=0*f;
            elseif  strcmpi(input,'Poisson')
                % Poisson Stimuation
                f=logspace(log10(f_min),log10(0.5*f_in),obj.num);
                res=0*f+1/f_in^3;
            elseif isnumeric(input)&&size(input,2)==1
                % all others
                L_input=length(input);
                f_input=linspace(0,0.5*f_in,L_input);
                f=logspace(log10(f_min),log10(0.5*f_in),obj.num);                
                res=interp1(f_input,input,f)*0.5;
            else
                error('input must be "periodical", "Poisson" or a numeric array');
            end
        end  
        
        function UpdateModel(obj,~,~)
                   
        f_in=obj.params_obj.f_in ;%[Hz]
        M=obj.params_obj.M; 
        N=obj.params_obj.N;
        alpha=obj.params_obj.alpha;
        
        epsilon=obj.params_obj.epsilon;
        epsilon_v= (epsilon.^(0:(M-1)))';
        switch obj.params_obj.model_type
            case{'HHS','HHMS','HHSTM','Coupled_HHS','MHHMS'}
                N_v=N*epsilon_v.^alpha;
            case{'HHSIP','HHMSIP'}
                N_v=N*[epsilon_v.^alpha;1];
                M=M+1; %system size increases to M+1 in this case
            otherwise
                error('unknown model type')
        end
        
        prob =Model.Get_prob(obj.params_obj);
        [rates_0 rates_m rates_p ] = Model.Get_rates(obj.params_obj);
        [rates_star,p_star,s_star,w_v] =Model.Get_fixed_point(rates_0,rates_m,rates_p,obj.params_obj.f_in,prob,obj.params_obj.model_type);

        tau_AP=rates_p.tau_AP;
        A_0=get_A(rates_0);b_0=get_b(rates_0); 
        A_m=get_A(rates_m);b_m=get_b(rates_m);
        A_p=get_A(rates_p);b_p=get_b(rates_p);
        A_star=get_A(rates_star);
        
        b_1=b_p-b_m; A_1=A_p-A_m;
        
        a=tau_AP*(A_1*s_star-b_1);
        d=A_0*s_star-b_0;
        
        Sigma_n=get_D(rates_star,s_star,N_v)/f_in; %true only for HHMS
        sigma_e=p_star-p_star^2;  % AP generation noise %-w*u_v'*Sigma_s*u_v*w

        % Define model structure according to :
        % s_{k+1}=F*s_k+Gamma*u_k+w_k
        % y_k = H.'*s_k+v_k
        % with <w_k*w_k.'>=Q,<w_k*v_k.'>=S,<v_k*v_k.'>=R
        % with intitial conditions s0, and a few other auxilary variables
        % (M,A_star,a, f_in)

        A_star_tag=a*w_v'*f_in+A_star;

        F=eye(M)+A_star_tag/f_in; %note that this is not the same F from the paper!
        H=w_v;
        Q=a*a'*sigma_e+Sigma_n;
        S=a*sigma_e;
        R=sigma_e;
        G=eye(M);
        Gamma=d;
        if (p_star>0)&&(p_star<1)
            [P,junk2,K] = dare(F.',H,Q,R,S,G); 
            K=K.';
        else
            K=a*0;  %if neuron is deterministic, no need for innovations...
            P=Sigma_n;  %no information is gained form innovations, therefore variance=error
        end
        
        
       
        obj.F=F;
        obj.H=H;
        obj.Q=Q;
        obj.S=S;
        obj.R=R;
        obj.Gamma=Gamma;
        obj.K=K;
        obj.s0=zeros(M,1);
        obj.M=M;
        obj.p_star=p_star;
        
        %Internal parameters
        obj.NoiseVar=H.'*P*H+R; %innovation noise variance! (?)
        obj.a=a;
        obj.A_star=A_star;
        obj.w=w_v(1);
            
            
        end      
         

    end
    
    methods (Static)
    
        function prob =Get_prob(params)                   
          %this function loads the @prob function from data file
          %corressponding to model, N and dt
          
            addpath(fullfile(pwd,'p_s_functions')); %go to directory where all p_AP(s) functions are saved  
          
            N=params.N;
            I0=params.I0;
            dt=params.dt;
            model_type =params.model_type;
            
            switch model_type
                case {'HHMS','MHHMS'}
                    rapid_system='HHS';
                case 'HHMSIP'
                    rapid_system='HHSIP';
                case 'Coupled_HHS'                
                   rapid_system='Coupled_HHS_I_8.3';
                case {'HHS','HHSTM','HHSIP'};              
                    rapid_system= model_type;
                otherwise
                    error('No P(s) for this system! Need another simulation here')
            end       
            
            
            % find Phi function fit of AP distriution
            switch model_type
                case {'HHMS','MHHMS','HHS'}  
                    load([rapid_system '_p(S)_dt5e' num2str(log10(dt/5)) '_N1e' num2str(log10(N),2) '.mat'],'p_s','I_array');
                    ii=find(I_array==I0);
                    a=p_s(ii,1); b=p_s(ii,2);
                    prob=@(x) 0.5*(1+erf((x-a)/(sqrt(2)*b)));                    
                case {'Coupled_HHS'}
                    disp('I_0=8.3 microAmpere used for p(s) - Coupled_HHS')
                    load([rapid_system '_p(S)_dt5e' num2str(log10(dt/5)) '_N1e' num2str(log10(N),2) '.mat'],'p_s','I_array');
                    a=p_s(1); b=p_s(2);
                    prob=@(x) 0.5*(1+erf((x-a)/(sqrt(2)*b)));
                case {'HHSTM'}
                    disp('I_0=70 microAmpere used for p(s) - HHSTM')
                    load(['HHSTM_p(S)_dt5e' num2str(log10(dt/5)) '_Ns_1e6_N1e6.mat'],'p_s','I_array');  
                    a=p_s(1); b=p_s(2);
                    prob=@(x) 0.5*(1+erf((x-a)/(sqrt(2)*b)));
                case {'HHSIP','HHMSIP'}
                    load([rapid_system '_p(S)_dt5e' num2str(log10(dt/5)) '_N1e6.mat'],'p_s','I_array');
                    N_ratio=N/1e6; %change width of prob() if N is different than 1e6
                    ii=find(I_array==I0);
                    q1=p_s(ii,1); q2=p_s(ii,2); theta=p_s(ii,3);
                    prob=@(x,y) 0.5*(1+erf((q1*x+q2*y-theta)/sqrt(2*N_ratio)));
                otherwise
                    error('Average rates for this system! Need another simulation here')
            end    
            
            if strcmpi(model_type,'MHHMS')&&I0==9 %special case
                  load('HHS_p(S)_N1e6_I_9','p_s');
                  a=p_s(1); b=p_s(2);
                  prob=@(x) 0.5*(1+erf((x-a)/(sqrt(2)*b)));
            end
            
            
        %     load('params_7.5_7.7_7.9_8.1_8.3.mat'); %gives worse results
          
        end
        
        function [rates_0 rates_m rates_p] =Get_rates(params)

      %this function loads the average rates from data file
      %corressponding to model, and dt, and returns them via rates objects
      
            addpath(fullfile(pwd,'Averaged_Rates')); %go to directory where all average rates are saved  
      
            model_type =params.model_type;       
            I0=params.I0;
            dt=params.dt;
            epsilon=params.epsilon;
            M=params.M;
            k_v= (epsilon.^(0:(M-1)))';  %#ok<*PROP>
            
            rates_0=Rates;
            rates_m=Rates;
            rates_p=Rates;
            
            switch model_type
                case {'HHMS','MHHMS','HHS','HHSTM','Coupled_HHS'}
                    load('HHS_params_6.9-0.1-12.mat','params','I0_array'); %what dt is used here?
                    
%                 if problematic, consider switching to the files belows                    
%                     if dt==5e-3
%                         load('HHS_params_7.5_7.7_7.9_8.1_8.3_dt5e-3.mat','params','I0_array'); 
%                     elseif dt==5e-4
%                         load('HHS_params_7.5_7.7_7.9_8.1_8.3_dt5e-4.mat','params','I0_array'); 
%                     else
%                         error('No Average rates for this dt! Need another simulation here')
%                     end

                    if  strcmpi(model_type,'HHSTM')
                        ii=5; % did not do a rate simulation for I0=70 yet, but it doesn't seem very important since the rates are mostly insensitive to voltage
                    else
                        ii=Model.GetIndex(I0,I0_array);
                    end
                    params_I0=params(ii,:);  %#ok<*FNDSB> %get params for this I0                    
                   
                    delta=params_I0(3);
                    rates_p.delta=delta.*k_v;
                    rates_m.delta=delta.*k_v;
                    rates_0.delta=delta.*k_v;
                    
                    rates_p.gamma=params_I0(4).*k_v;
                    rates_m. gamma=params_I0(5).*k_v;
                    rates_0.gamma=params_I0(6).*k_v;      
                    
                    tau_AP=1e-3*params_I0(1); % [sec] timescale averaging on a single AP
                    
                    rates_p.tau_AP=tau_AP;
                    rates_m.tau_AP=tau_AP;
                    rates_0.tau_AP=tau_AP;
                                        
                case {'HHMSIP','HHSIP'}
                    if dt==5e-3
                         load('HHSIP_params_7.5_7.7_7.9_8.1_8.3_dt5e-3.mat','cell_params','I0_array');  %load this even though it is of higher resolution... no other file exist now
                    elseif dt==5e-4
                         load('HHSIP_params_7.5_7.7_7.9_8.1_8.3_dt5e-4.mat','cell_params','I0_array');  %load this even though it is of higher resolution... no other file exist now
                    else
                        error('No Average rates for this dt! Need another simulation here')
                    end
                   
                    
                    ii=Model.GetIndex(I0,I0_array);
                    params_I0=cell2mat(cell_params(ii));
                    % params_I0=[T_H,w1,delta1,gamma1_H,gamma1_M,gamma1_L,gamma1_plus,gamma1_minus ;
                    % theta,w2,delta2,gamma2_H,gamma2_M,gamma2_L,gamma2_plus,gamma2_minus];

                    delta1=params_I0(1,3);
                    delta2=params_I0(2,3);
                    
                    rates_p.delta=[delta1.*k_v;delta2];
                    rates_m.delta=rates_p.delta;
                    rates_0.delta=rates_p.delta;

                    gamma_p1=params_I0(1,4);
                    gamma_m1=params_I0(1,5);
                    gamma_01=params_I0(1,6);  

                    gamma_p2=params_I0(2,4);
                    gamma_m2=params_I0(2,5);
                    gamma_02=params_I0(2,6);  
                    
                    rates_p.gamma=[gamma_p1.*k_v;gamma_p2];
                    rates_m. gamma=[gamma_m1.*k_v;gamma_m2];
                    rates_0.gamma=[gamma_01.*k_v;gamma_02];  
                    
                    
                    tau_AP=1e-3*params_I0(1); % [sec] timescale averaging on a single AP
                    
                    rates_p.tau_AP=tau_AP;
                    rates_m.tau_AP=tau_AP;
                    rates_0.tau_AP=tau_AP;
                    

                otherwise
                    error('No Average rates for this system! Need another simulation here')
            end       
            
        end   %rates=[A_0 A_m A_p]?
        
        function [rates_star,p_star,s_star,w_v] =Get_fixed_point(rates_0,rates_m,rates_p,f_in,prob,model_type)
            %This function calculates the location of fixed point for the
            %linear approximation, and outputs p_star,s_star,w
            
            tau_AP=rates_p.tau_AP; %arbitrary - I could have used rates_m or rate_0 instead
            syms p_AP
            rates_AP= (rates_p-rates_m)*f_in*tau_AP*p_AP+rates_m*f_in*tau_AP+rates_0*(1-tau_AP*f_in); %average rates for a given p_AP
            s_inf=simplify(get_steadyState(rates_AP)); %steady state for these average rates
            M=length(s_inf);            
            u_v=1+0*s_inf; % a vector of 1's
            w_v=0*u_v; %initialize w_v
            delta_s=1e-9; %some step for derivative calculation - for some reason cannot be arbtirarly small
            %Check source model_type
            switch model_type
            %If case= HHS,HHMS,HHSTM or Coupled HHS then
                case {'HHS','HHSTM','Coupled_HHS'}
                    func=@(p_AP) p_AP-prob(subs(mean(s_inf))) ;
                    p_star=fzero(func,0.999); %p0=solve(p_AP-prob(subs(s_inf(1))))  % symbolic solution is slower
                    s_star=subs(s_inf,p_AP,p_star);
                    w= (prob(mean(s_star)+delta_s)-prob(mean(s_star)-delta_s))/(2*delta_s); %slope of p(s) at s_star
                    w_v=w;
                case {'HHMS'}
                    func=@(p_AP) p_AP-prob(subs(mean(s_inf))) ;
                    p_star=fzero(func,0.999); %p0=solve(p_AP-prob(subs(s_inf(1))))  % symbolic solution is slower
                    s_star=subs(s_inf,p_AP,p_star);
                    w= (prob(mean(s_star)+delta_s)-prob(mean(s_star)-delta_s))/(2*delta_s); %slope of p(s) at s_star
                    w_v=w*u_v/M; %notice the 1/M
            %If case = MHHMS then
                case {'MHHMS'}                    
                    func=@(p_AP) p_AP-prob(subs(prod(s_inf))) ; %notice prod for MHHMS
                    p_star=fzero(func,0.999); %p0=solve(p_AP-prob(subs(s_inf(1))))  % symbolic solution is slower
                    s_star=subs(s_inf,p_AP,p_star);
                    w= (prob(prod(s_star)+delta_s)-prob(prod(s_star)-delta_s))/(2*delta_s); %slope of p(s) at s_star
                    w_v=w*u_v*s_star(1)^(M-1); %notice the prod(s_star(1:(M-1)))
            %If case = HHSIP or HHMSIP then
                case {'HHSIP'}
                    func=@(p_AP) p_AP-prob(mean(subs(s_inf(1:end-1))),subs(s_inf(end))) ; %notice prod for MHHMS
                    p_star=fzero(func,0.999); 
                    s_star=subs(s_inf,p_AP,p_star);
                    w1= (prob(mean(s_star(1:end-1))+delta_s,s_star(end))-prob(mean(s_star(1:end-1))-delta_s,s_star(end)))/(2*delta_s); %slope of p(s) at s_star
                    w2=  (prob(mean(s_star(1:end-1)),s_star(end)+delta_s)-prob(mean(s_star(1:end-1)),s_star(end)-delta_s))/(2*delta_s); %slope of p(s) at s_star
                    w_v(1)=w1;
                    w_v(2)=w2;
                case {'HHMSIP'}
                    func=@(p_AP) p_AP-prob(mean(subs(s_inf(1:end-1))),subs(s_inf(end))) ; %notice prod for MHHMS
                    p_star=fzero(func,0.999); 
                    s_star=subs(s_inf,p_AP,p_star);
                    w1= (prob(mean(s_star(1:end-1))+delta_s,s_star(end))-prob(mean(s_star(1:end-1))-delta_s,s_star(end)))/(2*delta_s); %slope of p(s) at s_star
                    w2=  (prob(mean(s_star(1:end-1)),s_star(end)+delta_s)-prob(mean(s_star(1:end-1)),s_star(end)-delta_s))/(2*delta_s); %slope of p(s) at s_star
                    w_v(1:(M-1))=w1*u_v/(M-1); %notice the 1/M
                    w_v(M)=w2;
                otherwise
                    error('unknown source type');
            end
            
            w_v=subs(w_v); %return numeric value for w_v
            rates_star=subs(rates_AP,p_AP,p_star);
            
        end
        
        function [p_star] =Get_deterministic_p_star(params,model_type)
        %This function calculates the location of fixed point for the
        %linear approximation, and outputs p_star,s_star,w
        
        prob =Model.Get_prob(params);
        
        [rates_0 rates_m rates_p ] = Model.Get_rates(params);        

        tau_AP=rates_p.tau_AP; %arbitrary - I could have used rates_m or rate_0 instead
        syms p_AP
        rates_AP= (rates_p-rates_m)*f_in*tau_AP*p_AP+rates_m*f_in*tau_AP+rates_0*(1-tau_AP*f_in); %average rates for a given p_AP
        s_inf=simplify(get_steadyState(rates_AP)); %steady state for these average rates
        M=length(s_inf);            
        u_v=1+0*s_inf; % a vector of 1's
        w_v=0*u_v; %initialize w_v
        delta_s=1e-15; %some arbtirarly small step for derivative calculation
        %Check source model_type
        switch model_type
        %If case= HHS,HHMS,HHSTM or Coupled HHS then
            case {'HHS','HHSTM','Coupled_HHS'}
                func=@(p_AP) p_AP-prob(subs(mean(s_inf))) ;
                p_star=fzero(func,0.999); %p0=solve(p_AP-prob(subs(s_inf(1))))  % symbolic solution is slower
                s_star=subs(s_inf,p_AP,p_star);
                w= (prob(mean(s_star)+delta_s)-prob(mean(s_star)-delta_s))/(2*delta_s); %slope of p(s) at s_star
                w_v=w;
            case {'HHMS'}
                func=@(p_AP) p_AP-prob(subs(mean(s_inf))) ;
                p_star=fzero(func,0.999); %p0=solve(p_AP-prob(subs(s_inf(1))))  % symbolic solution is slower
                s_star=subs(s_inf,p_AP,p_star);
                w= (prob(mean(s_star)+delta_s)-prob(mean(s_star)-delta_s))/(2*delta_s); %slope of p(s) at s_star
                w_v=w*u_v/M; %notice the 1/M
        %If case = MHHMS then
            case {'MHHMS'}                    
                func=@(p_AP) p_AP-prob(subs(prod(s_inf))) ; %notice prod for MHHMS
                p_star=fzero(func,0.999); %p0=solve(p_AP-prob(subs(s_inf(1))))  % symbolic solution is slower
                s_star=subs(s_inf,p_AP,p_star);
                w= (prob(prod(s_star)+delta_s)-prob(prod(s_star)-delta_s))/(2*delta_s); %slope of p(s) at s_star
                w_v=w*u_v*s_star(1)^(M-1); %notice the prod(s_star(1:(M-1)))
        %If case = HHSIP or HHMSIP then
            case {'HHSIP'}
                func=@(p_AP) p_AP-prob(mean(subs(s_inf(1:end-1))),subs(s_inf(end))) ; %notice prod for MHHMS
                p_star=fzero(func,0.999); 
                s_star=subs(s_inf,p_AP,p_star);
                w1= (prob(mean(s_star(1:end-1))+delta_s,s_star(end))-prob(mean(s_star(1:end-1))-delta_s,s_star(end)))/(2*delta_s); %slope of p(s) at s_star
                w2=  (prob(mean(s_star(1:end-1)),s_star(end)+delta_s)-prob(mean(s_star(1:end-1)),s_star(end)-delta_s))/(2*delta_s); %slope of p(s) at s_star
                w_v(1)=w1;
                w_v(2)=w2;
            case {'HHMSIP'}
                func=@(p_AP) p_AP-prob(mean(subs(s_inf(1:end-1))),subs(s_inf(end))) ; %notice prod for MHHMS
                p_star=fzero(func,0.999); 
                s_star=subs(s_inf,p_AP,p_star);
                w1= (prob(mean(s_star(1:end-1))+delta_s,s_star(end))-prob(mean(s_star(1:end-1))-delta_s,s_star(end)))/(2*delta_s); %slope of p(s) at s_star
                w2=  (prob(mean(s_star(1:end-1)),s_star(end)+delta_s)-prob(mean(s_star(1:end-1)),s_star(end)-delta_s))/(2*delta_s); %slope of p(s) at s_star
                w_v(1:(M-1))=w1*u_v/(M-1); %notice the 1/M
                w_v(M)=w2;
            otherwise
                error('unknown source type');
        end

        w_v=subs(w_v); %return numeric value for w_v
        rates_star=subs(rates_AP,p_AP,p_star);

        end
    
        function index=GetIndex(number,array)
            tolerance=1e-5;
            index=find(abs(number-array)<tolerance);
        end
                 
    end
        
end
    


