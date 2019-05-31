classdef Predictor
    %PREDICTOR Summary of this class goes here
    % This class generates objects which are various
    % predictors from the system identification tool box
    % each on of these objects gives estimation error and PSD
    
    %   Detailed explanation goes here
    
    properties  (SetAccess = private)  

        error_Prob=0;
        error_Std=0;
        ident_model=0; %identified model
        
    end
    
    methods (Static)        
        % Get_Filter
        function idpoly_model=GetFilter(model) %data array object and model object
                %Obtain expression for system discrete-time filters- needs "Get_pred_Params" being run first
                % note This file generates unstable filters for HHMS model (due to numerical instability?)
                syms zz

                %% used parameters
                F=model.F;
                H=model.H;
                K=model.K;
                M=model.M;
                Gamma=model.Gamma;
                T= model.params_obj.T;
%                 A_star=model.A_star;
                
%                 a=model.a;
           
                NoiseVar= GetInternalParam(model,'NoiseVar');
                inv_std=sqrt(NoiseVar); %the innovation noise std

                %% Input (T_hat) -> Output Filter

                H_inv_tag=(zz+1)*eye(M)-F;  % we use z=zz+1 for numerical stability (I think... if I don't do this I get nonsense)

                T2AP_H=simplify((H.'/H_inv_tag)*Gamma );
                [num den]=numden(T2AP_H);
                poly_num=sym2poly(num);
                poly_den=sym2poly(den);
                norm=poly_num(1); %normalization factor
                T2AP_num=poly_num/norm;
                T2AP_num_roots=roots(T2AP_num)'+1; % z=zz+1
                T2AP_den=poly_den/norm;
                T2AP_den_roots=roots(T2AP_den)'+1;  % z=zz+1
                % [r_T2AP p_T2AP k_T2AP]=residue(T2AP_num,T2AP_den)

                T2AP_sys=minreal(zpk(T2AP_num_roots,T2AP_den_roots,1/T2AP_den(1),T)); 

                %% Noise (n_hat and e) -> Output Filter
                H_inv_tag=(zz+1)*eye(M)-F;   

                noise2AP_H=simplify(((H.'/H_inv_tag)*K+1));

                [num den]=numden(noise2AP_H);
                poly_num=sym2poly(num);
                poly_den=sym2poly(den);
                norm=poly_num(1); %normalization factor
                noise2AP_num=poly_num/norm;
                noise2AP_num_roots=roots(noise2AP_num)'+1; % z=zz+1
                noise2AP_den=poly_den/norm;
                noise2AP_den_roots=roots(noise2AP_den)' +1; % z=zz+1

                % [r_noise2AP p_noise2AP k_noise2AP]=residue(noise2AP_num,noise2AP_den)
                noise2AP_sys=zpk(noise2AP_num_roots,noise2AP_den_roots,inv_std,T); 

                %% Full Filter
                AP_sys=[T2AP_sys , noise2AP_sys];

                AP_sys.InputGroup.u=1;
                AP_sys.InputGroup.Noise=2;

                idpoly_model= idpoly(AP_sys); %m = idpoly(A,B,C,D,F,NoiseVariance,Ts)
        end
        
        function [s_est s_hat AP_est AP t]=KalmanFilter(model,data)
            %% Kalman Filter - Needs AP series and model parameters - according to anderson1979optimal, page 108
            flag1=1; %compensate for innacuracy in estimation of mean(AP)? 
            flag2=0; %Intialize using estimated steady state values?

            if strcmpi(data.source,'experimental');            	
                flag1=0;
            end
            
            %% Initialize
            M=model.M;
            F=model.F;
            H=model.H;
            Q=model.Q;
            S=model.S;
            R=model.R;
            Gamma=model.Gamma;
            G=eye(M);
            
            Y_hat=get_hat(data,'AP');
            T_hat=get_hat(data,'Intervals');
            s_hat=get_hat(data,'s');
            t=cumsum(data.Intervals);

            L=length(Y_hat);

            if flag1==0 
                Y_star=model.p_star;
            else 
                Y_star=mean(data.AP);
            end

            %% Initialize
            s_est=zeros(M,L);  % estimator = \hat{x}_{k+1/k} in book
            Sigma_s=GetSigma_s(model,var(T_hat));

            if flag2==0 
                s_est(:,1)=0; % is this cheating?
                P_est= Sigma_s ;   % error covariance = Sigma_{k+1/k} in book. Initialized with the steady state covariance matrix of s
            else 
                s_est(:,1)=s_hat(1,:)'; % is this cheating?
                P_est= 0*Sigma_s ;   % error covariance = Sigma_{k+1/k} in book. Initialized with the steady state covariance matrix of s
            end

            tic 
            for ii=1:(L-1)
                K=(F*P_est*H+G*S)*(H.'*P_est*H+R)^(-1);
                s_est(:,ii+1)=F*s_est(:,ii)+Gamma*T_hat(ii)+K*(Y_hat(ii)-H.'*s_est(:,ii));
                P_est=(F-K*H.')*P_est*(F-K*H.').'+G*Q*G.'+K*R*K.'-G*S*K.'-K*S.'*G ;
                % P_est=F*P_est*F.'-K*(F*P_est*H+G*S).'+G*Q*G.' ;
            end
            toc
            
            
            Y_est=(H.'*s_est)';            
            s_est=s_est';
            
            AP_est=(Y_est+Y_star)>0.5;            
            AP=data.AP;
    
            %Maybe I should use mean(AP) instead of p0?
%             Kalman_error_Prob=mean(abs(AP_est(end/2:end)-AP(end/2:end)))
            % Kalman_error_std=std(AP_est(end/2:end)-AP(end/2:end))


            %% Test Kalma'ns Relaibility - see page 300 in "Optimal state estimation" book by Simon
            % z=AP-AP_est; %Innovations
            % 1-abs(mean(z))/std(z)  %zero mean test
            % (H.'*P_est*H+R)/var(z) %variance test
            % psd_z_hat=var(z)/pi; %(H.'*P_est*H+R)
            % 
            % nd=20;
            % L_window=floor(length(AP)/nd);
            % psd_z=pwelch(z,rectwin(L_window)); 
            % plot(psd_z);

            % x=psd_z./psd_z_hat;
            % Nf=length(x);
            % chi_sqr=nd*sum(log(x).^2);
            % p=chi2cdf(chi_sqr,Nf)

            % [emp_dist, bins]=hist(log(x),10);
            % pred_dist=normpdf(bins,0,sqrt(1/nd))*(Nf*(bins(2)-bins(1)));
            % figure
            % plot(bins,emp_dist,bins,pred_dist)


            % loglog(psd_z)
            % psd_zt=cpsd(z,u);
            % loglog(psd_zt)
        end       
 
    end
    
    methods (Static, Hidden)

        function  [A,B,C,D,K,x0]=GreyBox(par,T,model) %data array object and model object
        %GREYBOX Summary of this function goes here
        %   Detailed explanation goes here 
            a=GetInternalParam(model,'a');
            A_star=GetInternalParam(model,'A_star');
            M=model.M;
            B = model.Gamma;
            C = par(1)*ones(1,M);
            A = eye(M)+a*C+A_star*T;
            D = 0;
            [junk1,junk2,K] = dare(A.',C.',model.Q,model.R,model.S,eye(M)); 
            K=K.';
            x0 =model.s0;

        end  
        
    end
    
    methods
        
        % Get_PSD
%         function obj=GetPSD(obj,predictorType,PSD_type) %data array object and model object
% 
%         end        
        
        
        % Constructor
        function obj=Predictor(data,model,predictor_type) %data array object and model object
            % Give data,model,predictor_type and get predition errors. Also, can get PSD.           
            % model  - for black box models with initial-guess, and grey box model
            
            M=model.M;
            params=model.params_obj;
            source=data.source;
            T=params.T; 
            AP=data.AP;
            Y_star=mean(AP);
            u=get_hat(data,'Intervals');
            y=get_hat(data,'AP');
            
            fit_data=iddata(y(1:round(end/2)),u(1:round(end/2)),T);       % ident object containing AP and Intervals for model fitting     
            
            switch predictor_type
%               cases 'Mean','Moving Average', 'Oracle','Model - Kalman','Model - predict', 'StateSpace','ARMAx', 'AR','BJ','GreyBox','StateSpace0', 'ARMAx0' 
              
                case 'Mean'
                    Pred=round(mean(AP))+AP*0;  
                    
                case 'Moving Average'                    
                    
                    window_size=101;
                    window=ones(1,window_size)/(window_size-1);
                    window(1)=0; % don't use current AP
                    p_AP=filter(window,1,AP);
                    Pred=round(p_AP);

                
                case 'Oracle'
                    s=data.s;
                    
                    if any(any(isnan(s)))  %if experimental data, Oracle prediction is not possible
                         Pred=NaN+AP*0;  
                    else
                        if M==size(s,2)
                            prob= model.Get_prob(params);
                            if strcmpi(source.model_type,'MHHMS')
                                Pred=round(prob(prod(s,2)));
                            elseif strcmpi(source.model_type,'HHSIP') || strcmpi(source.model_type,'HHMSIP')
                                Pred=round(prob(mean(s(:,1:end-1),2),s(:,end) ) ); 
                            else
                                Pred=round(prob(mean(s,2)));
    %                         s_hat=get_hat(data,'s');
                                % pred_oracle=round(model.H.'*s_hat+Y_star); %use linear oracle?
                            end
                        else %if not all s components are available,, Oracle prediction is not possible
                            error('Oracle predictor must have M==size(data.s,2)')
                        end
                    end
                                      
                case 'Model - Kalman'    %original model (unchanged), Kalman filter estimation
                    [~, ~, Pred, ~, ~]=Predictor.KalmanFilter(model,data) ;                 
                    
              case 'Model - predict'  % original model, unchanged, Predict function estimation
                    
                    ident_model = idss(model.F,model.Gamma,model.H.',0,model.K,model.s0,T); %#ok<*PROP> %HHMS model - note that information on noise variance is not given here, so model is incomplete....

                    
                case 'AR'
                    
                    order =[M M 0]; %[na nb nk]
                   ident_model=arx(fit_data,order);
                     
                case 'StateSpace' 
                    
                    order =5; 
                    ident_model= n4sid(fit_data,order); %maybe use other parameters here?
                    
                case 'ARMAx'
                    
                    orders = [M M M 0]; %[na nb nc nk]
                    ident_model=armax(fit_data,orders);%     
                    
                case 'BJ'
                    
                    orders = [M M M M 0];  %[nb nc nd nf nk]
                    ident_model= bj(fit_data,orders);
                        
                case 'GreyBox'
                    par=GetInternalParam(model,'w'); % set initial value of param here
                    ident_model = idgrey('Predictor.GreyBox',par,'d',model,T);
                    ident_model=pem(fit_data,ident_model);                    
              
                case 'StateSpace0'  %State spcae fitting With Initial Guess of original model
                    
                    ident_model = idss(model.F,model.Gamma,model.H.',0,model.K,model.s0,T); %HHMS model - note that information on noise variance is not given here, so model is incomplete....
                    ident_model=pem(fit_data,ident_model); 
                
                case 'ARMAx0'   %With Initial Guess - does not work well, since model is unstable
                    
                    idpoly_model=Predictor.GetFilter(model);  %fill this up
                    ident_model=init(idpoly_model);
                    ident_model=pem(fit_data,ident_model);
                    
                otherwise
                    error('unknown predictor type')
                    
            end
                
            if exist('ident_model','var');
                validate=iddata(y(round(end/2):end),u(round(end/2):end),T);  %ident object containing AP and Intervals for model validation   
                [y_hat,~,~] = predict(ident_model,validate);
                Pred=(y_hat.y+Y_star)>0.5; %round
                obj.ident_model=ident_model;
            else
                 Pred= Pred(round(end/2):end);
            end
            
            obj.error_Prob=mean(abs(Pred-AP(round(end/2):end))); %start_cut to remove initial Kalman transients 
            obj.error_Std=std(Pred-AP(round(end/2):end));
            
        end
        
    end
    
end

