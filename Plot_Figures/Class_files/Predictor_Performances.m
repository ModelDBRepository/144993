% This script is used to generate data files of performances of various
% predictors used. This is later used in figures 

close all;
clear all;
clc

tic


%% HHMS - different N

% source=Data_source('simulation','HHMS','Periodical','short') ;
% % source.start_cut=0;
% % source.end_cut=1;
% data=Data_array(source);
% N_array=[6 8 10 12];
% predictor_types={ 'Mean','Moving Average', 'Oracle','Model - Kalman','Model - predict','StateSpace','ARMAx', ...
%     'AR','BJ','GreyBox','StateSpace0', 'ARMAx0' }';
% % predictor_types={ 'Mean','Oracle','Model - Kalman','ARMAx'}';
% 
% pred_error=zeros(length(predictor_types),length(N_array));
% pred_std=zeros(length(predictor_types),length(N_array));
% 
% for nn=1:length(N_array)
% 
% SetSourceAux(data,['N=1e' num2str(N_array(nn))]) 
% params=data.data_Parameters;
% model=Model(params);
% 
% for ii=1:length(predictor_types)
%     predictor=Predictor(data,model,predictor_types{ii}) ;
%     pred_error(ii,nn)=predictor.error_Prob;
%     pred_std(ii,nn)=predictor.error_Std; 
% end
% 
% end
% toc
% 
% save('Predictor_Performances_HHMS_short_N_1e6_1e8_1e10_1e12.mat','N_array','predictor_types','pred_error','pred_std');


%% HHS - different N
% source=Data_source('simulation','HHS','Periodical','short') ;
% source.start_cut=0.3;
% source.end_cut=1;
% data=Data_array(source);
% N_array=[4 6 8 10 12];
% predictor_types={ 'Mean','Moving Average', 'Oracle','Model - Kalman','Model - predict','StateSpace','ARMAx', ...
%     'AR','BJ','GreyBox','StateSpace0', 'ARMAx0' }';
% % predictor_types={ 'Mean','Oracle','Model - Kalman','ARMAx'}';
% 
% pred_error=zeros(length(predictor_types),length(N_array));
% pred_std=zeros(length(predictor_types),length(N_array));
% tic
% for nn=1:length(N_array)
% 
% SetSourceAux(data,['N=1e' num2str(N_array(nn))]) 
% params=data.data_Parameters;
% model=Model(params);
% 
% for ii=1:length(predictor_types)
%     predictor=Predictor(data,model,predictor_types{ii}) ;
%     pred_error(ii,nn)=predictor.error_Prob;
%     pred_std(ii,nn)=predictor.error_Std;
% end
% 
% end
% toc
% 
% save('Predictor_Performances_HHS_short_N_1e4_1e6_1e8_1e10_1e12.mat','N_array','predictor_types','pred_error','pred_std');


%% HHS - different I0 and f_in
% source=Data_source('simulation','HHS','Periodical','short') ;
% % source.start_cut=0;
% % source.end_cut=1;
% data=Data_array(source);
% f_array=[1 5 10 12.5 15 17.5 20 22.5 25 27.5 30 35 40 45]; %Hz
% I0_array=[7.5 7.7 7.9 8.1 8.3]; %microAmpere
% predictor_types={ 'Mean','Moving Average', 'Oracle','Model - Kalman','Model - predict','StateSpace','ARMAx', ...
%     'AR','BJ','GreyBox','StateSpace0', 'ARMAx0' }';
% % predictor_types={ 'Mean','Oracle','Model - Kalman','ARMAx'}';
% 
% pred_error=zeros(length(predictor_types),length(I0_array),length(f_array));
% pred_std=zeros(length(predictor_types),length(I0_array),length(f_array));
% min_ff=5;%for ff<5, intermittent mode is not reached for I0=8.3
% 
% tic
% for pp=1:length(I0_array) 
%     for ff=min_ff:length(f_array)
% 
%     SetSourceAux(data,[pp ff]) 
%     params=data.data_Parameters;
%     model=Model(params);
% 
%     for ii=1:length(predictor_types)
%         predictor=Predictor(data,model,predictor_types{ii}) ;
%         pred_error(ii,pp,ff)=predictor.error_Prob;
%         pred_std(ii,pp,ff)=predictor.error_Std;
%     end
%     end
% end
% toc
% 
% save('Predictor_Performances_HHS_short_I0-f_scan.mat','f_array','I0_array','min_ff','predictor_types','pred_error','pred_std');

%% HHSIP - different I0 and f_in
% source=Data_source('simulation','HHSIP','Periodical','short') ;
% % source.start_cut=0;
% % source.end_cut=1;
% data=Data_array(source);
% f_array=[10 20 30 40 50]; %Hz
% I0_array=[7.5 7.7 7.9 8.1 8.3]; %microAmpere
% predictor_types={ 'Mean','Moving Average', 'Oracle','Model - Kalman','Model - predict','StateSpace','ARMAx', ...
%     'AR','BJ','GreyBox','StateSpace0', 'ARMAx0' }';
% % predictor_types={ 'Mean','Oracle','Model - Kalman','ARMAx'}';
% 
% pred_error=zeros(length(predictor_types),length(I0_array),length(f_array));
% pred_std=zeros(length(predictor_types),length(I0_array),length(f_array));
% min_ff=1;%for ff<5, intermittent mode is not reached for I0=8.3
% 
% tic
% for pp=1:length(I0_array) 
%     for ff=min_ff:length(f_array)
% 
%     SetSourceAux(data,[pp ff]) 
%     params=data.data_Parameters;
%     model=Model(params);
% 
%     for ii=1:length(predictor_types)
%         predictor=Predictor(data,model,predictor_types{ii}) ;
%         pred_error(ii,pp,ff)=predictor.error_Prob;
%         pred_std(ii,pp,ff)=predictor.error_Std;
%     end
%     end
% end
% toc
% 
% save('Predictor_Performances_HHSIP_short_I0-f_scan.mat','f_array','I0_array','min_ff','predictor_types','pred_error','pred_std');

%% HHSIP - different N
source=Data_source('simulation','HHSIP','Periodical','short') ;
% source.start_cut=0;
% source.end_cut=1;
data=Data_array(source);
predictor_types={ 'Mean','Moving Average', 'Oracle','Model - Kalman','Model - predict','StateSpace','ARMAx', ...
    'AR','BJ','GreyBox','StateSpace0', 'ARMAx0' }';
% predictor_types={ 'Mean','Oracle','Model - Kalman','ARMAx'}';

N_array=[6 8 10 12];
pred_error=zeros(length(predictor_types),length(N_array));
pred_std=zeros(length(predictor_types),length(N_array));


tic
for nn=1:length(N_array)

    SetSourceAux(data,['N=1e' num2str(N_array(nn))]);
    params=data.data_Parameters;
    model=Model(params);

    for ii=1:length(predictor_types)
        predictor=Predictor(data,model,predictor_types{ii}) ;
        pred_error(ii,nn)=predictor.error_Prob;
        pred_std(ii,nn)=predictor.error_Std;
    end
end
toc

save('Predictor_Performances_HHSIP_short_N_1e6_1e8_1e10_1e12.mat','N_array','predictor_types','pred_error','pred_std');

%% HHMS - for 1/f stimulation
% source=Data_source('simulation','HHMS','1/f','long') ;
% source.start_cut=0.3;
% source.end_cut=1;
% data=Data_array(source);
% params=data.data_Parameters;
% model=Model(params);
% 
% predictor_types={ 'Mean','Moving Average', 'Oracle','Model - Kalman','Model - predict','StateSpace','ARMAx', ...
%     'AR','BJ','GreyBox','StateSpace0', 'ARMAx0' }';
% % predictor_types={ 'Mean','Oracle','Model - Kalman','ARMAx'}';
% 
% pred_error=zeros(length(predictor_types),1);
% pred_std=zeros(length(predictor_types),1);
% 
% tic
% 
% for ii=1:length(predictor_types)
%     predictor=Predictor(data,model,predictor_types{ii}) ;
%     pred_error(ii)=predictor.error_Prob;
%     pred_std(ii)=predictor.error_Std;
% end
% 
% toc
% 
% save('Predictor_Performances_HHMS_1_over_f.mat','predictor_types','pred_error','pred_std');
