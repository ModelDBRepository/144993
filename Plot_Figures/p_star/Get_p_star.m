clear all
close all
clc

cd ..
cd Class_files; %go to directory where all p_AP(s) functions are saved  

fontsize=16;


%% HHS


% source=Data_source('simulation','HHS','periodical','short') ;
% source.start_cut=0.3;
% source.end_cut=1;
% 
% I_array=[7.5 7.7 7.9 8.1 8.3] ;% [microamper]
% f_array=[1 5 10 12.5 15 17.5 20 22.5 25 27.5 30 35 40 45]; % [Hz]
% L_I=length(I_array); L_f=length(f_array);
% p_out=zeros(L_I,L_f);
% pred_p_out=zeros(L_I,L_f);
% 
% for ii=1:L_I
%     for ff=1:L_f
%         tic
%         Data=Data_array(source);
%         SetSourceAux(Data,[ii ff]) 
%         params=Data.data_Parameters;
%         model=Model(params);
%         p_out(ii,ff)=mean(Data.AP);
%         pred_p_out(ii,ff)=model.p_star;
%         toc
%     end
% end
% 
% f_mat=repmat(f_array,L_I,1);
% f_out=p_out./f_mat;
% pred_f_out=pred_p_out.*f_mat;
% 
% 
% save('HHS f_out','f_array','I_array','p_out','pred_p_out','f_out','pred_f_out');

%% HHSIP
% 
% 
% source=Data_source('simulation','HHSIP','periodical','short') ;
% source.start_cut=0.3;
% source.end_cut=1;
% 
% f_array= [10 20 30 40 50];
% I_array= [7.5 7.7 7.9 8.1 8.3];
% L_I=length(I_array); L_f=length(f_array);
% p_out=zeros(L_I,L_f);
% pred_p_out=zeros(L_I,L_f);
% 
% for ii=1:L_I
%     for ff=1:L_f
%         tic
%         Data=Data_array(source);
%         SetSourceAux(Data,[ii ff]) 
%         params=Data.data_Parameters;
%         model=Model(params);
%         p_out(ii,ff)=mean(Data.AP);
%         pred_p_out(ii,ff)=model.p_star;
%         toc
%     end
% end
% 
% f_mat=repmat(f_array,L_I,1);
% f_out=p_out./f_mat;
% pred_f_out=pred_p_out.*f_mat;
% 
% 
% save('HHSIP f_out','f_array','I_array','p_out','pred_p_out','f_out','pred_f_out');

%% HHMSIP


source=Data_source('simulation','HHMSIP','periodical','long') ;
source.start_cut=0.3;
source.end_cut=1;

f_array= [10 20 30 40 50];
L_f=length(f_array);
p_out=zeros(1,L_f);
pred_p_out=zeros(1,L_f);

for ff=1:L_f
    tic
    Data=Data_array(source);
    SetSourceAux(Data,[1 ff]) 
    params=Data.data_Parameters;
    model=Model(params);
    p_out(1,ff)=mean(Data.AP);
    p_out_std(1,ff)=std(Data.AP)/sqrt(length(Data.AP));
    pred_p_out(1,ff)=model.p_star;
    toc
end

f_mat=repmat(f_array,L_I,1);
f_out=p_out.*f_mat;
pred_f_out=pred_p_out.*f_mat;

cd ..
cd p_star
% save('HHMSIP f_out','f_array','p_out','p_out_std','pred_p_out','f_out','pred_f_out');

%% plot
L_I=1;
for ii=1:L_I
    plot(f_array,p_out(ii,:),'o',f_array,pred_p_out(ii,:),'-');
    hold on
end

xlabel('$f_{\mathrm{in}}$  [Hz]','interpreter','latex','Fontsize',fontsize);
ylabel('$p_{*}$  [Hz]','interpreter','latex','Fontsize',fontsize);
% ylim([0 1.2 *max(f_out)]);
% xlim([10, 40]);

width=0.84; height=0.8; %set width of margins
set(gca, 'Units','normalized','Position',[(1-width)*0.8 (1-height)*0.9 width height]);

