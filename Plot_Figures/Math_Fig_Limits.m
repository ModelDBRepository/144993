% This code plots Figure 3 - "Comparing mathematical results with CBM
% simulation when the assumptions fail to hold"

clear all;
close all;
clc

addpath(fullfile(pwd,'p_star')); %go to directory where all p_AP(s) functions are saved  
cd Class_files;

%% Define graphic details
Gspecs1=Graphic_Specs(1);
Gspecs2=Graphic_Specs(2);
Gspecs2.Label_fontsize=0.7*Gspecs1.Label_fontsize;
Gspecs2.Axes_fontsize=0.7*Gspecs1.Label_fontsize;
Gspecs2.Legend_fontsize=0.6*Gspecs1.Legend_fontsize;
Gspecs1=Gspecs2;

signal='AP';

title_position=[-0.15 1.06]; %for big figures
TitleFont=16;
%% Declare figure
figure
SetLines(Gspecs1);
SetFigureSize(Gspecs1) ;

%% p_star HHSIP
   
subplot(4,2,1) % S_Y
load('HHSIP f_out','f_array','I_array','p_out','pred_p_out','f_out','pred_f_out');
% T_array=1e3./f_array; %[ms]
plot(f_array,p_out','o')
hold off
hold on
plot(f_array,pred_p_out','-');

ylabel('$p_{*}$');
xlabel('$T^{-1}_{*}$ [Hz]');


lgnd=legend('Simulation','Approximation','Location','NorthEast');
set(lgnd,'box','off','Fontsize',Gspecs1.Legend_fontsize);
lgnd_children=get(lgnd, 'Children');
set(lgnd_children, 'Color', 'k');  %change color to black
set(lgnd_children(1), 'Marker','none','LineStyle', '-');  %lgnd_children(1) is the handle of only the last marker. lgnd_children(4) is the second to last, etc. 
set(lgnd_children(2), 'LineStyle', '-'); %lgnd_children(2) is the handle of only the last line. lgnd_children(5) is the second to last, etc. 
SetFonts(Gspecs1)
title('\bf{A}', 'Units','normalized','Position', title_position...
    ,'Fontsize',TitleFont);
set(gca,'box','on');



%% s_est HHSIP
subplot(4,2,3) % S_Y
source=Data_source('simulation','HHSIP','periodical','short');
source.start_cut=0.9;
Data=Data_array(source);
model_params=Data.data_Parameters;
model=Model(model_params);
cut=0.2;
num=1;


Y_star=mean(Data.AP);
[s_est s_hat AP_est AP t]=Predictor.KalmanFilter(model,Data);
Get_Plot.State_estimation(t,s_est,s_hat,num,cut,Gspecs1) 
ylim([-3.5 6.5]*1e-3);
ylabel('$\hat{s}_1$');
 
title('\bf{C}', 'Units','normalized','Position', title_position...
    ,'Fontsize',TitleFont);

%% HHSIP - timescale seperation

subplot(4,2,2)

source=Data_source('simulation','HHSIP','Periodical','short') ;
source.start_cut=0.3;
source.end_cut=1;
Data=Data_array(source);
params=Data.data_Parameters;
model=Model(params);
SetSourceAux(Data,[5 1]) 

Get_Plot.PSD(model,Data,signal,Gspecs2)
SetFonts(Gspecs2);

h_text=text(0,0,'$T^{-1}_{*}=10$ [Hz]', 'Units','normalized','Fontsize',Gspecs2.Legend_fontsize);

set(h_text,'Position', [0.65 0.15]);


lgnd=legend('Sim.','Approx.','Location','NorthWest');
set(lgnd,'box','off','Fontsize',Gspecs2.Legend_fontsize);


title('\bf{B}', 'Units','normalized','Position', title_position...
    ,'Fontsize',TitleFont);

%% HHSIP - attractor

subplot(4,2,4)

source=Data_source('simulation','HHSIP','Periodical','long') ;
source.start_cut=0.3;
source.end_cut=1;
Data=Data_array(source);
params=Data.data_Parameters;
model=Model(params);
% SetSourceAux(Data,[5 3]) 

Get_Plot.PSD(model,Data,signal,Gspecs2)
SetFonts(Gspecs2);



h_text=text(0,0,'$30$ [Hz]', 'Units','normalized','Fontsize',Gspecs2.Legend_fontsize);

set(h_text,'Position', [0.825 0.15]);

xlabel('$f$ [Hz]');



%% p_star HHMSIP
subplot(4,2,5) % S_Y
load('HHMSIP f_out','f_array','p_out','pred_p_out','f_out','pred_f_out');
% T_array=1e3./f_array; %[ms]
plot(f_array,p_out','ko')
hold off
hold on
plot(f_array,pred_p_out','k-');
% set(gco,'Color',color);

ylabel('$p_{*}$');
xlim([10 50]);
xlabel('$T^{-1}_{*}$ [Hz]');



lgnd=legend('Simulation','Approximation','Location','NorthEast');
set(lgnd,'box','off','Fontsize',Gspecs1.Legend_fontsize);
lgnd_children=get(lgnd, 'Children');
set(lgnd_children, 'Color', 'k');  %change color to black
set(lgnd_children(1), 'Marker','none','LineStyle', '-');  %lgnd_children(1) is the handle of only the last marker. lgnd_children(4) is the second to last, etc. 
set(lgnd_children(2), 'LineStyle', '-'); %lgnd_children(2) is the handle of only the last line. lgnd_children(5) is the second to last, etc. 
SetFonts(Gspecs1)
title('\bf{D}', 'Units','normalized','Position', title_position...
    ,'Fontsize',TitleFont);
set(gca,'box','on');


%% s_est HHMSIP
subplot(4,2,7) % S_Y
source=Data_source('simulation','HHMSIP','periodical','long');
source.start_cut=0.99;
Data=Data_array(source);
model_params=Data.data_Parameters;
model=Model(model_params);
cut=0.5;
num=1;


Y_star=mean(Data.AP);
[s_est s_hat AP_est AP t]=Predictor.KalmanFilter(model,Data);
Get_Plot.State_estimation(t,s_est,s_hat,num,cut,Gspecs1) 
ylim([-0.02 0.06])
 
title('\bf{E}', 'Units','normalized','Position', title_position...
    ,'Fontsize',TitleFont);
ylabel('$\hat{s}_1$');

%% HHMSIP - timescale seperation

subplot(4,2,6)

source=Data_source('simulation','HHMSIP','Periodical','long') ;
source.start_cut=0.3;
source.end_cut=1;
Data=Data_array(source);
params=Data.data_Parameters;
model=Model(params);
SetSourceAux(Data,[5 1]) 

Get_Plot.PSD(model,Data,'AP',Gspecs2)
SetFonts(Gspecs2);
 
title('\bf{F}', 'Units','normalized','Position', title_position...
    ,'Fontsize',TitleFont);


lgnd=legend('Sim.','Approx.','Location','NorthEast');
set(lgnd,'box','off','Fontsize',Gspecs2.Legend_fontsize);


h_text=text(0,0,'$10$ [Hz]', 'Units','normalized','Fontsize',Gspecs2.Legend_fontsize);
set(h_text,'Position', [0.02 0.15]);


%% HHMSIP - attractor

subplot(4,2,8)

source=Data_source('simulation','HHMSIP','Periodical','long') ;
source.start_cut=0.3;
source.end_cut=1;
Data=Data_array(source);
params=Data.data_Parameters;
model=Model(params);

Get_Plot.PSD(model,Data,'AP',Gspecs2)
SetFonts(Gspecs2);

h_text=text(0,0,'$30$ [Hz]', 'Units','normalized','Fontsize',Gspecs2.Legend_fontsize);
set(h_text,'Position', [0.02 0.15]);

xlabel('$f$ [Hz]');

%% Export
% 
% figure_name= ['Math_Fig_Limits.eps'];
% Get_Plot.Export2Folder(figure_name);

%%
cd ..