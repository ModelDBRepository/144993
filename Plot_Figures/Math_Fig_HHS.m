% This script plot Figure 2 - "Comparisons of simulations 
% with the mathematical results, for the stochastic HHS model."

clear all
close all
clc

addpath(fullfile(pwd,'p_star')); 
addpath(fullfile(pwd,'Predictor_Performances')); %go to directory where all p_AP(s) functions are saved  

cd Class_files;   

%%
source=Data_source('simulation','HHS','periodical','long');
Data=Data_array(source);
model_params=Data.data_Parameters;
model=Model(model_params);

source2=Data_source('reduced','HHS','periodical','long');
Data_reduced=Data_array(source2);

Gspecs1=Graphic_Specs(3);
Gspecs1.HaveXTicks=0; 
Gspecs1.Label_fontsize=15;
Gspecs1.Axes_fontsize=Gspecs1.Label_fontsize-5;

Gspecs2=Graphic_Specs(3);
h_fig=figure;
SetFonts(Gspecs1)
Gspecs2.Label_fontsize=18;
Gspecs2.Axes_fontsize=Gspecs2.Label_fontsize-5;
% SetFigureSize(Gspecs1)

%             data_type=Data.source.data_type;
%             lgnd=legend(data_type,'Approx.');

title_position=[-0.2 0.98];

%% HHS Firing rates
%go to directory where all p_AP(s) functions are saved  

subplot(3,2,1) % S_Y
load('HHS f_out','f_array','I_array','p_out','pred_p_out','f_out','pred_f_out');
% T_array=1e3./f_array; %[ms]
plot(f_array,p_out','o')
hold off
hold on
plot(f_array,pred_p_out','-');


ylabel('$p_{*}$');
xlabel('$T^{-1}_{*}$ [Hz]');


lgnd=legend('Simulation','Approximation','Location','NorthEast');
set(lgnd,'box','off','Fontsize',Gspecs1.Legend_fontsize-7);
lgnd_children=get(lgnd, 'Children');
set(lgnd_children, 'Color', 'k');  %change color to black
set(lgnd_children(1), 'Marker','none','LineStyle', '-');  %lgnd_children(1) is the handle of only the last marker. lgnd_children(4) is the second to last, etc. 
set(lgnd_children(2), 'LineStyle', '-'); %lgnd_children(2) is the handle of only the last line. lgnd_children(5) is the second to last, etc. 
SetFonts(Gspecs2)
title('\bf{A}', 'Units','normalized','Position', title_position...
    ,'Fontsize',Gspecs2.Label_fontsize);
set(gca,'box','on');

%% HHS Spectra


subplot(4,2,2) % S_Y

Get_Plot.PSD(model,Data,'AP',Gspecs1,Data_reduced); 


title('\bf{B}', 'Units','normalized','Position', title_position...
    ,'Fontsize',Gspecs2.Label_fontsize);

SetFonts(Gspecs1)

lgnd=legend('Sim.','Map','Approx.','Location','NorthWest');
set(lgnd,'box','off','Fontsize',Gspecs2.Legend_fontsize-7);

subplot(4,2,4)% S_s1  


Get_Plot.PSD(model,Data,'s',Gspecs2,Data_reduced);
xlabel('$f$ [Hz]');

SetFonts(Gspecs1)

%% Cross-Spectra 
source=Data_source('simulation','HHS','Poisson','long');
Data=Data_array(source);
model_params=Data.data_Parameters;
model=Model(model_params);

subplot(4,2,6)

% abs(S_YT)



Get_Plot.CPSD(model,Data,'AP','abs',Gspecs1)

% low=0.15;
% high=low;
% line([low high] ,[1e-5,1e-2],'Color','k','LineStyle','-.');
SetFonts(Gspecs1)
[CPSD f ~]=get_CPSD(Data,'AP','Intervals');
set(gca,'Yscale','Linear');
xlim([f(2) 0.15])
ylim([0 max(abs(CPSD))*1.2]) 
set(gca,'YTick',[0 2 4 6]*1e-3);

title('\bf{D}', 'Units','normalized','Position', title_position...
    ,'Fontsize',Gspecs2.Label_fontsize);

lgnd=legend('Sim.','Approx.','Location','NorthWest');
set(lgnd,'box','off','Fontsize',Gspecs2.Legend_fontsize-7,'Orientation','horizontal');

subplot(4,2,8)
% angle(S_YT)

Get_Plot.CPSD(model,Data,'AP','phase',Gspecs2)
xlim([f(2) 0.15])
xlabel('$f$ [Hz]');
SetFonts(Gspecs1)
ylim([-1 1]) 
% line([low high]  ,[-pi,pi],'Color','k','LineStyle','-.');


%% Predictors - state estimation
source=Data_source('simulation','HHS','periodical','short');
source.start_cut=0.2;
Data=Data_array(source);
model_params=Data.data_Parameters;
model=Model(model_params);
Gspecs1=Graphic_Specs(1);
Gspecs1.Label_fontsize=18;
Gspecs1.Axes_fontsize=Gspecs1.Label_fontsize-5;
Gspecs1.Legend_fontsize=Gspecs1.Label_fontsize-4;
cut=0.1;
num=1;

subplot(3,2,3)
% subplot(4,2,[5 7])

Y_star=mean(Data.AP)
[s_est s_hat AP_est AP t]=Predictor.KalmanFilter(model,Data);

Get_Plot.State_estimation(t,s_est,s_hat,num,cut,Gspecs1) 

title('\bf{C}', 'Units','normalized','Position', title_position...
    ,'Fontsize',Gspecs2.Label_fontsize);

% win=100; %smoothing window            
% Get_Plot.AP_estimation(t,AP_est,AP,win,cut,Gspecs1)

%% Predictors Performances     
load('Predictor_Performances_HHS_short_N_1e4_1e6_1e8_1e10_1e12');  

subplot(3,2,5)
show=[1 3 4 6 7];
x_axis=10.^N_array;
semilogx(x_axis,pred_error(show,:)','-o');
xlabel('$N$');
ylabel('$P_{error}$')
xlim([x_axis(1) x_axis(end)]);
set(gca,'XTick',x_axis)
lgnd=legend(predictor_types(show),'Location','SouthWest');
set(lgnd,'box','off','Fontsize',Gspecs2.Legend_fontsize-7);
ylim([0 0.55]);
set(gca,'YTick',(0:5)*0.1);
title('\bf{E}', 'Units','normalized','Position', title_position...
    ,'Fontsize',Gspecs2.Label_fontsize);

SetFonts(Gspecs2)

%% Export
% figure_name= ['Math_Fig2.eps'];
% Get_Plot.Export2Folder(figure_name);

%%
cd ..