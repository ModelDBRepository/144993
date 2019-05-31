% This code plot Figure S2 - "Estimation noise in 
% the cross-power spectral density."

clear all
close all
clc

cd Class_files;

%%
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

%


%% Cross-Spectra 
source=Data_source('simulation','HHS','Poisson','long');
Data=Data_array(source);
model_params=Data.data_Parameters;
model=Model(model_params);

subplot(2,1,1)

% abs(S_YT)



Get_Plot.CPSD(model,Data,'AP','abs',Gspecs1,'shuffled')

low=0.15;
high=low;
line([low high] ,[1e-6,1e-1],'Color','k','LineStyle','-.');
SetFonts(Gspecs1)
[CPSD f ~]=get_CPSD(Data,'AP','Intervals');
% set(gca,'Yscale','Linear');
% xlim([f(2) 0.15])
ylim([1e-6 1e-1]) 
% set(gca,'YTick',[0 2 4 6]*1e-3);


lgnd=legend('Sim.','Shuffled','Approx.','Location','North');
set(lgnd,'box','off','Fontsize',Gspecs2.Legend_fontsize-7,'Orientation','horizontal');

subplot(2,1,2)
% angle(S_YT)

Get_Plot.CPSD(model,Data,'AP','phase',Gspecs2)

xlabel('$f$ [Hz]');
SetFonts(Gspecs1)
% ylim([-1 1]) 
line([low high]  ,[-pi,pi],'Color','k','LineStyle','-.');

%% Export
% figure_name= ['Math_Fig_SYT_full.eps'];
% Get_Plot.Export2Folder(figure_name);

%%
cd ..