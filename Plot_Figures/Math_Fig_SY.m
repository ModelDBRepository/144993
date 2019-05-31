% This code plot Figure 4 - "Comparing mathematical results (green) with
% CBM simulation (blue) for various models."
clear all;
close all;
clc

cd Class_files;

%% Define graphic details
Gspecs=Graphic_Specs(1);
signal='AP';

title_position=[-0.15 1.06]; %for big figures
TitleFont=16;
%% Declare figure
figure
SetLines(Gspecs);
SetFigureSize(Gspecs) ;

%% Coupled HHS
subplot(2,2,1)

source=Data_source('simulation','Coupled_HHS','periodical','short') ;
source.start_cut=0.1;
source.end_cut=0.3;
Data=Data_array(source);
params=Data.data_Parameters;
model=Model(params);

Get_Plot.PSD(model,Data,signal,Gspecs)
SetFonts(Gspecs);

title('\bf{A}', 'Units','normalized','Position', title_position...
    ,'Fontsize',TitleFont);

xlabel('$f$ [Hz]');
%% HHMS

subplot(2,2,2)
source=Data_source('simulation','HHMS','periodical','long') ;
source.start_cut=0.3;
source.end_cut=1;
Data=Data_array(source);
SetSourceAux(Data,'alpha=1');
params=Data.data_Parameters;
model=Model(params);

Get_Plot.PSD(model,Data,signal,Gspecs)
SetFonts(Gspecs);

title('\bf{B}', 'Units','normalized','Position', title_position...
    ,'Fontsize',TitleFont);

xlabel('$f$ [Hz]');

%% HHSTM

subplot(2,2,3)
source=Data_source('simulation','HHSTM','Periodical','long') ;
source.start_cut=0.3;
source.end_cut=1;
Data=Data_array(source);
params=Data.data_Parameters;
model=Model(params);

Get_Plot.PSD(model,Data,signal,Gspecs)
SetFonts(Gspecs);

title('\bf{C}', 'Units','normalized','Position', title_position...
    ,'Fontsize',TitleFont);


xlabel('$f$ [Hz]');

%% Multiplicative HHMS

subplot(2,2,4)

source=Data_source('simulation','MHHMS','Periodical','long') ;
source.start_cut=0;
source.end_cut=1;
Data=Data_array(source);
% SetSourceAux(Data,'alpha=2');
params=Data.data_Parameters;
model=Model(params);

Get_Plot.PSD(model,Data,signal,Gspecs)
SetFonts(Gspecs);

title('\bf{D}', 'Units','normalized','Position', title_position...
    ,'Fontsize',TitleFont);

xlabel('$f$ [Hz]');


%% Export
% 
% figure_name= ['Math_Fig_SY.eps'];
% Get_Plot.Export2Folder(figure_name);

%%
cd ..