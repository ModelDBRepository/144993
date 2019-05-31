% This code plot Figure S1 - "Fitting of p_{\mathrm{AP}}\left(s\right)=\Phi\left(\left(s-a\right)/b\right)
%  in the stochastic HHS model.."

clear all
close all
clc

cd Class_files;

%% p_s_Figure
scrsz = get(0,'ScreenSize'); %get screen size for figures  
% fontsize=16;
% smallfontsize=fontsize-8;
LineWidth=2;
set(0,'DefaultTextInterpreter', 'latex');
width=0.8; height=0.7; %set width of margins
FigureSize=[scrsz(3)+100 100 1200 600];

BigSpecs=Graphic_Specs(1);
BigSpecs.Axes_fontsize=12;
BigSpecs.Label_fontsize=20;
BigSpecs.Legend_fontsize=12;

SmallSpecs=Graphic_Specs(1);
SmallSpecs.Axes_fontsize=12;
SmallSpecs.Label_fontsize=14;
SmallSpecs.Legend_fontsize=12;

title_position_big=[-0.05 1.03]; %for big figures
title_position_small=[-0.22 1.15]; %for small figures
TitleFont=16;

addpath(fullfile(pwd,'Class_files\p_s_functions')); %go to directory where all p_AP(s) functions are saved  
          


%% Plot p(s) for different currents
N=1e6;

load(['HHS_p(S)_dt5e-4_N1e' num2str(log10(N)) '.mat']);
L_I=length(I_array);
figure('Position',FigureSize);
hold all
a=0*I_array;
b=a;

subplot(4,3,[1,2,4,5])

for kk=1:L_I
    a(kk)=p_s(kk,1); b(kk)=p_s(kk,2);
    prob=@(x) 0.5*(1+erf((x-a(kk))/(sqrt(2)*b(kk))));
    plot(s1_array,prob(s1_array));
    hold all
end

hold off
hold on

for kk=1:L_I
    plot(s1_array,squeeze(AP_dist(:,:,kk)),':');
    hold all
end

xlim([0.84,0.97]);
set(gca,'XTick',[]);
ylabel('$p_{\mathrm{AP}}$');

lgnd=legend([repmat('$I_0=',length(I_array),1) num2str(I_array',2) ...
    repmat('\mu$A',length(I_array),1)] );
set(lgnd,'Location','SouthEast','Box','off','Fontsize',BigSpecs.Legend_fontsize);
SetFonts(BigSpecs);
title('\bf{A}', 'Units','normalized','Position', title_position_big...
    ,'Fontsize',TitleFont);

%% information added from SHHS_p(S)_N1e6_I_9. file
% I_array(end+1)=9;
% a(end+1)=0.7888;
% b(end+1)=0.0039;

%----------
subplot(4,3,3)
a_fit=a(1)+(I_array-I_array(1))*((a(end)-a(1))/(I_array(end)-I_array(1)));
plot(I_array,a_fit,':')
lgnd=legend('linear fit');
set(lgnd,'Location','NorthEast','Box','off','Fontsize',SmallSpecs.Legend_fontsize);
hold all
plot(I_array,a,'bO')
ylabel('a');
xlim([I_array(1)-0.1 I_array(end)+0.1]);
title('\bf{B}', 'Units','normalized','Position', title_position_small...
    ,'Fontsize',TitleFont);




set(gca,'XTick',[]);
SetFonts(SmallSpecs);

subplot(4,3,6)
plot(I_array,b,'O')
ylabel('b');
xlabel('$I_0$ [$\mu$A]');
xlim([I_array(1)-0.1 I_array(end)+0.1]);
SetFonts(SmallSpecs);


%% Plot p(s),a,b for different N

subplot(4,3,[7,8,10,11]);

N_array=[1e4 1e5 1e6 1e8 1e10 1e12];
ii=3; %I0=7.9
a=0*N_array;
b=a;

for kk=1:length(N_array)
    load(['HHS_p(S)_dt5e-4_N1e' num2str(log10(N_array(kk))) '.mat']);
    if size(s1_array,1)==1
        %ok
    elseif size(s1_array,5)==1
        s1_array=s1_array(ii,:);
    else
        error(['problem with s1_array in fil "SHHS_p(S)_N1e' num2str(log10(N_array(kk))) '.mat"'])
    end
    plot(s1_array,squeeze(AP_dist(:,:,ii)),':');
    hold all
end

hold off
hold on

for kk=1:length(N_array)
    load(['HHS_p(S)_dt5e-4_N1e' num2str(log10(N_array(kk))) '.mat']);
    a(kk)=p_s(ii,1); b(kk)=p_s(ii,2);
    prob=@(x) 0.5*(1+erf((x-a(kk))/(sqrt(2)*b(kk))));
    if size(s1_array,1)==1
        %ok
    elseif size(s1_array,5)==1
        s1_array=s1_array(ii,:);
    else
        error(['problem with s1_array in fil "SHHS_p(S)_N1e' num2str(log10(N_array(kk))) '.mat"'])
    end
    plot(s1_array,prob(s1_array));
    hold all
end


xlim([0.84,0.97]);
xlabel('$s$');
ylabel('$p_{\mathrm{AP}}$');

lgnd=legend([repmat('$\log_{10}N=$',length(N_array),1) num2str(log10(N_array'))]);
set(lgnd,'Location','SouthEast','Box','off','Fontsize',BigSpecs.Legend_fontsize);
SetFonts(BigSpecs);
title('\bf{C}', 'Units','normalized','Position', title_position_big...
    ,'Fontsize',TitleFont);

subplot(4,3,9);
semilogx(N_array,a,'O')
ylabel('a');
xlim([N_array(1)*0.1 N_array(end)*10]);
set(gca,'XTick',[]);

SetFonts(SmallSpecs);
title('\bf{D}', 'Units','normalized','Position', title_position_small...
    ,'Fontsize',TitleFont);

subplot(4,3,12);

loglog(N_array,b(1)./sqrt(N_array/N_array(1)),':')
lgnd=legend('$b_0/\sqrt N$ fit');
set(lgnd,'Location','NorthEast','Box','off','Fontsize',SmallSpecs.Legend_fontsize);
hold all
loglog(N_array,b,'bO')
xlabel('N');
ylabel('b');

min_x=-1+floor(log10(N_array(1))); % x(1) is sometimes problamtic in PSDs...
max_x=1+ceil(log10(N_array(end)));               
set(gca,'XLim',10.^[min_x max_x]);
set(gca,'XTick',10.^(min_x:2:max_x));

min_y=floor(log10(min(b)));
max_y=ceil(log10(max(b)));
set(gca,'YLim',10.^[min_y, max_y+1],'YTick',10.^(min_y:2:max_y+1));

              

SetFonts(SmallSpecs);

% lgnd=legend('$b \propto 1/\sqrt N$ fit');
% set(lgnd,'Location','NorthEast','Box','off','Fontsize',smallfontsize);

%% Export
% 
% figure_name= ['MathFigS1.eps'];
% 
% Get_Plot.Export2Folder(figure_name);


%%
cd ..