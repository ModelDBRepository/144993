classdef Get_Plot
    %PLOT_GENERATOR Summary of this class goes here
    
    % Just a big bunch of functions used to generate plot from "model" and
    % "data_array" objects
    
    % change location in Export2Folder as needed
    
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
            
        function PSD(model,Data,signal,Gspecs,varargin)
            %varargin can be used to add another simulation
            
        [data_PSD f1]=get_PSD(Data,signal); 
        stim=Data.source.stim_type;

        if (strcmpi(stim,'periodical'))%||(strcmpi(stim,'Poisson'))
            [model_PSD f2]=get_PSD(model,signal,stim,f1(2));
        else
            [Intervals_PSD junk]=get_PSD(Data,'Intervals');
            [model_PSD f2]=get_PSD(model,signal,Intervals_PSD,f1(2));
        end
        
       switch signal
            case 'AP'
                letter='Y';
            case 's'
                letter='s';
           otherwise
               error('unfamiliar signal');
       end
        
        if size(varargin,2)==0
            h_plot=loglog(f1,data_PSD,f2,model_PSD);
        elseif size(varargin,2)==1
            Data2=varargin{1};
            if ~isa( Data2, 'Data_array')
                error('last input must be Data_array object!');
            end
            [data_PSD2 f3]=get_PSD(Data2,signal); 
            h_plot=loglog(f1,data_PSD,f2,model_PSD,f3,data_PSD2,':');
            set(gca,'Children',[h_plot(2) h_plot(3) h_plot(1)]);
        elseif size(varargin,2)>1
            error('too many inputs!');
        end    
       
       if strcmpi(Data.source.data_type,'experiment')
            set(h_plot(1),{'Color'},{'k'}')
       end
        ylabel(['$S_{' letter '}$ [sec]']); 
        SetAxes(Gspecs,f2,model_PSD);
        SetLines(Gspecs);

        end

        function CPSD(model,Data,signal,plot_type,Gspecs,varargin)
        % Data2 can be 'shuffled', or another data object
        
        if size(varargin,2)==0            
            Data2_type='empty';
        elseif size(varargin,2)==1
            Data2=varargin{1};
            if  isa(Data2,'Data_array')
                [data2_CPSD f3 ~]=get_CPSD(Data2,signal,'Intervals');
                Data2_type='moreData'; 
            elseif strcmpi(Data2,'shuffled')
                Data2_type='shuffled';
            else
                error('wrong last input - can only be "shuffled" or data array');
            end
        else
            error('too many inputs!');            
        end
            
        
        [data_CPSD f1 shuffled_CPSD]=get_CPSD(Data,signal,'Intervals');
        stim=Data.source.stim_type;

        if (strcmpi(stim,'periodical'))%||(strcmpi(stim,'Poisson'))
            [model_CPSD f2]=get_CPSD(model,signal,stim,f1(2));                
        else
            [Intervals_PSD junk]=get_PSD(Data,'Intervals');
            [model_CPSD f2]=get_CPSD(model,signal,Intervals_PSD,f1(2));
        end
        
        switch signal
            case 'AP'
                letter='Y';
            case 's'
                letter='s';
        end

        if strcmpi(plot_type,'abs')
            y1=abs(data_CPSD);
            y2=abs(model_CPSD);
            
            if strcmpi(Data2_type,'empty')
                h_plot=loglog(f1,y1,f2,y2);
            elseif strcmpi(Data2_type,'moreData')
                y3=abs(data2_CPSD);
                h_plot=loglog(f1,y1,f2,y2,f3,y3,':');
                set(gca,'Children',[h_plot(2) h_plot(3) h_plot(1)]);
            elseif strcmpi(Data2_type,'shuffled')
                y3=abs(shuffled_CPSD);
                h_plot=loglog(f1,y1,f2,y2,f1,y3,':');
                set(gca,'Children',[h_plot(2) h_plot(3) h_plot(1)]);
            end
            
            ylabel(['$|S_{' letter 'T}|$ $[\mathrm{sec}^2]$']);
            SetAxes(Gspecs,f2,y2);            
        elseif strcmpi(plot_type,'phase')
            y1=angle(data_CPSD);
            y2=angle(model_CPSD);
            if strcmpi(Data2_type,'moreData')
                y3=abs(data2_CPSD);
                h_plot=semilogx(f1,y1,f2,y2,f3,y3,':');
                set(gca,'Children',[h_plot(2) h_plot(3) h_plot(1)]);
            else
                h_plot=semilogx(f1,y1,f2,y2);
            end
            ylabel(['$\angle S_{' letter 'T}$']);
            SetAxes(Gspecs,f2,y2);
            ylim([-pi pi]);
        else
            error('plot type can be "abs" or "phase" only');
        end

        SetLines(Gspecs);
        
       if strcmpi(Data.source.data_type,'experiment') %color experimental data in black
            set(h_plot(1),{'Color'},{'k'}')
       end
        




        end                
            
        function State_estimation(t,s_est,s_hat,num,cut,Gspecs)
            if size(s_hat,2)>1
                s_hat=s_hat(:,1);
                s_est=s_est(:,1);
            end
            index=round(length(t)/3)+(1:round(cut*length(t)/3)); % show only part of simulation, use only second half
            
            plot(t,s_hat(:,num),t,s_est(:,num),'--')
            R_square= 1- var(s_hat(:,num)-s_est(:,num))/var(s_hat(:,num));
                        
            ylabel('$\hat{s}$');
            h_text=text(0,0,['$R^2 = $' num2str(R_square,2) ],'interpreter','latex','Fontsize',Gspecs.Legend_fontsize);    
            set(h_text, 'Units','normalized', 'Position', [0.15 0.1]);
            xlabel('time [sec]');
            lgnd=legend('Simulation','Estimate','Orientation','horizantal');
            set(lgnd,'Location','North','Box','off','Fontsize',Gspecs.Legend_fontsize);

            xlim([t(index(1)) t(index(end))]);
            ylim([0.9*min(s_hat) 1.5*max(s_hat(index))]);
            SetFonts(Gspecs);

%             width=0.85; height=0.62 %set width of margins
%             set(gca, 'Fontsize',smallfontsize,'Units','normalized','Position',[(1-width)*0.75 (1-height)*0.65 width height]);

        end
        
        function AP_estimation(t,AP_est,AP,win,cut,Gspecs)
            index=round(length(t)/2)+(1:round(cut*length(t)/2)); % show only part of simulation
            
            t=t(index);
            rate_est=smooth(AP_est,win);
            rate=smooth(AP,win);
            
            rate_est=rate_est(index);
            rate=rate(index);
            
            plot(t,rate,t,rate_est,'--')
            R_square= 1- var(rate-rate_est)/var(rate);  
            
            ylabel('$\hat{Y}$');
%             text(t(index(round(end/3))),max(rate(index))*0.85,['$R^2 = $' num2str(R_square) ],'interpreter','latex','Fontsize',Gspecs.legend_fontsize);    
    
            xlabel('time [sec]');
            lgnd=legend('Simulation','Estimate','Orientation','horizantal');
            set(lgnd,'Location','Northwest','Box','off','Fontsize',Gspecs.legend_fontsize);
% 
%             xlim([t(index(1)) t(index(end))]);
%             ylim([1.5*min(rate_est(index)) 1.5*max(rate_est(index))]);
            SetFonts(Gspecs);

%             width=0.85; height=0.62 %set width of margins
%             set(gca, 'Fontsize',smallfontsize,'Units','normalized','Position',[(1-width)*0.75 (1-height)*0.65 width height]);

        end
                
        function h=NewScalingStatistics(AP,Intervals)
            % Plot figure with all rate statistics used in Gal2010 (Fig. 6) simualation results from data file
            % Note that all axes are set just as in Gal2010 (Fig. 6) .
            % Also not that that if simulation is shorter than 55 hours, then use shorter
            % L_trend or smaller maximal window_array

            scrsz = get(0,'ScreenSize'); %get screen size for figures 
            title_position=[-0.15 0.98];
            
            h=figure('Position',[scrsz(3)*0.01 scrsz(4)*0.1 scrsz(3)*0.6 scrsz(4)*0.8]);  
             
            fontsize=14;
            legendfontsize=fontsize-1;
            Linewidth=2;
            
            T=mean(Intervals);
            t=cumsum(Intervals);
            t_AP=t(AP>0.5);
            

            set(0,'DefaultTextInterpreter', 'latex');
            set(0,'Defaultlinelinewidth', Linewidth );
            Specs_obj=Graphic_Specs(1);
            % type='power' %power law fitting

            %% Start plotting

            

            %% Rate Fluctuations
            window_array=round([10 30 100 300 ]/T);
            for ww=1:length(window_array)

            win=window_array(ww);
            [t_trend,rate_fluc]=Get_Plot.GetRateFluctuations(AP,T,win);    
            
            subplot(3,3,[1 2]);

            hold all;
            plot(t_trend,rate_fluc);

            end
            lgnd=legend('10s','30s','100s','300s','1000s','Orientation','horizontal');
            set(lgnd,'box','off','Location','South','Fontsize',legendfontsize);
            box('off');
            xlabel('Normalized time','Fontsize',fontsize);
            ylabel('Rate Fluct.[Hz]','Fontsize',fontsize);
%             title('\textbf{A}', 'Fontsize',fontsize, 'Units','normalized','Position', title_position);
            
%             
            h_title=title('\textbf{A}','Fontsize',fontsize, 'Units','normalized'); 
            set(h_title,'Position', [-0.065 0.98]);
            set(h_title,'Position', [-0.065 0.98]); %strange Matlab bug
            %% Find RTPDFs
            
            [bins_A A_dist bins_I I_dist ]= Get_Plot.GetRTPDFs(AP);

            %% Plot A_dist
            subplot(3,3,8)
            hold all;
            semilogy(bins_A,A_dist,'.');

            xlabel('Sequence length','Fontsize',fontsize);
            ylabel('Prob','Fontsize',fontsize);
            set(gca,'Xlim',[0 35],'XTick', [0 10 20 30]); 
            set(gca,'Ylim',[1e-6 1e-1],'YTick', [1e-5 1e-4 1e-3 1e-2],'YScale','log'); 
            box('off');
            title('\textbf{G}', 'Fontsize',fontsize, 'Units','normalized','Position', title_position);
            SetAxes(Specs_obj,bins_A,A_dist)
            SetLineColor(Specs_obj);
            %fitting  
%             cutoff=inf;
%           [alpha,R2]=Getfitting(bins_A,A_dist,'exp',cutoff);



            %% Plot I_dist
            subplot(3,3,9)
            hold all;
            semilogy(bins_I,I_dist,'.');
            box('off');
            xlabel('Sequence length','Fontsize',fontsize);
            ylabel('Prob','Fontsize',fontsize);
            
%             set(gca,'Xlim',[1e1 2e3],'XTick', [1e1 1e2 1e3],'XScale','log'); 
            set(gca,'Ylim',[1e-6 1e-1],'YTick', [1e-5 1e-4 1e-3 1e-2],'YScale','log'); 
%             title('\textbf{H}', 'Position', [3.5 1e-1],'Fontsize',fontsize, 'Units','normalized');
            title('\textbf{H}', 'Fontsize',fontsize, 'Units','normalized','Position', title_position);
            SetAxes(Specs_obj,bins_I,I_dist);
            SetLineColor(Specs_obj);
            %fitting  
%             cutoff=inf;
%           [alpha,R2]=Getfitting(bins_I,I_dist,type,cutoff);


            %%  Fano Factor, Allan Factor, CV
            [bin_size cv ff af]=Get_Plot.Get_CV_FF_AF(t_AP);

            subplot(3,3,3)
            hold all;
            semilogx(bin_size,cv)
            xlabel('T [sec]','Fontsize',fontsize);
            ylabel('CV','Fontsize',fontsize);
            set(gca,'Ylim',[0 1.5],'YTick', [0 0.5 1 1.5]); 
            set(gca,'Xlim',[1e-1 1e5],'XTick', [1e-1 1 1e1 1e2 1e3 1e4],'XScale','log');  
            title('\textbf{B}', 'Fontsize',fontsize, 'Units','normalized','Position', title_position);
            box('off');
            SetAxes(Specs_obj,bin_size,cv);
            SetLineColor(Specs_obj);

            subplot(3,3,6)
            hold all;
            loglog(bin_size,ff)
            xlabel('T [sec]','Fontsize',fontsize);
            ylabel('FF','Fontsize',fontsize);
            set(gca,'Xlim',[1e-1 1e5],'XTick', [1e-1 1 1e1 1e2 1e3 1e4],'XScale','log');  
            set(gca,'Ylim',[1e-1 1e5],'YTick', [1e-1 1 1e1 1e2 1e3 1e4],'YScale','log');  
            title('\textbf{E}', 'Fontsize',fontsize, 'Units','normalized','Position', title_position);
            box('off');
            SetAxes(Specs_obj,bin_size,ff);
            SetLineColor(Specs_obj);

            %fitting
%             cutoff=1e0;
%           [alpha,R2]=Getfitting(bin_size,ff,type,cutoff);

            subplot(3,3,7)
            hold all;
            loglog(bin_size,af)
            xlabel('T [sec]','Fontsize',fontsize);
            ylabel('AF','Fontsize',fontsize);
            set(gca,'Xlim',[1e-1 1e5],'XTick', [1e-1 1 1e1 1e2 1e3 1e4],'XScale','log');  
            set(gca,'Ylim',[1e-1 1e5],'YTick', [1e-1 1 1e1 1e2 1e3 1e4],'YScale','log');  
            title('\textbf{F}', 'Fontsize',fontsize, 'Units','normalized','Position', title_position);
            box('off');
            SetAxes(Specs_obj,bin_size,af);
            SetLineColor(Specs_obj);

            %fitting
%             cutoff=1e2;
%           [alpha,R2]=Getfitting(bin_size,af,type,cutoff);

            %% PSD
            subplot(3,3,4)
            hold all;
            [PSD f]=Get_Plot.BinnedPSD(t_AP);
            loglog(f,PSD);
            set(gca,'Xlim',[f(2) 2e-2],'XTick', [1e-4 1e-3 1e-2],'XScale','log');  
            set(gca,'Ylim',[1e-1 1e4],'YTick', [1e-3 1e-1 1e1 1e3],'YScale','log');  
            xlabel('f [Hz]','Fontsize',fontsize);
            ylabel('PSD','Fontsize',fontsize);
            title('\textbf{C}', 'Fontsize',fontsize, 'Units','normalized','Position', title_position);
            box('off');
            SetLineColor(Specs_obj);
            
%             SetAxes(Specs_obj,f,PSD)
            

            %fitting  
%             cutoff=1e-2;
%            [alpha,R2]=Getfitting(f,PSD,type,cutoff);


            %% DFA
            [bin_size,DFA]=Get_Plot.GetDFA(AP,T);

            subplot(3,3,5)
            hold all;
            loglog(bin_size,DFA)
            xlabel('T [sec]','Fontsize',fontsize);
            ylabel('DFA','Fontsize',fontsize);
            set(gca,'Xlim',[1e1 1e5],'XTick', [1e1 1e2 1e3 1e4],'XScale','log');  
            set(gca,'Ylim',[5e-1 1e4],'YTick', [1 1e1 1e2 1e3 1e4],'YScale','log');  
            title('\textbf{D}', 'Fontsize',fontsize, 'Units','normalized','Position', title_position);
            box('off');
            SetAxes(Specs_obj,bin_size,DFA);
            SetLineColor(Specs_obj);

            %fitting 
%             cutoff=1e3;
%             [alpha,R2]=Getfitting(bin_size,DFA,type,cutoff);

                        
       
        end  %Graph from Gal2010, inner parameters
        
        function AddScalingStatistics(AP,Intervals,h)
            % Adds more data to existing scaling statistics figure

            figure(h);
                       
            T=mean(Intervals);
            t=cumsum(Intervals);
            t_AP=t(AP>0.5);
            
            Specs_obj=Graphic_Specs(1);

             %% Rate Fluctuations
            window_array=round([10 30 100 300 ]/T);
            subplot(3,3,[1 2]);
            hold off
            hold on
            
            for ww=1:length(window_array)

            window=window_array(ww);
            [t_trend,rate_fluc]=Get_Plot.GetRateFluctuations(AP,T,window);    
                       
            plot(t_trend,rate_fluc,':');
            hold all;

            end
        
            %% Find RTPDFs
            
            [bins_A A_dist bins_I I_dist ]= Get_Plot.GetRTPDFs(AP);

            %% Plot A_dist
            subplot(3,3,8)
            hold all;
            semilogy(bins_A,A_dist,':');

                        %% Plot I_dist
            subplot(3,3,9)
            hold all;
            loglog(bins_I,I_dist,':');
            
            %%  Fano Factor, Allan Factor, CV
            [bin_size cv ff af]=Get_Plot.Get_CV_FF_AF(t_AP);

            subplot(3,3,3)
            hold all;
            semilogx(bin_size,cv,':')
            
            subplot(3,3,6)
            hold all;
            loglog(bin_size,ff,':')            

            subplot(3,3,7)
            hold all;
            loglog(bin_size,af,':')

            %% PSD
            subplot(3,3,4)
            hold all;
            [PSD f]=Get_Plot.BinnedPSD(t_AP);
            loglog(f,PSD,':');
            

            %% DFA
            [bin_size,DFA]=Get_Plot.GetDFA(AP,T);

            subplot(3,3,5)
            hold all;
            loglog(bin_size,DFA,':')           
     

        end  %Graph from Gal2010, inner parameters
        
        function Export2Folder(figure_name) %exports figure to articale figure folder
            set(gcf, 'Color', 'w');
            
            location=Data_array.GetLocation('export');

            export_fig(figure_name)  %cool function - check internet for support

            movefile(figure_name,...
            [location 'Research\Neuron\PNAS\Figures\' figure_name]);
            
        end
     
    end
    
    methods (Access=private,Hidden,Static)   
         
        function [t_trend,rate_fluc]=GetRateFluctuations(AP,T,window)
            
            L_trend=100;
            initial_index=1;
            f_in=1/T;
            rate_fluc=zeros(L_trend,1);

            for ii=1:L_trend
            rate_fluc(ii) =f_in*mean(AP(initial_index+(((ii-1)*window+1):(ii*window))))...
            -f_in*mean(AP(initial_index+((1:(L_trend*window)))));     
            end
            t_trend=1:L_trend;            
            
        end
                
        function [bins_A A_dist bins_I I_dist ]= GetRTPDFs(AP)
            markers=diff(AP);
            A2I=find(markers<0);
            I2A=find(markers>0);

            flag=1-abs(length(A2I)-length(I2A)); %are lengths equal?

            if AP(1)  
            if flag
            A_times=(A2I(2:end)-I2A(1:(end-1)));
            I_times=(I2A-A2I);
            else
            A_times=(A2I(2:end)-I2A);
            I_times=(I2A-A2I(1:(end-1)));
            end
            else
            if flag
            A_times=(A2I-I2A);
            I_times=(I2A(2:end)-A2I(1:(end-1)));
            else
            A_times=(A2I-I2A(1:(end-1)));
            I_times=(I2A(2:end)-A2I);
            end
            end
            
            bins_A=unique(A_times);
            bins_I= unique(I_times);

            A_dist=histc(A_times,bins_A);
            I_dist=histc(I_times,bins_I);

            A_dist=A_dist/sum(A_dist);
            I_dist=I_dist/sum(I_dist);
        end   
        
        function [bin_size cv ff af]=Get_CV_FF_AF(t_AP)            
            L_bin_size=100;
            bin_size = logspace(-1,3.5,L_bin_size);
            ff = zeros(L_bin_size,1);
            cv = zeros(L_bin_size,1);
            af = zeros(L_bin_size,1);

            for ii=1:L_bin_size
            bins =[0:bin_size(ii):max(t_AP)]+bin_size(ii)/2;
            Z = hist(t_AP,bins);
            Z = Z(1:end-1);% throw last bin
            dZ = diff(Z);

            af(ii) = mean(dZ.^2)/mean(Z)/2;    
            ff(ii) = var(Z,1)/mean(Z);
            cv(ii) = std(Z,1)/mean(Z);
            end
        end
        
        function [PSD f]=BinnedPSD(t_AP)  %as in gal2010
            bin_size=1;% sec;
            fs=1/bin_size; %Hz
            bins=0:bin_size:max(t_AP);
            Z = hist(t_AP,bins);
            [PSD f]=periodogram(Z-mean(Z),[],[],fs);            
        end
        
        function [bin_size,DFA]=GetDFA(AP,T)
            L_bin_size=20;
            L_AP=length(AP);
            bin_size = logspace(0,4,L_bin_size);
            DFA=zeros(L_bin_size,1);

            for ww=1:L_bin_size

            window=ceil(bin_size(ww)/T);
            L_calc=floor(L_AP/window);
            int_AP=cumsum(AP-mean(AP)); %intergrated AP timeseries
            detrended_int_AP=0*(1:L_calc*window);

            for ii=1:L_calc  %detrend
                window_int=int_AP((((ii-1)*window+1):(ii*window)))';
                AP_trend=window_int/[1:window ; ones(1,window)];      %least square fitting
                detrended_int_AP((((ii-1)*window+1):(ii*window)))=...  %remove trend
                int_AP((((ii-1)*window+1):(ii*window)))'-AP_trend*[1:window ; ones(1,window)];
            end

            DFA(ww)=std(detrended_int_AP);

            end
        end
        
        function [alpha,R2]=Getfitting(x,y,type,cutoff)

        %  power law fitting
        if strcmpi(type,'power')
        opt = fitoptions('Method','NonlinearLeastSquares',...
           'Lower',[-inf, -inf,0 ],...
           'Startpoint',[0 1 1]);
        ftype = fittype('b+c*x^a','options',opt);

        %  exponential fitting
        elseif strcmpi(type,'exp')
        opt= fitoptions('Method','NonlinearLeastSquares',...
           'Lower',[0,0 ],...
           'Startpoint',[1 1]);
        ftype = fittype('b*exp(-x*a)', 'options',opt);
        end

        indices=(x>0)&&(x>cutoff);
        xdata=x(indices)';
        ydata=(y(indices) );
        [cfun gof]= fit(xdata,ydata,ftype);
        alpha=cfun.a;
        R2=gof.rsquare;

        end
        

            
         
     end
    
end

