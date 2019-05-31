classdef Graphic_Specs 
    %UNTITLED Summary of this class goes here
    
    % object used to set figure parameters
    % probably was not a good idea - don't do this next time
    
    %   Detailed explanation goes here
   
    properties
        
        Label_fontsize; 
        Axes_fontsize;
        Legend_fontsize;
        
        HaveXTicks
        YTicksres=1;
        
        LineWidth;
        LineStyle;
        LineColor;
        
        width; 
        height; %set width of margins        
       
        FigureSize;
              
    end      
    
    methods
        
       function obj=Graphic_Specs(flag)
                         
            set(0,'DefaultTextInterpreter', 'latex');   

            obj.Label_fontsize=25;
            obj.Axes_fontsize=obj.Label_fontsize-5;
            obj.Legend_fontsize=obj.Label_fontsize-5;
            obj.HaveXTicks=1;
            obj.YTicksres=1;

            obj.LineWidth=2;
            obj.LineStyle='-';
            obj.LineColor='k';
            
            obj.width=0.8; 
            obj.height=0.7; %set width of margins
            
            if computer==1
                x_shift=1;
            else
                x_shift=0;
            end
            
            %check units in PNAS
            scrsz = get(0,'ScreenSize'); %get screen size for figures
            
            switch flag %choose figure size and position
                case 1
                    obj.FigureSize=[x_shift*scrsz(3)+10 200 600 300];
                case 2
                    obj.FigureSize=[x_shift*scrsz(3)+10 200 800 300];
                case 3
                    obj.FigureSize=scrsz+[x_shift 0 0 0 ];
                otherwise
                    error('Incorrect input')
            end
            
       end %get default parameter and set default interpeter to Latex
        
       function ApplyPlot(obj,x,y,lgnd)
                         
             SetFonts(obj);
             SetLegend(obj,lgnd,y);
             SetAxes(obj,x,y,flag);
             SetLines(obj);
                        
       end
        
       function ApplyFigure(obj,x,y,lgnd)
                         
             SetFonts(obj);
             SetLegend(obj,lgnd,y);
             SetAxes(obj,x,y,flag);
             SetPlotSize(obj);
             SetFigureSize(obj);
             SetLines(obj);
                        
        end
        
       function SetFonts(obj)
             
             set(get(gca,'XLabel'),'Fontsize',obj.Label_fontsize);
             set(get(gca,'YLabel'),'Fontsize',obj.Label_fontsize);
             set(gca, 'Fontsize',obj.Axes_fontsize);             
                        
       end        
                
       function SetLegend(obj,lgnd,y)
            
             set(0,'DefaultTextInterpreter', 'latex'); 
             
             
             set(lgnd,'Box','off','Fontsize',obj.Legend_fontsize);
             
             if y(2)>y(end)  % y(1) is sometimes problamtic in PSDs...
                set(lgnd,'Location','Northeast');         
             else
                set(lgnd,'Location','Northwest');
             end
                        
       end
        
       function SetAxes(obj,x,y)           
           
           if strcmpi(get(gca,'Xscale'),'log')                             
               min_x=floor(log10(x(2))); % x(1) is sometimes problamtic in PSDs...
               max_x=ceil(log10(x(end)));               
               set(gca,'XLim',[x(2) x(end)]);
               set(gca,'XTick',10.^(min_x:max_x));
           end
           
           if ~obj.HaveXTicks
               set(gca,'XTick',[]);
           end
           
           if strcmpi(get(gca,'Yscale'),'log')                              
               min_y=floor(log10(min(y)));
               max_y=ceil(log10(max(y)));
               set(gca,'YLim',10.^([min_y, max_y]));
               vector=10.^(min_y:obj.YTicksres:(max_y-1));
               if length(vector)<4
                set(gca,'YTick',vector);
               else 
                set(gca,'YTick',vector(round(linspace(1,length(vector),3))));
               end
           end
           
           
       end
       
       function SetPlotSize(obj)
            position=[(1-obj.width)*0.825 (1-obj.height)*0.8 obj.width obj.height];
            set(gca, 'Units','normalized','Position',position);
       end

       function SetFigureSize(obj)            
            set(gcf, 'Position',obj.FigureSize);
       end
       
       function SetLines(obj)
            set(0,'Defaultlinelinewidth', obj.LineWidth );
            set(0,'Defaultlinelinestyle', obj.LineStyle );            
       end
       
       function SetLineColor(obj)
            h=get(gca,'children');
            set(h(end),'Color',obj.LineColor);
       end

    end 
      
end

