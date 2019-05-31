classdef Data_source < hgsetget 
    %UNTITLED Summary of this class goes here
    
    % source of simulation data - seed for "data array"
    
    %   Detailed explanation goes here
   
    properties
        data_type='simulation';% 'experiment' or 'simulation' or 'reduced'
        model_type='HHS'; %'HHS' or 'HHMS' 
        stim_type='periodical'; % 'periodical' or 'Poisson' or '1/f' or 'sinusoids'        
        length_type='long'; %'long' or 'short'
        
        start_cut=0;  % remove part of beginning  of data ?
        end_cut=1;    % remove part of end of data  ?
    end    
    
    properties (Hidden)
        aux='Null'; %choose a specific kind of data, should be changed only by Data_array.SetSourceAux function
    end
           
    events
        SourceChanged;
    end
    
    methods
        
        %constructor
        function obj = Data_source (data_type,model_type,stim_type,length_type)
            obj.data_type=data_type;
            obj.model_type=model_type;
            obj.stim_type=stim_type;
            obj.length_type=length_type;
            
            if (strcmpi(data_type,'experiment')&&...      
                strcmpi(stim_type,'periodical'))

                obj.start_cut=0.01;  %remove part of beginning to reduce effects of transient 
                obj.end_cut=0.36;  %remove part of end where big fluctuations occur
                
            elseif (strcmpi(data_type,'experiment')&&...
                    strcmpi(stim_type,'Poisson'))
                obj.start_cut=0.6;  %remove part of beginning to reduce effects of transient 
                obj.end_cut=1;  %remove part of beginning to reduce effects of transient 
            else
                obj.start_cut=0.1;  %remove part of beginning to reduce effects of transient 
                obj.end_cut=1;  %remove part of beginning to reduce effects of transient 
            end
                
        end
        
        function res=check_source(obj,data_type,model_type,stim_type,length_type)
            res=strcmpi(data_type,obj.data_type)&&...
                        strcmpi(model_type,obj.model_type)&&...
                        strcmpi(stim_type,obj.stim_type)&&...
                         strcmpi(length_type,obj.length_type);
        end
        
        %check if input is correct and notify changes
        function set.data_type(obj,val)
         if ~(strcmpi(val,'experiment') ||... 
            strcmpi(val,'simulation') ||... 
            strcmpi(val,'reduced'))
            error('data_type must be experiment or simulation or reduced')
         end
            obj.data_type=val;
            notify(obj,'SourceChanged');
        end
        
        function set.model_type(obj,val)
         if ~(strcmpi(val,'HHS') ||... 
            strcmpi(val,'HHMS')||... 
            strcmpi(val,'Coupled_HHS')||... 
            strcmpi(val,'HHSIP')||... 
            strcmpi(val,'MHHMS')||...  %Multiplicative HHSIP
            strcmpi(val,'HHMSIP')||...  
            strcmpi(val,'HHSTM'))
            error('model_type must be one of the existing models')
         end
            obj.model_type=val;
            notify(obj,'SourceChanged');
        end
        
        function set.stim_type(obj,val)
         if ~(strcmpi(val,'periodical') ||... 
            strcmpi(val,'Poisson') ||... 
            strcmpi(val,'sinusoids') ||... 
            strcmpi(val,'1/f'))
            error('stim_type must be periodical or Poisson or 1/f or sinusoids ')
         end
            obj.stim_type=val;
            notify(obj,'SourceChanged');
        end
        
        function set.length_type(obj,val)
         if ~(strcmpi(val,'long') ||... 
            strcmpi(val,'short') )
            error('stim_type must be long or short')
         end
            obj.length_type=val;
            notify(obj,'SourceChanged');
        end
        
        function set.start_cut(obj,val)
         if ~((val<=1) && (val>=0))
             error('start_cut must between 0 and 1')
         end
            obj.start_cut=val;
            notify(obj,'SourceChanged');
        end
        
        function set.end_cut(obj,val)
         if  ~((val<=1) && (val>=0))
             error('end_cut must between 0 and 1')
         end
            obj.end_cut=val;
            notify(obj,'SourceChanged');
        end
 

        
    end
    
end

