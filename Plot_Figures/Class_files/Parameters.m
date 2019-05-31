classdef Parameters < hgsetget 
    %UNTITLED Summary of this class goes here
    
    % Model parameters (also used to select data)
    
    %   Detailed explanation goes here
   
    properties (SetObservable)
        f_in=20 %[Hz]
        N=1e6
        M=1 
        epsilon=0.2
        alpha=0
        dt=5e-3; %[ms]
        I0=7.9; %[microampere]
        model_type='HHS';
        
    end
    
    properties (Dependent = true, SetAccess = private,Hidden)
        T   
    end
    
    events
        ParametersChanged
    end
    
    methods
        
       function obj =  Parameters(f_in,I0,N,M,epsilon,alpha,dt,model_type)
           obj.f_in=f_in;
           obj.I0=I0;
           obj.N=N;
           obj.M=M;
           obj.epsilon=epsilon;
           obj.alpha=alpha;
           obj.dt=dt;
           obj.model_type=model_type;       
      end
       
       function T = get.T(obj)
        T=1/obj.f_in;
       end
    end % Modulus get method
      
end

