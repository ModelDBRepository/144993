classdef Rates
    %RATES Summary of this class goes here
    
%     This objects from class holds together all the relevant (average)
%     rates a specific case (0,p or m).
    
    %   Detailed explanation goes here
    
    % the methods of this class generate  % A,b,B,D for the relevant SDE 
    % with the reduced state space 
    
    % Note: this version assumes uncoupled first order kinetics
    
    properties
        gamma=[];
        delta=[];

    end
    
    properties (Hidden)
        tau_AP=[];
    end
    
            
    
    methods
        %constuctor
        function obj=Rates(varargin)  %gamma,delta
            if size(varargin,2)==0
                obj.gamma=[];
                obj.delta=[];
                obj.tau_AP=[];
            elseif size(varargin,2)==2
                obj.gamma=varargin{1};
                obj.delta=varargin{2};
                obj.tau_AP=[];
            elseif size(varargin,2)==3
                obj.gamma=varargin{1};
                obj.delta=varargin{2};
                obj.tau_AP=varargin{3};
            else
                error('wrong number of input!');
            end
        end
        
        % overloaded functions
        function obj=plus(obj1,obj2)  % +  addition
            if ~(isa(obj1, 'Rates')&& isa(obj2, 'Rates'))
                error('both input must be "Rates" objects')
            end
            
            obj=Rates;
            prop_list=properties(obj);
            for ii=1:length(prop_list)
                prop=cell2mat(prop_list(ii));
                obj.(prop)=obj1.(prop)+obj2.(prop);    
            end 
        end
        
        function obj=minus(obj1,obj2)  % - substraction
            if ~(isa(obj1, 'Rates')&& isa(obj2, 'Rates'))
                error('both input must be "Rates" objects')
            end
            
            obj=Rates;
            prop_list=properties(obj);
            for ii=1:length(prop_list)
                prop=cell2mat(prop_list(ii));
                obj.(prop)=minus(obj1.(prop),obj2.(prop));    
            end 
        end
        
        function obj=mtimes(obj1,x)  % * multiplication
            obj=Rates;
            prop_list=properties(obj);
            for ii=1:length(prop_list)
                prop=cell2mat(prop_list(ii));
                obj.(prop)=obj1.(prop)*x;    
            end 
        end
        
        function obj=times(obj1,x)  % .* multiplication
            obj=Rates;
            prop_list=properties(obj);
            for ii=1:length(prop_list)
                prop=cell2mat(prop_list(ii));
                obj.(prop)=obj1.(prop)*x;    
            end 
        end
        
        function obj=subs(obj1,x,y)  % subs(S, old, new)
            obj=Rates;
            prop_list=properties(obj);
            for ii=1:length(prop_list)
                prop=cell2mat(prop_list(ii));
                obj.(prop)=subs(obj1.(prop),x,y);    
            end 
        end
        
        
        % Get relevant rates-based matrices
        function A=get_A(obj)
            A=diag(-obj.gamma-obj.delta);
        end
        
        function b=get_b(obj)            
            b=-obj.delta;
            if size(b,2)>1
                b=b.';    %make b a colmun vector
            end
        end        

        function D=get_D(obj,s,N) 
            D=diag(  (obj.gamma.*s+obj.delta.*(1-s) )./N);
        end
        
        function B=get_B(obj,s,N) 
            B=diag(  sqrt((obj.gamma.*s+obj.delta.*(1-s) )./N) );
        end        
        
        function x_inf=get_steadyState(obj) %steady state distribution in reduced state space
            A=get_A(obj);
            b=get_b(obj);
            x_inf=A\b;          
        end
        
    end
    
end

