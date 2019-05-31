function res = beta( V, x )

phi_s=1/20; %s1 is slower, to get relaxation as in Gal 2010, and 
amp=0.0034*3;sigma=1/0.3; % makes inactivation faster to make 

if x=='n'
        res= 0.125.*exp(-(V+65)./80);
    elseif x=='m'
        res= 4.*exp(-(V+65)./18);
    elseif x=='h'
        res= 1/(exp(-0.1*(V+35))+1);
    elseif x=='s'
        res= phi_s*amp./(exp(-(V+17))/sigma+1); %- orignally from Fleidervish1996    
    else
        error('Input is Wrong!');
end 
    
