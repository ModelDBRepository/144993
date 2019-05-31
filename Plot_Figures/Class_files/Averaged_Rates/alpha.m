function res = alpha( V, x)

phi_s=1/20;

    if x=='n'
        res= 0.01*(V+55)./(1-exp(-0.1*(V+55)));
    elseif x=='m'
        res= 0.1*(V+40)./(1-exp(-0.1.*(V+40)));
    elseif x=='h'
        res= 0.07.*exp(-(V+65)./20);
    elseif x=='s'
        res= phi_s*0.001.*exp(-(V+85)./30); % 0.001.*exp(-(V+85)./30); - changed from Fleidervish1996    
    else
        error('Input is Wrong!');
    end 
    
