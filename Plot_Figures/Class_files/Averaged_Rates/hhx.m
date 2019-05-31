%Classic HH parameters + s1 variable parameters are based on Lundstrom2008 (supporting information+paper) 
% The s variable parameters are orignally from Fleidervish1996

function dy = hhx(y,I)

    global A dt phi_HH
    
    
    VNa=50; VK=-77;VL=-54; %[mV]
    gNa=120;gK=36;gL=0.3; %[mS/cm^2]
    Cm=1/phi_HH; %[microFarad/cm^2]
      
      V=y(1);
      beta(1)= phi_HH*0.125.*exp(-(V+65)./80); %n
      beta(2)= phi_HH*4.*exp(-(V+65)./18); %m
      beta(3)= phi_HH*1/(exp(-0.1*(V+35))+1); %h

      alpha(1)=phi_HH*0.01*(V+55)./(1-exp(-0.1*(V+55)));%n
      alpha(2)=phi_HH*0.1*(V+40)./(1-exp(-0.1.*(V+40)));%m
      alpha(3)=phi_HH*0.07.*exp(-(V+65)./20); %h


    dy2dt(1)=(gNa.*(y(2).^3).*y(4).*A.*(VNa-y(1))+gK.*(y(3).^4).*(VK-y(1))+gL.*(VL-y(1))+I)./Cm;% dVm/dt
    dy2dt(2)=(1-y(2)).*alpha(2)-y(2).*beta(2);  % dm/dt
    dy2dt(3)=(1-y(3)).*alpha(1)-y(3).*beta(1); % dn/dt
    dy2dt(4)=(1-y(4)).*alpha(3)-y(4).*beta(3); % dh/dt
        
    dy=dt*dy2dt;