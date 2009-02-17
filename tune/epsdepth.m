function[Tbv,Tbh,rsah,rsav] = epsdepth (frequency,epsi1,epsii1,epsi2,epsii2)
  
% air:  theta1 (1,1)
% --------------- R1 
% snow:  theta2 (epsi1,epsii1)
%----------------- R2
% ice:  theta3 (epsi2,epsii2)
% e_{snow}: 
%  
%   epsi_step =30;
%   epsi1 = 1.4;
%   % epsi1 = linspace(1,3.2,epsi_step);
%   epsii1 = 0.01;
%   epsi2 = 3.15;
%   epsii2 = 0.1;
  c = 2.99793;
  
  d = linspace(0.01,0.2,100); % [m]
  T1 = 0;
  
  T2 = 250;
  T3 = 270;

  freq=frequency;
    lamd=c/(10*freq);


    for n = 1:100
      theta1 = 55/180*pi;
      theta2 = asin(sin(theta1)/sqrt(epsi1));
      theta3 = asin(sin(theta1)/sqrt(epsi2));
      n1 = 1;
      n2 = sqrt(epsi1);
      n3 = sqrt(epsi2);
      
      Rh1 = (n1*cos(theta1)-n2*cos(theta2))/(n1*cos(theta1)+n2* ...
                                                cos(theta2));
      Rh1 = Rh1*Rh1;
      Rv1 = (n2*cos(theta1)-n1*cos(theta2))/(n2*cos(theta1)+n1* ...
                                                cos(theta2));
      Rv1 = Rv1*Rv1;

      
      Rh2 = (n2*cos(theta2)-n3*cos(theta3))/(n2*cos(theta2)+n3* ...
                                                cos(theta3));
      Rh2 = Rh2*Rh2;
      Rv2 = (n3*cos(theta2)-n2*cos(theta3))/(n3*cos(theta2)+n2* ...
                                                cos(theta3));
      Rv2 = Rv2*Rv2;

      t(n) =exp(-4*pi/lamd*(imag(sqrt(epsi1+i*epsii1)))*d(n)* ...
                  abs(sec(theta2)));
      
      rsah(n) = Rh1 + (1-Rh1)^2*Rh2*t(n)^2/(1-Rh1*Rh2);
      rsav(n) = Rv1 + (1-Rv1)^2*Rv2*t(n)^2/(1-Rv1*Rv2);
      tsah(n)= t(n)*(1-Rh1)*(1-Rh2)/(1-t(n)^2*Rh1*Rh2);
      tsav(n)= t(n)*(1-Rv1)*(1-Rv2)/(1-t(n)^2*Rv1*Rv2);
      esah(n) = 1 - rsah(n) -tsah(n);
      esav(n) = 1 - rsav(n) -tsav(n);
      Tbh(n) = rsah(n)*T1+esah(n)*T2+ tsah(n)*T3;
      Tbv(n) = rsav(n)*T1+esav(n)*T2+ tsav(n)*T3;
  end
