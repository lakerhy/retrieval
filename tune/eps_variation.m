function[Tbv,Tbh] = eps_variation (freq,epsi1,epsii1,epsi2,epsii2)
  
% air:  theta1 (1,1)
% --------------- R1 
% snow:  theta2 (epsi1,epsii1)
%----------------- R2
% ice:  theta3 (epsi2,epsii2)
  step= 100
%   epsi2 = 3.2;
%   epsii2 = 0.2;
  c = 2.99793;
  % frq = dlmread ('Freq-Memls-in.txt');
  
  lamd=c/(10*freq);
  
  d = 0.2 ; % [m]
  T1 = 0;
  
T2 = 250;
T3 = 270;


  for n = 1:step
    theta1 = 55/180*pi;
    theta2 = asin(sin(theta1)/sqrt(epsi1(n)));
    theta3 = asin(sin(theta1)/sqrt(epsi2));
    n1 = 1;
    n2 = sqrt(epsi1(n));
    n3 = sqrt(epsi2);
    
    Rh1(n) = (n1*cos(theta1)-n2*cos(theta2))/(n1*cos(theta1)+n2* ...
                                              cos(theta2));
    Rh1(n) = Rh1(n)*Rh1(n);
    Rv1(n) = (n2*cos(theta1)-n1*cos(theta2))/(n2*cos(theta1)+n1* ...
                                             cos(theta2));
    Rv1(n) = Rv1(n)*Rv1(n);

    
    Rh2(n) = (n2*cos(theta2)-n3*cos(theta3))/(n2*cos(theta2)+n3* ...
                                              cos(theta3));
    Rh2(n) = Rh2(n)*Rh2(n);
    Rv2(n) = (n3*cos(theta2)-n2*cos(theta3))/(n3*cos(theta2)+n2* ...
                                             cos(theta3));
    Rv2(n) = Rv2(n)*Rv2(n);
    
    t =exp(-4*pi/lamd*(imag(sqrt(epsi1(n)+i*epsii1(n))))*d* ...
               abs(sec(theta2)))
    
    rsah(n) = Rh1(n) + (1-Rh1(n))^2*Rh2(n)*t^2/(1-Rh1(n)*Rh2(n));
    rsav(n) = Rv1(n) + (1-Rv1(n))^2*Rv2(n)*t^2/(1-Rv1(n)*Rv2(n));
    tsah(n)= t*(1-Rh1(n))*(1-Rh2(n))/(1-t^2*Rh1(n)*Rh2(n));
        tsav(n)= t*(1-Rv1(n))*(1-Rv2(n))/(1-t^2*Rv1(n)*Rv2(n));
        esah(n) = 1 - rsah(n) -tsah(n);
         esav(n) = 1 - rsav(n) -tsav(n);
         
         Tbh(n) = rsah(n)*T1+esah(n)*T2+ tsah(n)*T3;
         Tbv(n) = rsav(n)*T1+esav(n)*T2+ tsav(n)*T3;
  end
  size(Tbv)
