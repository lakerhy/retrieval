function[] = epstune5 ()
  
% air:  theta1 (1,1)
% --------------- R1 
% snow:  theta2 (epsi1,epsii1)
%----------------- R2
% ice:  theta3 (epsi2,epsii2)
  epsi_step =30;
  epsi1 = linspace(1,3.2,epsi_step);
  epsii1 = 0.05;
  epsi2 = 3.2;
  epsii2 = 0.5;
  c = 2.99793;
  % frq = dlmread ('Freq-Memls-in.txt');
  freq = 37;    
  lamd=c/(10*freq);

  
  d = 0.2 ; % [m]
  T1 = 0;
  
T2 = 250;
T3 = 270;


  for n = 1:epsi_step
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
    
    t =exp(-4*pi/lamd*(imag(sqrt(epsi1(n)+i*epsii1)))*d* ...
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
  
  
%     figure
%   hold on
%   grid on
%   plot(epsi1,rsav,'r');
%   plot(epsi1,rsah,'g');
%   title('Total reflectvity at air-snow interface vs. dielectric constant');
%   xlabel('Real part of dielectric constant of snow');
%   ylabel('Total reflectvity');
%   legend('v','h');
%   hold off


%   figure
%   hold on
%   grid on
%   plot(epsi1,Rv2,'r');
%   plot(epsi1,Rh2,'g');
%   title('M. Reflectivity');
%   hold off
figure
hold on
plot(epsi1,Tbv,'r');
  plot(epsi1,Tbh,'g');
  title('Upwelling brightness temperature');
  xlabel('permitivity of snow');
  ylabel('TB');
  legend('V(impedance)','H(impedance)','V(MEMLS)','H(MEMLS)');
  hold off

  
    
%     figure
%   hold on
%   grid on
%   plot(epsi1,tsav,'r');
%   plot(epsi1,tsah,'g');
% %   title('Total reflectvity at air-snow interface vs. dielectric constant');
% %   xlabel('Real part of dielectric constant of snow');
% %   ylabel('Total reflectvity');
%   legend('v','h');
%   hold off

