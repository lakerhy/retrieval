function[] = epstune6 ()
  
% air:  theta1 (1,1)
% --------------- R1 
% snow:  theta2 (epsi1,epsii1) T1 t1 e1 r1
%                        Rsa  
%----------------- R2           
% slushy:  theta3 (epsi2,epsii3) T2 t2 e2 r2*T2
%----------------- R3  
%                        Rsb     
% ice:  theta4 (epsi3,epsii3) T3
  epsi1 = 1.4;
  epsii1 = 0.01;
  epsi_step =30;
    

   epsi2 = linspace(1.5,3.2,epsi_step);
   epsii2 = 0.02;
 
  epsi3 = 3.2;
  epsii3 = 0.2;
  c = 2.99793;
  % frq = dlmread ('Freq-Memls-in.txt');
  freq = 37;    
  lamd=c/(10*freq);

  
  d1 = 0.2 ; % [m]
  d2 = 0.05; % [m]
  T_sky = 0;
  
  T1 = 250;
  T2 = 260;
  T3 = 264;


  
  
  for n = 1:epsi_step
     theta1 = 55/180*pi;
    theta2 = asin(sin(theta1)/sqrt(epsi1));
    theta3 = asin(sin(theta1)/sqrt(epsi2(n)));
    theta4 = asin(sin(theta1)/sqrt(epsi3));
    n1 = 1 ;
    n2 = sqrt(epsi1);
    n3 = sqrt(epsi2(n));
    n4 = sqrt(epsi3);
    
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
    
    Rh3(n) = (n3*cos(theta3)-n4*cos(theta4))/(n3*cos(theta3)+n4* ...
                                              cos(theta4));
    Rh3(n) = Rh3(n)*Rh3(n);
    Rv3(n) = (n4*cos(theta3)-n3*cos(theta4))/(n4*cos(theta3)+n3* ...
                                              cos(theta4));
    Rv3(n) = Rv3(n)*Rv3(n);
    
    t1(n) =exp(-4*pi/lamd*(imag(sqrt(epsi1+i*epsii1)))*d1* ...
            abs(sec(theta2)));
    
    
    t2(n) =exp(-4*pi/lamd*(imag(sqrt(epsi2(n)+i*epsii2)))*d2* ...
            abs(sec(theta3)));
    
       rsah(n) = Rh2(n) + (1-Rh2(n))^2*Rh3(n)*t2(n)^2/(1-Rh2(n)*Rh3(n));
     rsav(n) = Rv2(n) + (1-Rv2(n))^2*Rv3(n)*t2(n)^2/(1-Rv2(n)*Rv3(n));
%     rsah(n) = Rh2(n) + (1-Rh2(n))^2*Rh3(n)/(1-Rh2(n)*Rh3(n));
%     rsav(n) = Rv2(n) + (1-Rv2(n))^2*Rv3(n)/(1-Rv2(n)*Rv3(n));
    
    tsah(n)= t2(n)*(1-Rh2(n))*(1-Rh3(n))/(1-t2(n)^2*Rh2(n)*Rh3(n));
    tsav(n)= t2(n)*(1-Rv2(n))*(1-Rv3(n))/(1-t2(n)^2*Rv2(n)*Rv3(n));
    
    esah(n) = 1 - rsah(n) -tsah(n);
    esav(n) = 1 - rsav(n) -tsav(n);
    
    th_snow(n)= t1(n)*(1-Rh1(n))*(1-rsah(n))/(1-t1(n)^2*Rh1(n)*rsah(n));
    tv_snow(n)= t1(n)*(1-Rv1(n))*(1-rsav(n))/(1-t1(n)^2*Rv1(n)*rsav(n));
    
    
    rsh_snow(n) = Rh1(n) + (1-Rh1(n))^2*rsah(n)*t1(n)^2/(1-Rh1(n)*rsah(n));
    rsv_snow(n) = Rv1(n) + (1-Rv1(n))^2*rsav(n)*t1(n)^2/(1-Rv1(n)*rsav(n));
    
    esh_snow(n) = 1-th_snow(n)-rsh_snow(n);
    esv_snow(n) = 1-tv_snow(n)-rsv_snow(n);
    T1h_slushy(n) = T_sky*th_snow(n);
    T1v_slushy(n) = T_sky*tv_snow(n);
    
    Tbh_slushy(n) = rsah(n)*T1h_slushy(n)+esah(n)*T2+ tsah(n)*T3;
    Tbv_slushy(n) = rsav(n)*T1v_slushy(n)+esav(n)*T2+ tsav(n)*T3;
    
    Tbv(n) = T_sky*rsv_snow(n)+ (1+tv_snow(n)*rsav(n))*esv_snow(n)*T1+tv_snow(n)*Tbv_slushy(n);
    Tbh(n) = T_sky*rsh_snow(n)+ (1+th_snow(n)*rsah(n))*esh_snow(n)*T1+th_snow(n)*Tbh_slushy(n);
  end
  
  
%   figure
%   hold on
%   grid on
%   plot(epsi2,rsv_snow,'r');
%   plot(epsi2,rsh_snow,'g');
%   title('Total reflectvity at air-snow interface vs. dielectric constant');
%   xlabel('Real part of dielectric constant of slushy layer');
%   ylabel('Total reflectvity');
%   legend('v','h');
%   hold off
  
  
%   figure
%   hold on
%   grid on
%   plot(epsi2,rsav,'r');
%   plot(epsi2,rsah,'g');
%   title('Total reflectvity at snow-slushy interface vs. dielectric constant');
%   xlabel('Real part of dielectric constant of slushy layer');
%   ylabel('Total reflectvity');
%   legend('v','h');
%   hold off


%   figure
%   hold on
%   grid on
%   plot(epsi2,Tbv_slushy,'r');
%   plot(epsi2,Tbh_slushy,'g');
%   title('Upwelling brightness temperature at snow-slushy interface vs. dielectric constant');
%   xlabel('Real part of dielectric constant of slushy layer');
%   ylabel('TB');
%   legend('v','h');
%   hold off

  
 
  plot(epsi2,Tbv,'r--');
  plot(epsi2,Tbh,'g--');
  title('Upwelling brightness temperature at air-snow interface vs. dielectric constant(MEMLS)');
  xlabel('permitivity of slushy layer');
  ylabel('TB');
  legend('V(impedance)','H(impedance)');
  hold off
%   figure
%   hold on 
%   grid on
%   plot(epsi2,tsah,'r');
%   plot(epsi2,tsav,'g');
%   legend('H','V')
%   hold off

  
    
%   figure
%   hold on 
%   grid on
%   plot(epsi2,th_snow,'r');
%   plot(epsi2,tv_snow,'g');
%   legend('H','V')
%   hold off
