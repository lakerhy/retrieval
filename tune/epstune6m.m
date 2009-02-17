function[] = epstune6m ()
 epsi1 = 1.4;
 epsii1 = 0.01;
 epsi_step =30;
 epsi2 = linspace(1.5,3.2,epsi_step);
 epsii2 = 0.02;
 epsi3 = 3.2;
 epsii3 = 0.2;
  c = 2.99793;
  % frq = dlmread ('Freq-Memls-in.txt');
  freq = 6.9;    
  lamd=c/(10*freq);

  
  d1 = 0.2 ; % [m]
  d2 = 0.05; % [m]
  T_sky = 0;
  
  T_snow = 250;
  T_slushy = 260;
  T_ice = 264;

  for n = 1:epsi_step
    theta4 = 55/180*pi;
    theta3 = asin(sin(theta4)/sqrt(epsi1));
    theta2 = asin(sin(theta4)/sqrt(epsi2(n)));
    theta1 = asin(sin(theta4)/sqrt(epsi3));
    n4 = 1 ;
    n3 = sqrt(epsi1);
    n2 = sqrt(epsi2(n));
    n1 = sqrt(epsi3);
    
    Rh3(n) = (n1*cos(theta4)-n2*cos(theta3))/(n1*cos(theta4)+n2* ...
                                              cos(theta3));
    Rh3(n) = Rh3(n)*Rh3(n);
    Rv3(n) = (n2*cos(theta4)-n1*cos(theta3))/(n2*cos(theta4)+n1* ...
                                              cos(theta3));
    Rv3(n) = Rv3(n)*Rv3(n);
    Rh2(n) = (n2*cos(theta3)-n3*cos(theta2))/(n2*cos(theta3)+n3* ...
                                              cos(theta2));
    Rh2(n) = Rh2(n)*Rh2(n);
    Rv2(n) = (n3*cos(theta3)-n2*cos(theta2))/(n3*cos(theta3)+n2* ...
                                              cos(theta2));
    Rv2(n) = Rv2(n)*Rv2(n);
    
    Rh1(n) = (n3*cos(theta2)-n4*cos(theta1))/(n3*cos(theta2)+n4* ...
                                              cos(theta1));
    Rh1(n) = Rh1(n)*Rh1(n);
    Rv1(n) = (n4*cos(theta2)-n3*cos(theta1))/(n4*cos(theta2)+n3* ...
                                              cos(theta1));
    Rv1(n) = Rv1(n)*Rv1(n);
    
    t2(n) =exp(-4*pi/lamd*(imag(sqrt(epsi1+i*epsii1)))*d1* ...
            abs(sec(theta3)));
    e2(n) = 1- t2(n);
    
    t1(n) =exp(-4*pi/lamd*(imag(sqrt(epsi2(n)+i*epsii2)))*d2* ...
            abs(sec(theta2)));
    e1(n) = 1-t1(n);
    
       rsah(n) = Rh2(n) + (1-Rh2(n))^2*Rh3(n)*t2(n)^2/(1-Rh2(n)*Rh3(n));
     rsav(n) = Rv2(n) + (1-Rv2(n))^2*Rv3(n)*t2(n)^2/(1-Rv2(n)*Rv3(n));
%     rsah(n) = Rh2(n) + (1-Rh2(n))^2*Rh3(n)/(1-Rh2(n)*Rh3(n));
%     rsav(n) = Rv2(n) + (1-Rv2(n))^2*Rv3(n)/(1-Rv2(n)*Rv3(n));
M1_h = 0;
M1_v = 0;
M2_h = [t1(n)*Rh2(n),0;0,t2(n)*Rh3(n)];
M2_v = [t1(n)*Rv2(n),0;0,t2(n)*Rv3(n)];
M3_h = [t1(n)*Rh1(n),0;0,t2(n)*Rh2(n)];
M3_v = [t1(n)*Rv1(n),0;0,t2(n)*Rv2(n)];
M4_h = 0;
M4_v = 0;
E_h = [e1(n)*T_slushy;e2(n)*T_snow];
E_v =  [e1(n)*T_slushy;e2(n)*T_snow];
F_h = [e1(n)*T_slushy+t1(n)*(1-Rh1(n))*T_ice;e2(n)*T_snow];
F_v = [e1(n)*T_slushy+t1(n)*(1-Rv1(n))*T_ice;e2(n)*T_snow];
I = eye(2);
M5_h = M3_h * (inv(I - M1_h) * M2_h) + M4_h;
M5_v = M3_v * (inv(I - M1_v) * M2_v) + M4_v;
D_h = inv(I - M5_h) * (M3_h * inv(I - M1_h) * E_h + F_h);
D_v = inv(I - M5_v) * (M3_v * inv(I - M1_v) * E_v + F_v);

Tbh(n) = (1-Rh1(n))*D_h(2);
Tbv(n) = (1-Rv1(n))*D_v(2);
    
  end
  
  
  figure
    hold on 
    grid on
    plot(epsi2,Tbv,'r');
    plot(epsi2,Tbh,'g');
  title('Upwelling brightness temperature (MEMLS)');
xlabel('permitivity of slushy layer');
ylabel('Tb [K]')


