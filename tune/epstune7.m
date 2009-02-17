function[] = epstune7 ()
% air: Z1 theta1 (1,1)
% --------------- R1 
% snow: T2,Z2 theta2 (epsi1,epsii1) ,L2
%----------------- R2
% slushy: T3, Z3 theta3 (epsi2,epsii2), L3
%----------------- R3
% ice: Z4 T4, theta4 (epsi3,epsii3)  
  
  epsi_step =30;
  
  epsi1 = 1.4;
  epsii1 = 0.01;

  
  epsi2 = linspace(1.5,3.2,epsi_step);
  epsii2 = 0.02;
  
  epsi3 = 3.2;
  epsii3 = 0.2;
  c = 2.99793;
  freq = 6.9;    
  lamd=c/(10*freq);

  d_snow = 0.2 ; % [m]
  d_slushy = 0.05 ; % [m]
  T2 = 250;
  T3 = 260;
  T4 = 264;
  for n = 1:epsi_step
    theta1 = 55/180*pi;
    theta2 = asin(sin(theta1)/sqrt(epsi1));
    theta3 = asin(sin(theta1)/sqrt(epsi2(n)));
    theta4 = asin(sin(theta1)/sqrt(epsi3));    
    
    Z1v = cos(theta1);
    Z1h = sec(theta1);
    Z2v = cos(theta2)/sqrt(epsi1);
    Z2h = sec(theta2)/sqrt(epsi1);
    Z3v = cos(theta3)/sqrt(epsi2(n));
    Z3h = sec(theta3)/sqrt(epsi2(n));
    Z4v = cos(theta4)/sqrt(epsi3);
    Z4h = sec(theta4)/sqrt(epsi3);

    R1v(n) = ((Z1v-Z2v)/(Z1v+Z2v))^2;
    R1h(n) = ((Z1h-Z2h)/(Z1h+Z2h))^2; 
    R2v(n) = ((Z2v-Z3v)/(Z2v+Z3v))^2;
    R2h(n)= ((Z2h-Z3h)/(Z2h+Z3h))^2;
    R3v(n) = ((Z3v-Z4v)/(Z3v+Z4v))^2;
    R3h(n)= ((Z3h-Z4h)/(Z3h+Z4h))^2;

    L2 = exp(4*pi/lamd*imag(sqrt(epsi1+i*epsii1))*d_snow*abs(sec(theta2)));
    L3 = exp(4*pi/lamd*imag(sqrt(epsi2(n)+i*epsii2))*d_slushy*abs(sec(theta3)))

    Sum_x1_h = 1/((1-R1h(n)*R2h(n)/L2^2));
    Sum_x1_v = 1/((1-R1v(n)*R2v(n)/L2^2))   ;
    Sum_x2_h = 1/((1-R2h(n)*R3h(n)/L3^2));
    Sum_x2_v = 1/((1-R2v(n)*R3v(n)/L3^2)); 
    
    Ts2 = T2*(1-1/L2);
    Ts3 = T3*(1-1/L3);    
    
    T2u_h(n) = (1-R1h(n))*Ts2*Sum_x1_h;
    T2u_v(n) = (1-R1v(n))*Ts2*Sum_x1_v;    
    T2D_h(n) = R2h(n)*(1-R1h(n))*Ts2*Sum_x1_h/L2;
    T2D_v(n) = R2v(n)*(1-R1v(n))*Ts2*Sum_x1_v/L2    
    
    T3u_h(n) = (1-R1h(n))* (1-R2h(n))*Ts3*Sum_x1_h*Sum_x2_h/L2;
    T3u_v(n) = (1-R1v(n))* (1-R2v(n))*Ts3*Sum_x1_v*Sum_x2_v/L2;
    T3D_h(n) = (1-R1h(n))* (1-R2h(n))*R3h(n)*Ts3*Sum_x1_h*Sum_x2_h/L2/L3;
    T3D_v(n) = (1-R1v(n))* (1-R2v(n))*R3v(n)*Ts3*Sum_x1_v*Sum_x2_v/L2/L3;
    
    T4u_h(n) = (1-R1h(n))* (1-R2h(n))*(1-R3h(n))*T4*Sum_x1_h*Sum_x2_h/L2/L3;
    T4u_v(n) = (1-R1v(n))* (1-R2v(n))*(1-R3v(n))*T4*Sum_x1_h*Sum_x2_h/L2/L3
    
    Tbh(n) = T2u_h(n)+T2D_h(n)+T3u_h(n)+T3D_h(n)+T4u_h(n);
    Tbv(n) = T2u_v(n)+T2D_v(n)+T3u_v(n)+T3D_v(n)+T4u_v(n);    
    
  end
%   figure
%   hold on
%   grid on
  plot(epsi2,Tbv,'r--'); 
  plot(epsi2,Tbh,'g--');
  legend('V','H');
  title('Upwelling brightness temperature');
  xlabel('permitivity of slushy layer');
ylabel('Tb [K]');
legend('V(MEMLS)','H(MEMLS)','V(Impedance)','H(Impedance)');

  hold off
