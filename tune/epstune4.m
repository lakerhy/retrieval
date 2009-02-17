function[] = epstune4 ()
% air: Z1 theta1 (1,1)
% --------------- R1 
% snow: Z2 theta2 (epsi1,epsii1)
%----------------- R2
% ice: Z3 theta3 (epsi2,epsii2)
  epsi_step =30;
  epsi1 = linspace(1,3.2,epsi_step);
  epsii1 = 0.1;
  epsi2 = 3.2;
  epsii2 = 0.5;
  c = 2.99793;
  % frq = dlmread ('Freq-Memls-in.txt');
  freq = 37;    
  lamd=c/(10*freq);

  d = 0.2 ; % [m]
  T2 = 250;
  T3 = 270;
  for n = 1:epsi_step
    theta1 = 55/180*pi;
    theta2 = asin(sin(theta1)/sqrt(epsi1(n)));
    theta3 = asin(sin(theta1)/sqrt(epsi2));
    Z1v = cos(theta1);
    Z1h = sec(theta1);
    Z2v = cos(theta2)/sqrt(epsi1(n));
    Z2h = sec(theta2)/sqrt(epsi1(n));
    Z3v = cos(theta3)/sqrt(epsi2);
    Z3h = sec(theta3)/sqrt(epsi2);

    R1v(n) = ((Z1v-Z2v)/(Z1v+Z2v))^2;
    R1h(n) = ((Z1h-Z2h)/(Z1h+Z2h))^2; 
    R2v(n) = ((Z2v-Z3v)/(Z2v+Z3v))^2;
    R2h(n)= ((Z2h-Z3h)/(Z2h+Z3h))^2;
4*pi/lamd*imag(sqrt(epsi1(n)+i*epsii1))*d*abs(sec(theta2))
    L2 = exp(4*pi/lamd*imag(sqrt(epsi1(n)+i*epsii1))*d*abs(sec(theta2)));
    
    Tbh(n) = (1-R1h(n))*((1+R2h(n)/L2)*(1-1/L2)*T2+(1-R2h(n))*T3/L2)/(1-R1h(n)*R2h(n)/L2^2);
    Tbv(n) = (1-R1v(n))*((1+R2v(n)/L2)*(1-1/L2)*T2+(1-R2v(n))*T3/L2)/(1-R1v(n)*R2v(n)/L2^2);
    
  end
  
  
%     figure
%   hold on
%   grid on
%   plot(epsi1,R2v,'r');
%   plot(epsi1,R2h,'g');
%   hold off

  
  figure
  hold on
  grid on
  plot(epsi1,Tbv,'r--');
  plot(epsi1,Tbh,'g--');

hold off  
  
