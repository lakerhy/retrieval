function[] = depthtune ()
  
% air:  theta1 (1,1)
% --------------- R1 
% snow:  theta2 (epsi1,epsii1)
%----------------- R2
% ice:  theta3 (epsi2,epsii2)
% e_{snow}: 
%  
  epsi_step =30;
  epsi1 = 1.53;
  % epsi1 = linspace(1,3.2,epsi_step);
  epsii1 = 0.073e-3;
  epsi2 = 3.4718;
  epsii2 = 0.0444;
  c = 2.99793;
  
  d = linspace(0,0.5,100); % [m
  T1 = 0;
  
  T2 = 250;
  T3 = 270;


  frequency = [6.9,10.7,18.7,23.8,36.5];
  for m=1:5
    freq=frequency(m);
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

%      t(n,i) =exp(-4*pi/lamd*(imag(sqrt(epsi1+i*epsii1)))*d(n)* ...
%                   abs(sec(theta2)));
      
%       rsah(n,i) = Rh1(n) + (1-Rh1(n))^2*Rh2(n)*t(n,i)^2/(1-Rh1(n)*Rh2(n));
%       rsav(n,i) = Rv1(n) + (1-Rv1(n))^2*Rv2(n)*t(n,i)^2/(1-Rv1(n)*Rv2(n));
%       tsah(n,i)= t(n,i)*(1-Rh1(n))*(1-Rh2(n))/(1-t(n,i)^2*Rh1(n)*Rh2(n));
%       tsav(n,i)= t(n,i)*(1-Rv1(n))*(1-Rv2(n))/(1-t(n,i)^2*Rv1(n)*Rv2(n));
%       esah(n,i) = 1 - rsah(n,i) -tsah(n,i);
%       esav(n,i) = 1 - rsav(n,i) -tsav(n,i);
      
%       Tbh(n,i) = rsah(n,i)*T1+esah(n,i)*T2+ tsah(n,i)*T3;
%       Tbv(n,i) = rsav(n,i)*T1+esav(n)*T2+ tsav(n,i)*T3;


      t(n,m) =exp(-4*pi/lamd*(imag(sqrt(epsi1+i*epsii1)))*d(n)* ...
                  abs(sec(theta2)));
      
      rsah(n,m) = Rh1 + (1-Rh1)^2*Rh2*t(n,m)^2/(1-Rh1*Rh2);
      rsav(n,m) = Rv1 + (1-Rv1)^2*Rv2*t(n,m)^2/(1-Rv1*Rv2);
      tsah(n,m)= t(n,m)*(1-Rh1)*(1-Rh2)/(1-t(n,m)^2*Rh1*Rh2);
      tsav(n,m)= t(n,m)*(1-Rv1)*(1-Rv2)/(1-t(n,m)^2*Rv1*Rv2);
      esah(n,m) = 1 - rsah(n,m) -tsah(n,m);
      esav(n,m) = 1 - rsav(n,m) -tsav(n,m);
      Tbh(n,m) = rsah(n,m)*T1+esah(n,m)*T2+ tsah(n,m)*T3;
      Tbv(n,m) = rsav(n,m)*T1+esav(n,m)*T2+ tsav(n,m)*T3;
      
%       t(n) =exp(-4*pi/lamd*(imag(sqrt(epsi1+i*epsii1)))*d(n)* ...
%                   abs(sec(theta2)));
      
%       rsah(n) = Rh1 + (1-Rh1)^2*Rh2*t(n)^2/(1-Rh1*Rh2);
%       rsav(n) = Rv1 + (1-Rv1)^2*Rv2*t(n)^2/(1-Rv1*Rv2);
%       tsah(n)= t(n)*(1-Rh1)*(1-Rh2)/(1-t(n)^2*Rh1*Rh2);
%       tsav(n)= t(n)*(1-Rv1)*(1-Rv2)/(1-t(n)^2*Rv1*Rv2);
%       esah(n) = 1 - rsah(n) -tsah(n);
%       esav(n) = 1 - rsav(n) -tsav(n);
%       Tbh(n) = rsah(n)*T1+esah(n)*T2+ tsah(n)*T3;
%       Tbv(n) = rsav(n)*T1+esav(n)*T2+ tsav(n)*T3;
        
end
  end

%   figure
%   plot(frequency,Tbv(25,:))
%   figure
%   grid on

%   hold on
%   plot(d,Tbv(:,1),'r');
%   plot(d,Tbh(:,1),'g');

%   title('Upwelling brightness temperature');
%   xlabel('depth of snow');
%   ylabel('TB');
%   %  legend('V(impedance)','H(impedance)','V(MEMLS)','H(MEMLS)');
%   hold off

cd ../MEMLS/Data
  buffer = load('area.33.2005.sort');
f1V = buffer(:,6); % 6
f1H = buffer(:,7);% 17
f2V = buffer(:,8);
f2H = buffer(:,9);
f3V = buffer(:,10); % 6
f3H = buffer(:,11);% 17
f4V = buffer(:,12);
f4H = buffer(:,13);
f5V = buffer(:,14); % 6
f5H = buffer(:,15);% 17

f1Vmean = mean(f1V);
f1Hmean = mean(f1H);
f2Vmean = mean(f2V);
f2Hmean = mean(f2H);
f3Vmean = mean(f3V);
f3Hmean = mean(f3H);
f4Vmean = mean(f4V);
f4Hmean = mean(f4H);
f5Vmean = mean(f5V);
f5Hmean = mean(f5H);


cd ../../tune

num = size(f1V);
f1pol = f1V-f1H;
f2pol = f2V-f2H;
f3pol = f3V-f3H;
f4pol = f4V-f4H;
f5pol = f5V-f5H;

pol1mean = mean(f1pol);
pol2mean = mean(f2pol);
pol3mean = mean(f3pol);
pol4mean = mean(f4pol);
pol5mean = mean(f5pol);

figure

  
  
subplot(2,1,1);  
  grid on
hold on
  plot(d,Tbv(:,1)-Tbh(:,1),'g');
  plot(d,Tbv(:,2)-Tbh(:,2),'b');
  plot(d,Tbv(:,3)-Tbh(:,3),'k');
  plot(d,Tbv(:,4)-Tbh(:,4),'m');
  plot(d,Tbv(:,5)-Tbh(:,5),'c');
  legend('6.9','11','19','24','37');
    title('Polariztion vs. depth of snow');
  xlabel('depth [m]');
  ylabel('TB_V-TB_H [K]');
hold off

subplot(2,1,2);
hold on
  grid on
plot(d,Tbv(:,1)-Tbh(:,1)-pol1mean,'g');
  plot(d,Tbv(:,2)-Tbh(:,2)-pol2mean,'b');
  plot(d,Tbv(:,3)-Tbh(:,3)-pol3mean,'k');
  plot(d,Tbv(:,4)-Tbh(:,4)-pol4mean,'m');
  plot(d,Tbv(:,5)-Tbh(:,5)-pol5mean,'c');
legend('6.9','11','19','24','37');
  
  title('Polariztion difference vs. depth of snow');
  xlabel('depth [m]');
  ylabel('Pol-Pol_{33} [K]');
 hold off

Tbv(:,1)-Tbh(:,1)
