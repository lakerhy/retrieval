function[] = epstune ()
  epsi_snow_num = 30;
  epsi_snow = linspace(1.5,3.2,epsi_snow_num);
  epsii_snow = linspace(0,0.2,30);
  epsi_ice = 3.5;
  epsii_ice = 0.2;
  di = 0.2 ;% [m] snow depth
    T = 250;
  Tsky=271.35+2.7;
  Tgnd=0;

  frq = dlmread ('Freq-Memls-in.txt');

  cd Subroutine/Memls
  format long

  for n=1:length(frq)
    for i=1:epsi_snow_num
        epsii_snow = 0;
        epsi = [epsi_ice; epsi_snow(i)];
      epsii = [epsii_ice; epsii_snow];
      ns = sqrt(epsi);

      freq=frq(n);
      teta = 55;
      teta = (teta * pi) / 180;
      tei = [asin(sin(teta)./ns);teta];

      [sih,siv] = fresnelc(tei,[epsi;1])
      sih = sih.^2;
      siv = siv.^2;
      gai = abscoeff(epsi_snow(i),epsii_snow,T,freq)
      dei = pfadi(tei,di);
      
      % Calculate the reflectivity of the snow layer
      gamma = gai;
      t0i = exp(gamma * dei .* (-1));
      r0i = zeros(size(t0i));
      ri =0 ;
      ti = t0i;
      
      % H channel
      roa_h = ri + sih(1)*ti^2/(1-ri*sih(1));
      rob_h = ri + sih(2)*ti^2/(1-ri*sih(2));
      rsa_h(n,i) = sih(2)+(1-sih(2))^2*roa_h./(1-sih(2)*ri);
      rsb_h(n,i) = sih(1)+(1-sih(1))^2*rob_h./(1-sih(1)*ri);
      ts_h=ti*(1-sih(1))*(1-sih(2))/((1-ri*sih(1))*(1-ri*sih(2))-ti^2*sih(1)* ...
                                     sih(2));
      esa_h = 1- ts_h - rsa_h(n,i);
      Tb_h(i) = rsa_h(n,i)*Tsky+esa_h*T+ts_h*Tgnd;  
      % V channel
      roa_v = ri + siv(1)*ti^2/(1-ri*siv(1));
      rob_v = ri + siv(2)*ti^2/(1-ri*siv(2));
      rsa_v(n,i) = siv(2)+(1-siv(2))^2*roa_v/(1-siv(2)*ri);
      rsb_v(n,i) = siv(1)+(1-siv(1))^2*rob_v/(1-siv(1)*ri);
      ts_v=ti*(1-siv(1))*(1-siv(2))/((1-ri*siv(1))*(1-ri*siv(2))-ti^2*siv(1)* ...
                                     siv(2));
      esa_v = 1- ts_v - rsa_v(n,i);
      Tb_v(i) = rsa_v(n,i)*Tsky+esa_v*T+ts_v*Tgnd;  
      
      
    end
    TB_h(n,:) = Tb_h;
    TB_v(n,:) = Tb_v;
  end
  
  cd ../..
      
      
      figure 
  hold on
  plot(epsi_snow,rsa_v(5,:),'g');
  plot(epsi_snow,rsa_h(5,:),'r');
  
  
  figure
  hold on
  grid on
  plot(epsi_snow,TB_h(1,:),'r');
  plot(epsi_snow,TB_h(2,:),'g');
  plot(epsi_snow,TB_h(3,:),'b');
  plot(epsi_snow,TB_h(4,:),'m');
  plot(epsi_snow,TB_h(5,:),'k');
  legend('6H','11H','19H','24H','37H');
  title('TB_H vs. real part of the dielectric constant of snow');
  xlabel('epsi');
  ylabel('TB_H');
  hold off

        
  figure
  hold on
  grid on
  plot(epsi_snow,TB_v(1,:),'r');
  plot(epsi_snow,TB_v(2,:),'g');
  plot(epsi_snow,TB_v(3,:),'b');
  plot(epsi_snow,TB_v(4,:),'m');
  plot(epsi_snow,TB_v(5,:),'k');
  legend('6V','11V','19V','24V','37V')
   title('TB_V vs. real part of the dielectric constant of snow');
  xlabel('epsi');
  ylabel('TB_V');
 
  hold off

  figure
  hold on
  grid on
  plot(epsi_snow,TB_v(1,:)-TB_h(1,:),'r');
  plot(epsi_snow,TB_v(2,:)-TB_h(2,:),'g');
  plot(epsi_snow,TB_v(3,:)-TB_h(3,:),'b');
  plot(epsi_snow,TB_v(4,:)-TB_h(4,:),'m');
  plot(epsi_snow,TB_v(5,:)-TB_h(5,:),'k');
  legend('6V-H','11V-H','19V-H','24V-H','37V-H')
    title('Polarizaion vs. real part of the dielectric constant of snow');
  xlabel('epsi');
  ylabel('TB_V-TB_H');

  hold off
