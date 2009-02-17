function[] = phy2eps (type1,type2)
  if type1 ==1
    roi_snow = linspace(150,500,100); %kg/m3
    T_snow = 250; % K
    freq = [6.9,10.7,18.7,23.8,36.5];

    W  = 0;
    W1 = 0.01 ;
    W2 = 0.05

    for m =1:5
      for n = 1: 100
        [epsi(m,n), epsii(m,n)] = ro2epsd(roi_snow(n)/1000,T_snow,freq(m));
        [epsi_w1(m,n),epsii_w1(m,n)] = mixmod(freq(m),T_snow,W1,epsi(m,n),epsii(m,n))      ;
        [epsi_w2(m,n),epsii_w2(m,n)] = mixmod(freq(m),T_snow,W2,epsi(m, ...
                                                          n),epsii(m,n));
                
      end
      [Tbv_w1_FY(m,:),Tbh_w1_FY(m,:)]=eps_variation(freq(m),epsi_w1(m,:),epsii_w1(m,:),3.2,0.2);
  [Tbv_w2_FY(m,:),Tbh_w2_FY(m,:)]=eps_variation(freq(m),epsi_w2(m,:),epsii_w2(m,:),3.2,0.2);
  [Tbv_w1_MY(m,:),Tbh_w1_MY(m,:)]=eps_variation(freq(m),epsi_w1(m,:),epsii_w1(m,:),3.14,0.01);
  [Tbv_w2_MY(m,:),Tbh_w2_MY(m,:)]=eps_variation(freq(m),epsi_w2(m,:),epsii_w2(m,:),3.14,0.01);


    
    end


    m = find(abs(roi_snow - 150)< 1.75);
    n = find(abs(roi_snow - 200) < 1.75);
    k = find(abs(roi_snow - 300) < 1.75);
    l = find(abs(roi_snow - 400) < 1.75);
    h = find(abs(roi_snow - 450) < 1.75);

    figure
  %  subplot(2,1,1)
    hold on
    grid on
    plot(roi_snow,Tbv_w1_FY(1,:)-Tbh_w1_FY(1,:),'r');

    plot(roi_snow,Tbv_w1_FY(2,:)-Tbh_w1_FY(2,:),'g');
    plot(roi_snow,Tbv_w1_FY(3,:)-Tbh_w1_FY(3,:),'b');
    plot(roi_snow,Tbv_w1_FY(4,:)-Tbh_w1_FY(4,:),'k');
    plot(roi_snow,Tbv_w1_FY(5,:)-Tbh_w1_FY(5,:),'m');
    legend('6.9','10.7','18.7','23.8','36.5');
    plot(roi_snow,Tbv_w2_FY(1,:)-Tbh_w2_FY(1,:),'--r');
    plot(roi_snow,Tbv_w2_FY(2,:)-Tbh_w2_FY(2,:),'--g');
    plot(roi_snow,Tbv_w2_FY(3,:)-Tbh_w2_FY(3,:),'--b');
    plot(roi_snow,Tbv_w2_FY(4,:)-Tbh_w2_FY(4,:),'--k');
    plot(roi_snow,Tbv_w2_FY(5,:)-Tbh_w2_FY(5,:),'--m');

    title('On FY:Polariztion vs. density of snow with different water content W:0.01 solid,W2:0.05 dashed');
    xlabel('density [kg/m^3]');
    ylabel('Tb_v-Tb_h [K]');
     hold off

         figure
  %  subplot(2,1,1)
    hold on
    grid on
    plot(roi_snow,Tbv_w1_MY(1,:)-Tbh_w1_MY(1,:),'r');
    plot(roi_snow,Tbv_w1_MY(2,:)-Tbh_w1_MY(2,:),'g');
    plot(roi_snow,Tbv_w1_MY(3,:)-Tbh_w1_MY(3,:),'b');
    plot(roi_snow,Tbv_w1_MY(4,:)-Tbh_w1_MY(4,:),'k');
    plot(roi_snow,Tbv_w1_MY(5,:)-Tbh_w1_MY(5,:),'m');
    legend('6.9','10.7','18.7','23.8','36.5');
    plot(roi_snow,Tbv_w2_MY(1,:)-Tbh_w2_MY(1,:),'--r');
    plot(roi_snow,Tbv_w2_MY(2,:)-Tbh_w2_MY(2,:),'--g');
    plot(roi_snow,Tbv_w2_MY(3,:)-Tbh_w2_MY(3,:),'--b');
    plot(roi_snow,Tbv_w2_MY(4,:)-Tbh_w2_MY(4,:),'--k');
    plot(roi_snow,Tbv_w2_MY(5,:)-Tbh_w2_MY(5,:),'--m');

    title('On MY:Polariztion vs. density of snow with different water content W:0.01 solid,W2:0.05 dashed');
    xlabel('density [kg/m^3]');
    ylabel('Tb_v-Tb_h [K]');
     hold off

      figure
  %  subplot(2,1,1)
    hold on
    grid on
    plot(roi_snow,Tbv_w1_MY(1,:),'r');
    plot(roi_snow,Tbv_w1_MY(2,:),'g');
    plot(roi_snow,Tbv_w1_MY(3,:),'b');
    plot(roi_snow,Tbv_w1_MY(4,:),'k');
    plot(roi_snow,Tbv_w1_MY(5,:),'m');
    legend('6.9','10.7','18.7','23.8','36.5');
   %  plot(roi_snow,Tbv_w2_MY(1,:),'--r');
%     plot(roi_snow,Tbv_w2_MY(2,:),'--g');
%     plot(roi_snow,Tbv_w2_MY(3,:),'--b');
%     plot(roi_snow,Tbv_w2_MY(4,:),'--k');
%     plot(roi_snow,Tbv_w2_MY(5,:),'--m');

    title('On MY: Tbv vs. density of snow ');
    xlabel('density [kg/m^3]');
    ylabel('Tb_v [K]');
     hold off
     
     
      figure
  %  subplot(2,1,1)
    hold on
    grid on
    plot(roi_snow,Tbh_w1_MY(1,:),'r');
    plot(roi_snow,Tbh_w1_MY(2,:),'g');
    plot(roi_snow,Tbh_w1_MY(3,:),'b');
    plot(roi_snow,Tbh_w1_MY(4,:),'k');
    plot(roi_snow,Tbh_w1_MY(5,:),'m');
    legend('6.9','10.7','18.7','23.8','36.5');
   %  plot(roi_snow,Tbv_w2_MY(1,:),'--r');
%     plot(roi_snow,Tbv_w2_MY(2,:),'--g');
%     plot(roi_snow,Tbv_w2_MY(3,:),'--b');
%     plot(roi_snow,Tbv_w2_MY(4,:),'--k');
%     plot(roi_snow,Tbv_w2_MY(5,:),'--m');

    title('On MY: Tbh vs. density of snow ');
    xlabel('density [kg/m^3]');
    ylabel('Tb_v [K]');
     hold off

    figure
    hold on
    grid on
    plot(roi_snow,epsi(1,:),'r');
    plot(roi_snow,epsi(2,:),'g');
    plot(roi_snow,epsi(3,:),'b');
    plot(roi_snow,epsi(4,:),'k');
    plot(roi_snow,epsi(5,:),'m');
    legend('6.9G','10.7G','18.7G','23.8G','36.5G');
    title('Permitivity of dry snow vs density')
    xlabel('roi [kg/m^3]');
    ylabel('epsi');
    hold off

    figure 
    hold on
    grid on

    plot(freq,epsii(:,m),'r')
    plot(freq,epsii(:,n),'g')
    plot(freq,epsii(:,k),'b');
    plot(freq,epsii(:,l),'m');
    plot(freq,epsii(:,h),'k');
    legend('150','200','300','400','450');
    title('Im(e) vs frequency with various density ');
    xlabel('Frequency [GHz]');
    ylabel('epsi');
    hold off


    figure
    hold on
    grid on
    plot(roi_snow,epsii(1,:),'r');
    plot(roi_snow,epsii(2,:),'g');
    plot(roi_snow,epsii(3,:),'b');
    plot(roi_snow,epsii(4,:),'k');
    plot(roi_snow,epsii(5,:),'m');
    legend('6.9G','10.7G','18.7G','23.8G','36.5G');
    title('Im(e) of dry snow vs density')
    xlabel('roi [kg/m^3]');
    ylabel('epsii');
    hold off

    figure 
    hold on
    grid on
    plot(freq,epsii_w1(:,k),'g');
    plot(freq,epsii_w2(:,k),'b');
    legend('1%','5%');
    title('Imaginary part of permitivity vs frequency with various water content (rho=300kg/m^3)');
    xlabel('Frequency [GHz]');
    ylabel('epsii');
    hold off


    figure 
    hold on
    grid on
    plot(roi_snow,epsi(1,:),'r');
    plot(roi_snow,epsi_w1(1,:),'g')
    plot(roi_snow,epsi_w2(1,:),'b');
    legend('no water','1%','5%');
    title('Real part of permitivity vs density with various water content')
    xlabel('roi [kg/m^3]');
    ylabel('epsi');
    hold off

    
  end
  
  
  if type1 ==2 % ice
    
    
    T_ice = 260; % K
    freq = [6.9,10.7,18.7,23.8,36.5];
    W = linspace(0.001,0.01,51);
    w1 = find(abs(W - 0.002)< 0.00009)
    w2 = find(abs(W - 0.008) < 0.00009)
    w3 = find(abs(W - 0.01) < 0.00009)

    sal = linspace(0.5,7.5,51);
    s1 = find(abs(sal - 0.5)< 0.07);
    s2 = find(abs(sal - 1.5) < 0.07) 
    s3 = find(abs(sal - 2.5) < 0.07);
    s4 = find(abs(sal - 5)< 0.07);
    s5 = find(abs(sal - 6) < 0.07) 
    s6 = find(abs(sal - 7) < 0.07);

    roi_ice = linspace(700,1000,102); %kg/m3
    r1 = find(abs(roi_ice - 800)< 1.5);
    r2 = find(abs(roi_ice - 850.1) < 1.5);
    r3 = find(abs(roi_ice - 900) < 1.5);
    r4 = find(abs(roi_ice - 920) < 1.5);
    %--------------------------- Density
    %---------------------------------------------
    
    if type2 ==1
    for m =1:5
      for n = 1: 102
        [epsi_FY(m,n), epsii_FY(m,n)] = ro2epsd(roi_ice(n)/1000,T_ice,freq(m));
        [epsi_MY(m,n), epsii_MY(m,n)] = ro2epsd(roi_ice(n)/1000,T_ice,freq(m));
        [epsi_FY_w(m,n),epsii_FY_w(m,n)] = mixmod(freq(m),T_ice,W(w1),epsi_FY(m,n),epsii_FY(m,n));
        [epsi_MY_w(m,n),epsii_MY_w(m,n)] = mixmod(freq(m),T_ice,W(w1),epsi_MY(m,n),epsii_MY(m,n));
        [epsi_FY_w1(m,n),epsii_FY_w1(m,n)] = mixmod(freq(m),T_ice,W(w2),epsi_FY(m,n),epsii_FY(m,n));
        [epsi_FY_w2(m,n),epsii_FY_w2(m,n)] = mixmod(freq(m),T_ice,W(w3),epsi_FY(m,n),epsii_FY(m,n));           
        [epsi_MY_w1(m,n),epsii_MY_w1(m,n)] = mixmod(freq(m),T_ice,W(w2),epsi_MY(m,n),epsii_MY(m,n));
        [epsi_MY_w2(m,n),epsii_MY_w2(m,n)] = mixmod(freq(m),T_ice,W(w3), epsi_MY(m,n),epsii_MY(m,n));
        
        [epsi_FY_s(m,n),epsii_FY_s(m,n)] = sie(1,sal(s5),T_ice,freq(m),epsi_FY_w(m,n),epsii_FY_w(m,n));
        [epsi_MY_s(m,n),epsii_MY_s(m,n)] = mysie(1,roi_ice(n)/1000,T_ice,sal(s2),freq(m),epsi_MY_w(m,n),epsii_MY_w(m,n));
        [epsi_FY_s2(m,n),epsii_FY_s2(m,n)] = sie(1,sal(s4),T_ice,freq(m),epsi_FY_w(m,n),epsii_FY_w(m,n));
        [epsi_MY_s2(m,n),epsii_MY_s2(m,n)] = mysie(1,roi_ice(n)/1000,T_ice,sal(s1),freq(m),epsi_MY_w(m,n),epsii_MY_w(m,n));
        [epsi_FY_s3(m,n),epsii_FY_s3(m,n)] = sie(1,sal(s6),T_ice,freq(m),epsi_FY_w(m,n),epsii_FY_w(m,n));
        [epsi_MY_s3(m,n),epsii_MY_s3(m,n)] = mysie(1,roi_ice(n)/1000,T_ice,sal(s3),freq(m),epsi_MY_w(m,n),epsii_MY_w(m,n));
        
      end
    end
    
    figure 
    grid on
    hold on
    plot(roi_ice,epsi_FY_s(1,:),'r');
    plot(roi_ice,epsi_MY_s(1,:),'g');
    legend('FY','MY');
    xlabel('density [kg/m^3]');
    ylabel('epsi');
    title('Real part of permitivity of ice vs density W=0.002 sal_{FY}=6  sal_{MY}=0.5');
    hold off

    figure 
    grid on
    hold on
    plot(freq,epsii_FY_s(:,r1),'r');
    xlabel('frequency [GHz]');
    ylabel('epsi');
    title('Imaginary part of permitivity of FY ice vs frequency under various density   W=0.002  sal_{FY} = 6');

    hold off

    figure
    grid on
    hold on
    plot(freq,epsii_MY_s(:,r1),'r');
    plot(freq,epsii_MY_s(:,r2),'g');
    plot(freq,epsii_MY_s(:,r3),'b');
    plot(freq,epsii_MY_s(:,r4),'k');
    legend('800','850','900','920');
    xlabel('frequency [GHz]' );
    ylabel('epsi');
    title('Imaginary part of permitivity of MY ice vs frequency under various density W=0.002 sal_{MY} = 1.5');
    hold off
    
    %     figure
    %     grid on
    %     plot(roi_ice,epsi_FY_w1(1,:),'r');
    %     plot(roi_ice,epsi_FY_w1(5,:),'r--');
    %    %  plot(roi_ice,epsi_FY_w2(1,:),'r--');
    % %     plot(roi_ice,epsi_FY_w2(5,:),'r--');
    % %            %  plot(roi_ice,epsi_MY_w1(1,:),'b');
    % % %                 plot(roi_ice,epsi_MY_w2(5,:),'b--');
    % %     legend('FY W:0.02 6.7G','FY W:0.02 37G','FY W:0.02 6.7G','MY W:0.02 6.7G');
    %     xlabel('density [kg/m^3]' );
    %     ylabel('epsi');
    %     title('Imaginary part of permitivity of MY ice vs density under various water content sal_{FY} = 6 sal_{MY} = 1.5');
    %     hold off


    figure
    grid on
    hold on
    plot(roi_ice,epsi_FY_s2(1,:),'r');
    plot(roi_ice,epsi_FY_s2(5,:),'r--');
    plot(roi_ice,epsi_FY_s3(1,:),'b');
    plot(roi_ice,epsi_FY_s3(5,:),'b--');
    plot(roi_ice,epsi_MY_s2(1,:),'m');
    plot(roi_ice,epsi_MY_s2(5,:),'m--');
    plot(roi_ice,epsi_MY_s3(1,:),'k');
    plot(roi_ice,epsi_MY_s3(5,:),'k--');

    legend('FY 5 6.9G','FY 5 37G','FY 7 6.9G','FY 7 37G','MY 0.5 6.9G','MY 0.5 37G','MY 2.5 6.9G','MY 2.5 37G');
    xlabel('density [kg/m^3]' );
    ylabel('epsi');
    title('Real part of permitivity of MY ice vs density under various salinity W=0.002 ');
    hold off
    
    
    figure
    grid on
    hold on
    plot(roi_ice,epsii_FY_s2(1,:),'r');
    plot(roi_ice,epsii_FY_s2(5,:),'r--');
    plot(roi_ice,epsii_FY_s3(1,:),'b');
    plot(roi_ice,epsii_FY_s3(5,:),'b--');

    legend('5 6.9G','5 37G','7 6.9G','7 37G'),
    xlabel('density [kg/m^3]' );
    ylabel('epsi');
    title('Imaginary part of permitivity of FY ice vs density under various salinity W=0.002 ');
    hold off
    
    figure
    grid on
    hold on
    
    plot(roi_ice,epsii_MY_s2(1,:),'r');
    plot(roi_ice,epsii_MY_s2(5,:),'r--');
    plot(roi_ice,epsii_MY_s3(1,:),'b');
    plot(roi_ice,epsii_MY_s3(5,:),'b--');
    legend('0.5 6.9G','0.5 37G','2.5 6.9G','2.5 37G');
    xlabel('density [kg/m^3]' );
    ylabel('epsi');
    title('Imaginary part of permitivity of MY ice vs density under various salinity W=0.002 ');
    hold off
    end

    
 %   ---- Water Content ---FY:900,6.5--MY:920,1.5--------------------------------------
 if type2 ==2   
 for m =1:5
      for n = 1: 51
        [epsi_FY(m,n),epsii_FY(m,n)] = ro2epsd(roi_ice(r4)/1000,T_ice,freq(m));
        [epsi_MY(m,n),epsii_MY(m,n)] = ro2epsd(roi_ice(r3)/1000,T_ice,freq(m));
        [epsi_FY_w(m,n),epsii_FY_w(m,n)] = mixmod(freq(m),T_ice,W(n),epsi_FY(m,n),epsii_FY(m,n));
        [epsi_MY_w(m,n),epsii_MY_w(m,n)] = mixmod(freq(m),T_ice,W(n),epsi_MY(m,n),epsii_MY(m,n));
        
        [epsi_FY_s(m,n),epsii_FY_s(m,n)] = sie(1,sal(s5),T_ice,freq(m),epsi_FY_w(m,n),epsii_FY_w(m,n));
        [epsi_MY_s(m,n),epsii_MY_s(m,n)] = mysie(1,roi_ice(r3)/1000,T_ice,sal(s2),freq(m),epsi_MY_w(m,n),epsii_MY_w(m,n));
        [epsi_FY_s2(m,n),epsii_FY_s2(m,n)] = sie(1,sal(s4),T_ice,freq(m),epsi_FY_w(m,n),epsii_FY_w(m,n));
        [epsi_MY_s2(m,n),epsii_MY_s2(m,n)] = mysie(1,roi_ice(r3)/1000,T_ice,sal(s1),freq(m),epsi_MY_w(m,n),epsii_MY_w(m,n));
        [epsi_FY_s3(m,n),epsii_FY_s3(m,n)] = sie(1,sal(s6),T_ice,freq(m),epsi_FY_w(m,n),epsii_FY_w(m,n));
        [epsi_MY_s3(m,n),epsii_MY_s3(m,n)] = mysie(1,roi_ice(r3)/1000,T_ice,sal(s3),freq(m),epsi_MY_w(m,n),epsii_MY_w(m,n));
      end
    end

    figure 
    grid on
    hold on
    epsi_FY_s(1,:)
    plot(W,epsi_FY_s(1,:),'r');
    plot(W,epsi_MY_s(1,:),'b');
    %         plot(W,epsi_FY(5,:),'r--');
    %     plot(W,epsi_MY(5,:),'b--');

    legend('FY','MY');
    xlabel('water content []');
    ylabel('epsi');
    title('Real part of permitivity of ice vs water content rho _{FY}=920 rho_{MY}=900 sal_{FY}=6  sal_{MY}=1.5');
    hold off

    figure 
    grid on
    hold on
    plot(W,epsii_FY_s(1,:),'r');
    plot(W,epsii_MY_s(1,:),'b');
    plot(W,epsii_FY_s(5,:),'r--');
    plot(W,epsii_MY_s(5,:),'b--');

    legend('FY 6.9G','MY 6.9G','FY 37G','MY 37G');
    xlabel('water content []');
    ylabel('epsi');
    title('Real part of permitivity of ice vs water content rho _{FY}=920 rho_{MY}=900 sal_{FY}=6  sal_{MY}=1.5');
    hold off
    
    figure 
    grid on
    hold on
    plot(freq,epsii_FY(:,m),'r');
    xlabel('frequency [GHz]');
    ylabel('epsi');
    title('Imaginary part of permitivity of ice vs frequency under various density(W=0.002)');
    hold off
    end
  %  --------------------------- sal
  %  ---------------------------------------------
  if type2 ==3
    for m =1:5
      for n = 1: 51
        [epsi_FY(m,n), epsii_FY(m,n)] = ro2epsd(roi_ice(r4)/1000,T_ice,freq(m));
        [epsi_MY(m,n), epsii_MY(m,n)] = ro2epsd(roi_ice(r3)/1000,T_ice,freq(m));
        [epsi_FY(m,n),epsii_FY(m,n)] = mixmod(freq(m),T_ice,W(w1),epsi_FY(m,n),epsii_FY(m,n));
        [epsi_MY(m,n),epsii_MY(m,n)] = mixmod(freq(m),T_ice,W(w1),epsi_MY(m,n),epsii_MY(m,n));

        [epsi_FY(m,n),epsii_FY(m,n)] = sie(1,sal(n),T_ice,freq(m),epsi_FY(m,n),epsii_FY(m,n));
        [epsi_MY(m,n),epsii_MY(m,n)] = mysie(1,roi_ice(r3)/1000,T_ice,sal(n),freq(m),epsi_MY(m,n),epsii_MY(m,n));
        
      end
    end
    
    figure 
    grid on
    hold on
    plot(sal,epsi_FY(1,:),'r');
    plot(sal,epsi_FY(5,:),'r--');

    plot(sal,epsi_MY(1,:),'b');
    plot(sal,epsi_MY(5,:),'b--');
    legend('FY 6.9G','FY 37G','MY 6.9G','MY 37G');
    xlabel('Salility [ppt]');
    ylabel('epsi');
    title('Real part of permitivity of ice vs sal. W=0.002 roi_{FY}=920 roi_{MY}=900');
    hold off
    
    
    figure 
    grid on
    hold on
    plot(sal,epsii_FY(1,:),'r');
    plot(sal,epsii_FY(5,:),'r--');

    plot(sal,epsii_MY(1,:),'b');
    plot(sal,epsii_MY(5,:),'b--');
    legend('FY 6.9G','FY 37G','MY 6.9G','MY 37G');
    xlabel('Salility [ppt]');
    ylabel('epsii');
    title('Imaginary part of permitivity of ice vs sal. W=0.002 roi_{FY}=920 roi_{MY}=900');
    hold off
    figure
    grid on
    hold on 
    plot(freq,epsii_FY(:,s5),'r');
    plot(freq,epsii_MY(:,s2),'b');
    plot(freq,epsii_FY(:,s6),'r--');
    plot(freq,epsii_MY(:,s3),'b--');

    legend('FY 6','MY 2.5','FY 7','MY 2.5 ');
    xlabel('Frequency [GHz]');
    ylabel('epsii');

    title('Imaginary part of permitivity of ice vs frequency W=0.002 roi_{FY}=920 roi_{MY}=900');
    hold off
  end
  end
  
