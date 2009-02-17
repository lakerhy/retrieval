function[] = memlscomp (type)

   frequency = [6.9,10.7,18.7,23.8,36.5]; 
  if type ==1 %FY
    cd ../MEMLS
    FY=load('FY.profile.1');
    T_snow=FY(2,2);
    T_ice=FY(1,2)
    W_ice=FY(1,3)
    roi_snow=FY(2,4);
    roi_ice=FY(1,4)
    pci_snow=FY(2,6);
    pci_ice=FY(1,6);
    sal=FY(1,7);
    di_snow=linspace(0.03,0.2,100);
    Tb_memls_V=zeros(100,5);
    Tb_memls_H=zeros(100,5);
    for n=1:100
      FY(2,5)=di_snow(n)*100;
   %    cd ../ying
%       Tb_FY_memls=icemain(55,0.5,0.5,FY,0,0,11);
   
      Tb_FY_memls=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', FY);
      Tb_memls_V(n,:) = Tb_FY_memls(:,1);
      Tb_memls_H(n,:) = Tb_FY_memls(:,2); 

    
      for i = 1:5
          cd ../tune
        freq=frequency(i);
        [epsi,epsii] = ro2epsd(roi_snow/1000,T_snow,freq);
        [epsi1(i),epsii1(i)] = mixmod(freq,T_snow,0,epsi,epsii);
        [epsi,epsii] = ro2epsd(roi_ice/1000,T_ice,freq);
        [epsi,epsii] = mixmod(freq,T_ice,W_ice,epsi,epsii);
        fy=1;
        [epsi2(i),epsii2(i)] = sie(fy,sal,T_ice,freq,epsi,epsii);
        
        [Tbv(:,i),Tbh(:,i),rsah,rsav] = epsdepth(freq,epsi1(i),epsii1(i),epsi2(i), ...
                                       epsii2(i));
        cd ../MEMLS
      end
    end
    
      figure

      grid on
      hold on
      title('FY:pol vs. snow depth using MEMLS(dashed line) and mine model(solid line) ')
      plot(di_snow,Tbv(:,1)-Tbh(:,1),'r');
      plot(di_snow,Tbv(:,2)-Tbh(:,2),'g');
      plot(di_snow,Tbv(:,3)-Tbh(:,3),'b');
      plot(di_snow,Tbv(:,4)-Tbh(:,4),'k');
      plot(di_snow,Tbv(:,5)-Tbh(:,5),'m');
      legend('6.9','10.7','18.7','23.8','36.5');



      plot(di_snow,abs(Tb_memls_V(:,1)-Tb_memls_H(:,1)),'r--');
      plot(di_snow,Tb_memls_V(:,2)-Tb_memls_H(:,2),'g--');
      plot(di_snow,Tb_memls_V(:,3)-Tb_memls_H(:,3),'b--');
      plot(di_snow,Tb_memls_V(:,4)-Tb_memls_H(:,4),'k--');
      plot(di_snow,Tb_memls_V(:,5)-Tb_memls_H(:,5),'m--');
      hold off
      ylabel('Tbv-Tbh');
      xlabel('snow depth [m]');
      cd ../tune
    
  end
  


    if type ==2 %MY
      
      cd ../MEMLS
      MY=load('MY.profile.1');
      T_snow=MY(2,2);
      T_ice=MY(1,2)
      W_ice=MY(1,3)
      roi_snow=MY(2,4);
      roi_ice=MY(1,4)
      pci_snow=MY(2,6);
      pci_ice=MY(1,6);
      sal=MY(1,7)
      di_snow=linspace(0.03,0.2,100);
      Tb_memls_V=zeros(100,5);
      Tb_memls_H=zeros(100,5);
      Tb_MY=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', MY);
T_MY_V=Tb_MY(:,1)
T_MY_H=Tb_MY(:,2)
      for n=1:100
        MY(2,5)=di_snow(n)*100;
        
             % cd ../ying
%      Tb_MY_memls=icemain(55,0.5,0.5,MY,0,0,11);
 Tb_MY_memls=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', MY);
%         Tb_MY_memls=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', MY,4);
        Tb_memls_V(n,:) = Tb_MY_memls(:,1);
        Tb_memls_H(n,:) = Tb_MY_memls(:,2); 
      end

      cd ../tune

      for l = 1:5
        freq=frequency(l);
        [epsi,epsii] = ro2epsd(roi_snow/1000,T_snow,freq);
        [epsi1(l),epsii1(l)] = mixmod(freq,T_snow,0,epsi,epsii);
        [epsi,epsii] = ro2epsd(roi_ice/1000,T_ice,freq);
        [epsi,epsii] = mixmod(freq,T_ice,W_ice,epsi,epsii);
        my=1;
        [epsi2(l),epsii2(l)] = mysie(my,roi_ice/1000,T_ice,sal,freq,epsi,epsii);
        
        [Tbv(:,l),Tbh(:,l)] = epsdepth(freq,epsi1(l),epsii1(l),epsi2(l),epsii2(l));
      end
      
      cd ../inversion
      Tb_amsr = amsr('area.33.2005.sort');
      Tb_amsr_V = Tb_amsr(1:2:10);
      Tb_amsr_H = Tb_amsr(2:2:10);
      cd ../tune
      
      figure
      subplot(2,1,1)
      hold on
      plot(frequency,T_MY_V,'r');
      plot(frequency,Tb_amsr_V,'g');
      legend('sim','amsr');
      ylabel('Tb_V');
      xlabel('frequency');
      title('AMSR and simulated');

            subplot(2,1,2)
      hold on
      plot(frequency,T_MY_H,'r');
      plot(frequency,Tb_amsr_H,'g');
      legend('sim','amsr');
      ylabel('Tb_H');
      xlabel('frequency');
      title('AMSR and simulated');

      figure
      subplot(3,2,1);
      grid on
      hold on
      title('MY:pol vs. snow depth using MEMLS(dashed line) and my model(solid line) ')
      plot(di_snow,Tbv(:,1)-Tbh(:,1));
      legend('6.9');
      plot(di_snow,Tb_memls_V(:,1)-Tb_memls_H(:,1),'--');
      ylabel('Tbv-Tbh [K]');
      xlabel('snow depth [m]');
      hold off

      subplot(3,2,2);
      grid on
      hold on
      title('MY:pol vs. snow depth using MEMLS(dashed line) and my model(solid line) ')
      plot(di_snow,Tbv(:,2)-Tbh(:,2));
      legend('10.7');
      plot(di_snow,Tb_memls_V(:,2)-Tb_memls_H(:,2),'--');
      ylabel('Tbv-Tbh [K]');
      xlabel('snow depth [m]');
      hold off

      subplot(3,2,3);
      grid on
      hold on
      title('MY:pol vs. snow depth using MEMLS(dashed line) and my model(solid line) ')
      plot(di_snow,Tbv(:,3)-Tbh(:,3));
      legend('18.7');
      plot(di_snow,Tb_memls_V(:,3)-Tb_memls_H(:,3),'--');
      ylabel('Tbv-Tbh [K]');
      xlabel('snow depth [m]');
      hold off

      subplot(3,2,4);
      grid on
      hold on
      title('MY:pol vs. snow depth using MEMLS(dashed line) and my model(solid line) ')
      plot(di_snow,Tbv(:,4)-Tbh(:,4));
      legend('23.8');
      plot(di_snow,Tb_memls_V(:,4)-Tb_memls_H(:,4),'--');
      ylabel('Tbv-Tbh [K]');
      xlabel('snow depth [m]');
      hold off

      subplot(3,2,5);
      grid on
      hold on
      title('MY:pol vs. snow depth using MEMLS(dashed line) and my model(solid line) ')
      plot(di_snow,Tbv(:,1)-Tbh(:,1));
      legend('36.5');
      plot(di_snow,Tb_memls_V(:,5)-Tb_memls_H(:,5),'--');
      ylabel('Tbv-Tbh [K]');
      xlabel('snow depth [m]');
      hold off
      hold off

    end

  end
