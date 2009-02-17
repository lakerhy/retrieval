
function[] = sens(type,partype)

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
    di_snow=FY(2,5);
     if partype==6 % T_s
      T_snow=linspace(235,270,100);
      x_axis =T_snow
      Tb_memls_V=zeros(100,5);
      Tb_memls_H=zeros(100,5);
      for n=1:100
        FY(2,2)=T_snow(n);
        Tb_FY_memls=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', FY);
        Tb_memls_V(n,:) = Tb_FY_memls(:,1);
        Tb_memls_H(n,:) = Tb_FY_memls(:,2); 
        delta_T_snow=(270-235)/1000*100; %[cm]
        FY(2,2)=T_snow(n)+delta_T_snow;
        Tb_FY_memls_delta=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', FY);
        Tb_memls_V_delta(n,:) = Tb_FY_memls_delta(:,1);
        Tb_memls_H_delta(n,:) = Tb_FY_memls_delta(:,2); 
        Dir_V(n,:)= (Tb_memls_V_delta(n,:)-Tb_memls_V(n,:))/delta_T_snow;%k/cm
        Dir_H(n,:)= (Tb_memls_H_delta(n,:)-Tb_memls_H(n,:))/delta_T_snow; ...
        %k/cm
        pol_delta(n,:)=Tb_memls_V_delta(n,:)-Tb_memls_H_delta(n,:);
        pol(n,:)=Tb_memls_V(n,:)-Tb_memls_H(n,:);
        Dir_pol(n,:)= (pol_delta(n,:)-pol(n,:))/delta_T_snow; 
        par = ' snow temperature';
        unit =  ' [K/K]';
        xaxis = '  snow temperature[K]';
      end
      cd ../inversion
    end
    
    
    if partype==1 % snow depth
      step=200
      di_snow=linspace(0.03,0.5,step);
      x_axis =di_snow
      Tb_memls_V=zeros(100,5);
      Tb_memls_H=zeros(100,5);
      for n=1:step
        FY(2,5)=di_snow(n)*100;
        Tb_FY_memls=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', FY);
        Tb_memls_V(n,:) = Tb_FY_memls(:,1);
        Tb_memls_H(n,:) = Tb_FY_memls(:,2); 
        delta_snow=(0.5-0.03)/2000*100; %[cm]
        FY(2,5)=di_snow(n)*100+delta_snow;
        Tb_FY_memls_delta=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', FY);
        Tb_memls_V_delta(n,:) = Tb_FY_memls_delta(:,1);
        Tb_memls_H_delta(n,:) = Tb_FY_memls_delta(:,2); 
        Dir_V(n,:)= (Tb_memls_V_delta(n,:)-Tb_memls_V(n,:))/delta_snow;%k/cm
        Dir_H(n,:)= (Tb_memls_H_delta(n,:)-Tb_memls_H(n,:))/delta_snow; ...
        %k/cm
        pol_delta(n,:)=Tb_memls_V_delta(n,:)-Tb_memls_H_delta(n,:);
        pol(n,:)=Tb_memls_V(n,:)-Tb_memls_H(n,:);
        Dir_pol(n,:)= (pol_delta(n,:)-pol(n,:))/delta_snow; 
        par = ' snow depth';
        unit =  ' [K/cm]';
        xaxis = ' snow depth [m]';
      end
      cd ../inversion
    end

    if partype==2 % snow density
      step=500;
      roi_snow=linspace(200,450,step);
      x_axis = roi_snow;
      for n=1:step
        FY(2,4)=roi_snow(n);
        
        Tb_FY_memls=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', FY);
        Tb_memls_V(n,:) = Tb_FY_memls(:,1);
        Tb_memls_H(n,:) = Tb_FY_memls(:,2); 
        delta_roi_snow=(450-200)/5000; 
        FY(2,4)=roi_snow(n)+delta_roi_snow;
        Tb_FY_memls_delta=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', FY);
        Tb_memls_V_delta(n,:) = Tb_FY_memls_delta(:,1);
        Tb_memls_H_delta(n,:) = Tb_FY_memls_delta(:,2); 
        Dir_V(n,:)= (Tb_memls_V_delta(n,:)-Tb_memls_V(n,:))/delta_roi_snow;%k/cm
        Dir_H(n,:)= (Tb_memls_H_delta(n,:)-Tb_memls_H(n,:))/delta_roi_snow; ...
        %k/cm
        pol_delta(n,:)=Tb_memls_V_delta(n,:)-Tb_memls_H_delta(n,:);
        pol(n,:)=Tb_memls_V(n,:)-Tb_memls_H(n,:);
        Dir_pol(n,:)= (pol_delta(n,:)-pol(n,:))/delta_roi_snow; 
        par=' snow density';
        unit = 'K/kg/m^3';
        xaxis=' snow density [kg/m^3]';
      end
      cd ../inversion
    end

    if partype==3 % grain size
      pci_snow=linspace(0.001,0.3,100);
      x_axis = pci_snow;
      for n=1:100
        FY(2,6)=pci_snow(n);
        Tb_FY_memls=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', FY);
        Tb_memls_V(n,:) = Tb_FY_memls(:,1);
        Tb_memls_H(n,:) = Tb_FY_memls(:,2); 
        delta_pci_snow=(0.3-0.001)/2500; 
        FY(2,6)=pci_snow(n)+delta_pci_snow;
        Tb_FY_memls_delta=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', FY);
        Tb_memls_V_delta(n,:) = Tb_FY_memls_delta(:,1);
        Tb_memls_H_delta(n,:) = Tb_FY_memls_delta(:,2); 
        Dir_V(n,:)= (Tb_memls_V_delta(n,:)-Tb_memls_V(n,:))/delta_pci_snow;%k/cm
        Dir_H(n,:)= (Tb_memls_H_delta(n,:)-Tb_memls_H(n,:))/delta_pci_snow; ...
        %k/cm
        pol_delta(n,:)=Tb_memls_V_delta(n,:)-Tb_memls_H_delta(n,:);
        pol(n,:)=Tb_memls_V(n,:)-Tb_memls_H(n,:);
        Dir_pol(n,:)= (pol_delta(n,:)-pol(n,:))/delta_pci_snow; 
        par=' grain size of snow';
        unit = 'K/mm';
        xaxis=' grain size of snow [mm]';
      end
      cd ../inversion
    end
    if partype==4 % salinity
      sal=linspace(0.05,8,100);
      x_axis =sal
      Tb_memls_V=zeros(100,5);
      Tb_memls_H=zeros(100,5);
      for n=1:100
        FY(1,7)=sal(n);
        Tb_FY_memls=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', FY);
        Tb_memls_V(n,:) = Tb_FY_memls(:,1);
        Tb_memls_H(n,:) = Tb_FY_memls(:,2); 
        delta_sal=(8-0.05)/1000; %[cm]
        FY(1,7)=sal(n)+delta_sal;
        Tb_FY_memls_delta=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', FY);
        Tb_memls_V_delta(n,:) = Tb_FY_memls_delta(:,1);
        Tb_memls_H_delta(n,:) = Tb_FY_memls_delta(:,2); 
        Dir_V(n,:)= (Tb_memls_V_delta(n,:)-Tb_memls_V(n,:))/delta_sal;%k/cm
        Dir_H(n,:)= (Tb_memls_H_delta(n,:)-Tb_memls_H(n,:))/delta_sal; ...
        %k/cm
        pol_delta(n,:)=Tb_memls_V_delta(n,:)-Tb_memls_H_delta(n,:);
        pol(n,:)=Tb_memls_V(n,:)-Tb_memls_H(n,:);
        Dir_pol(n,:)= (pol_delta(n,:)-pol(n,:))/delta_sal; 
        par = 'salinity';
        unit =  '?';
        xaxis = 'salinity';
      end
      cd ../inversion
    end
    
    if partype==5 %ice temperature 
      T_ice=linspace(250,275,100);
      x_axis =T_ice
      Tb_memls_V=zeros(100,5);
      Tb_memls_H=zeros(100,5);
      for n=1:100
        FY(1,2)=T_ice(n);
        Tb_FY_memls=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', FY);
        Tb_memls_V(n,:) = Tb_FY_memls(:,1);
        Tb_memls_H(n,:) = Tb_FY_memls(:,2); 
        delta_T_ice=(275-250)/500; 
        FY(1,2)=T_ice(n)+delta_T_ice;
        Tb_FY_memls_delta=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', FY);
        Tb_memls_V_delta(n,:) = Tb_FY_memls_delta(:,1);
        Tb_memls_H_delta(n,:) = Tb_FY_memls_delta(:,2); 
        Dir_V(n,:)= (Tb_memls_V_delta(n,:)-Tb_memls_V(n,:))/delta_T_ice;%k/cm
        Dir_H(n,:)= (Tb_memls_H_delta(n,:)-Tb_memls_H(n,:))/delta_T_ice; ...
        %k/cm
        pol_delta(n,:)=Tb_memls_V_delta(n,:)-Tb_memls_H_delta(n,:);
        pol(n,:)=Tb_memls_V(n,:)-Tb_memls_H(n,:);
        Dir_pol(n,:)= (pol_delta(n,:)-pol(n,:))/delta_T_ice; 
        par = 'T_ice';
        unit =  'K';
        xaxis = 'T_ice';
      end
      cd ../inversion
    end

    
    figure
    subplot(1,2,1)
    grid on
    hold on
    title (['Sensitivity of Tb_V according to ',par,'  (FY)']);
    plot(x_axis,Dir_V(:,1),'r--');
    plot(x_axis,Dir_V(:,2),'g--');
    plot(x_axis,Dir_V(:,3),'b--');
    plot(x_axis,Dir_V(:,4),'k--');
    plot(x_axis,Dir_V(:,5),'m--');
    hold off
    ylabel(['sensitivty of Tb_v', unit]);
    xlabel(xaxis);
    legend('6.9','10.7','18.7','23.8','36.5' );
    
    subplot(1,2,2)
    grid on
    hold on
    title (['Sensitivity of Tb_H according to',par,'  (FY)']);
    plot(x_axis,Dir_H(:,1),'r--');
    plot(x_axis,Dir_H(:,2),'g--');
    plot(x_axis,Dir_H(:,3),'b--');
    plot(x_axis,Dir_H(:,4),'k--');
    plot(x_axis,Dir_H(:,5),'m--');
    hold off
    ylabel(['sensitivty of Tb_h', unit]);
    xlabel(xaxis);
    legend('6.9','10.7','18.7','23.8','36.5' );


  figure
  grid on
  hold on
  
  title (['Sensitivity of polarization according to ',par, ' (FY)']);
  plot(x_axis,Dir_pol(:,1),'r--');
  plot(x_axis,Dir_pol(:,2),'g--');
  plot(x_axis,Dir_pol(:,3),'b--');
  plot(x_axis,Dir_pol(:,4),'k--');
  plot(x_axis,Dir_pol(:,5),'m--');
  hold off
  ylabel(['sensitivty of polarization ',unit]);
  xlabel(xaxis);
  legend('6.9','10.7','18.7','23.8','36.5' );
    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
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
    sal=MY(1,7);
    di_snow=MY(2,5);

    if partype==1    
      di_snow=linspace(0.03,0.2,100);
      x_axis = di_snow;

      for n=1:100
        MY(2,5)=di_snow(n)*100;
        Tb_MY_memls=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', MY);
        Tb_memls_V(n,:) = Tb_MY_memls(:,1);
        Tb_memls_H(n,:) = Tb_MY_memls(:,2); 
        delta_snow=(0.2-0.03)/1000*100; %[cm]
        MY(2,5)=di_snow(n)*100+delta_snow;
        Tb_MY_memls_delta=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', MY);
        Tb_memls_V_delta(n,:) = Tb_MY_memls_delta(:,1);
        Tb_memls_H_delta(n,:) = Tb_MY_memls_delta(:,2); 
        Dir_V(n,:)= (Tb_memls_V_delta(n,:)-Tb_memls_V(n,:))/delta_snow;%k/cm
        Dir_H(n,:)= (Tb_memls_H_delta(n,:)-Tb_memls_H(n,:))/delta_snow; ...
            pol_delta(n,:)=Tb_memls_V_delta(n,:)-Tb_memls_H_delta(n,:);
        pol(n,:)=Tb_memls_V(n,:)-Tb_memls_H(n,:);
        Dir_pol(n,:)= (pol_delta(n,:)-pol(n,:))/delta_snow; 
        %k/cm
        par=' snow depth';
        unit = 'K/cm';
        xaxis=' snow depth [m]';

      end
      
      cd ../inversion
    end

    if partype==2 % snow density
      roi_snow=linspace(200,450,100);
      x_axis = roi_snow;
      for n=1:100
        MY(2,4)=roi_snow(n);
        
        Tb_MY_memls=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', MY);
        Tb_memls_V(n,:) = Tb_MY_memls(:,1);
        Tb_memls_H(n,:) = Tb_MY_memls(:,2); 
        delta_roi_snow=(450-200)/1000; 
        MY(2,4)=roi_snow(n)+delta_roi_snow;
        Tb_MY_memls_delta=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', MY);
        Tb_memls_V_delta(n,:) = Tb_MY_memls_delta(:,1);
        Tb_memls_H_delta(n,:) = Tb_MY_memls_delta(:,2); 
        Dir_V(n,:)= (Tb_memls_V_delta(n,:)-Tb_memls_V(n,:))/delta_roi_snow;%k/cm
        Dir_H(n,:)= (Tb_memls_H_delta(n,:)-Tb_memls_H(n,:))/delta_roi_snow; ...
        %k/cm
        pol_delta(n,:)=Tb_memls_V_delta(n,:)-Tb_memls_H_delta(n,:);
        pol(n,:)=Tb_memls_V(n,:)-Tb_memls_H(n,:);
        Dir_pol(n,:)= (pol_delta(n,:)-pol(n,:))/delta_roi_snow; 
        par=' snow density';
        unit = 'K/kg/m^3';
        xaxis=' snow density [kg/m^3]';
      end
      cd ../inversion
    end

    if partype==3 % grain size
      pci_snow=linspace(0.001,0.3,100);
      x_axis = pci_snow;
      for n=1:100
        MY(2,6)=pci_snow(n);
        Tb_MY_memls=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', MY);
        Tb_memls_V(n,:) = Tb_MY_memls(:,1);
        Tb_memls_H(n,:) = Tb_MY_memls(:,2); 
        
        pol(n,:)=Tb_memls_V(n,:)-Tb_memls_H(n,:);
        delta_pci_snow=(0.3-0.001)/2500; 
        MY(2,6)=pci_snow(n)+delta_pci_snow;
        Tb_MY_memls_delta=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', MY);
        Tb_memls_V_delta(n,:) = Tb_MY_memls_delta(:,1);
        Tb_memls_H_delta(n,:) = Tb_MY_memls_delta(:,2); 
        pol_delta(n,:)=Tb_memls_V_delta(n,:)-Tb_memls_H_delta(n,:);
        Dir_V(n,:)= (Tb_memls_V_delta(n,:)-Tb_memls_V(n,:))/delta_pci_snow;%k/cm
        Dir_H(n,:)= (Tb_memls_H_delta(n,:)-Tb_memls_H(n,:))/delta_pci_snow; ...
        %k/cm
        Dir_pol(n,:)= (pol_delta(n,:)-pol(n,:))/delta_pci_snow; 
        par=' grain size of snow';
        unit = 'K/mm';
        xaxis=' grain size of snow [mm]';
      end
      cd ../inversion
    end

    if partype==4 % salinity
      sal=linspace(0.05,8,100);
      x_axis =sal
      Tb_memls_V=zeros(100,5);
      Tb_memls_H=zeros(100,5);
      for n=1:100
        MY(1,7)=sal(n);
        Tb_MY_memls=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', MY);
        delta_sal=(8-0.05)/1000; %[cm]
        Tb_memls_V(n,:) = Tb_MY_memls(:,1);
        Tb_memls_H(n,:) = Tb_MY_memls(:,2); 
        pol(n,:)=Tb_memls_V(n,:)-Tb_memls_H(n,:);

        MY(1,7)=sal(n)+delta_sal;
        Tb_MY_memls_delta=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', MY);
        Tb_memls_V_delta(n,:) = Tb_MY_memls_delta(:,1);
        Tb_memls_H_delta(n,:) = Tb_MY_memls_delta(:,2);
        pol_delta(n,:)=Tb_memls_V_delta(n,:)-Tb_memls_H_delta(n,:);
        Dir_V(n,:)= (Tb_memls_V_delta(n,:)-Tb_memls_V(n,:))/delta_sal;%k/cm
        Dir_H(n,:)= (Tb_memls_H_delta(n,:)-Tb_memls_H(n,:))/delta_sal; ...
        %k/cm
        Dir_pol(n,:)= (pol_delta(n,:)-pol(n,:))/delta_sal; ...
        par = ' salinity';
        unit =  '?';
        xaxis = 'salinity ';
      end
      cd ../inversion
    end

        if partype==5 % T_s
      T_snow=linspace(235,270,100);
      x_axis =T_snow
      Tb_memls_V=zeros(100,5);
      Tb_memls_H=zeros(100,5);
      for n=1:100
        MY(2,2)=T_snow(n);
        Tb_MY_memls=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', MY);
        Tb_memls_V(n,:) = Tb_MY_memls(:,1);
        Tb_memls_H(n,:) = Tb_MY_memls(:,2); 
        delta_T_snow=(270-235)/1000*100;
        MY(2,2)=T_snow(n)+delta_T_snow;
        Tb_MY_memls_delta=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', MY);
        Tb_memls_V_delta(n,:) = Tb_MY_memls_delta(:,1);
        Tb_memls_H_delta(n,:) = Tb_MY_memls_delta(:,2); 
        Dir_V(n,:)= (Tb_memls_V_delta(n,:)-Tb_memls_V(n,:))/delta_T_snow;%k/cm
        Dir_H(n,:)= (Tb_memls_H_delta(n,:)-Tb_memls_H(n,:))/delta_T_snow; ...
        %k/cm
        pol_delta(n,:)=Tb_memls_V_delta(n,:)-Tb_memls_H_delta(n,:);
        pol(n,:)=Tb_memls_V(n,:)-Tb_memls_H(n,:);
        Dir_pol(n,:)= (pol_delta(n,:)-pol(n,:))/delta_T_snow; 
        par = ' snow temperature';
        unit =  ' [K/K]';
        xaxis = '  snow temperature[K]';
      end
      cd ../inversion
    end

    figure
    subplot(1,2,1)
    grid on
    hold on
    
    title (['Sensitivity of Tb_V according to ',par, ' (MY)']);
    plot(x_axis,Dir_V(:,1),'r--');
    plot(x_axis,Dir_V(:,2),'g--');
    plot(x_axis,Dir_V(:,3),'b--');
    plot(x_axis,Dir_V(:,4),'k--');
    plot(x_axis,Dir_V(:,5),'m--');
    hold off
    ylabel(['sensitivty of Tb_v ',unit]);
    xlabel(xaxis);
    legend('6.9','10.7','18.7','23.8','36.5' );
    
    subplot(1,2,2)
    grid on
    hold on

    title (['Sensitivity of Tb_H according to ',par, ' (MY)']);
    plot(x_axis,Dir_H(:,1),'r--');
    plot(x_axis,Dir_H(:,2),'g--');
    plot(x_axis,Dir_H(:,3),'b--');
    plot(x_axis,Dir_H(:,4),'k--');
    plot(x_axis,Dir_H(:,5),'m--');
    hold off
    ylabel(['sensitivty of Tb_h',unit]);
    xlabel(xaxis);
    legend('6.9','10.7','18.7','23.8','36.5' );
    
    
    figure
    grid on
    hold on
    title (['Sensitivity of polarization according to ',par, ' (MY)']);
    plot(x_axis,Dir_pol(:,1),'r--');
    plot(x_axis,Dir_pol(:,2),'g--');
    plot(x_axis,Dir_pol(:,3),'b--');
    plot(x_axis,Dir_pol(:,4),'k--');
    plot(x_axis,Dir_pol(:,5),'m--');
    hold off
    ylabel(['sensitivty of polarization ',unit]);
    xlabel(xaxis);
    legend('6.9','10.7','18.7','23.8','36.5' );
  end
  
