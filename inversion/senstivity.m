
function[] = sens(type)

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

    di_snow=linspace(0.03,0.2,100);
    Tb_memls_V=zeros(100,5);
    Tb_memls_H=zeros(100,5);
    for n=1:100
      FY(2,5)=di_snow(n)*100;
      Tb_FY_memls=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', FY);
      Tb_memls_V(n,:) = Tb_FY_memls(:,1);
      Tb_memls_H(n,:) = Tb_FY_memls(:,2); 
      delta_snow=(0.2-0.03)/1000*100; %[cm]
      FY(2,5)=di_snow(n)*100+delta_snow;
      Tb_FY_memls_delta=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', FY);
      Tb_memls_V_dleta(n,:) = Tb_FY_memls_delta(:,1);
      Tb_memls_H_delta(n,:) = Tb_FY_memls_delta(:,2); 
      Dir_V_depth(n,:)= (Tb_memls_V_delta(n,:)-Tb_memls_V(n,:))/delta_snow;%k/cm
      Dir_H_depth(n,:)= (Tb_memls_H_delta(n,:)-Tb_memls_H(n,:))/delta_snow;%k/cm
    end
    
    figure
    grid on
    hold on
    title ('Sensitivity of Tb_V according to snow depth ');
    plot(di_snow,Dir_V_depth(:,1),'r--');
    plot(di_snow,Dir_V_depth(:,2),'g--');
    plot(di_snow,Dir_V_depth(:,3),'b--');
    plot(di_snow,Dir_V_depth(:,4),'k--');
    plot(di_snow,Dir_V_depth(:,5),'m--');
    hold off
    ylabel('sensitivty of Tb_v [K/cm]');
    xlabel('snow depth [m]');
    cd ../inversion
    
  end
  
